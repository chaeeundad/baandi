//! Default-off, bounded collective arrival capture. No MPI or I/O here.
//! Epoch anchors are bracketed by SystemTime on ONE host; monotonic elapsed
//! times supply event stamps. Post-run drift is reported, never hidden.
use std::cell::{Cell, RefCell};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

pub const WORDS: usize = 18;
pub const CAPACITY: usize = 131_072;
thread_local! { static ENABLED: Cell<bool> = const { Cell::new(false) }; }
thread_local! { static TRACE: RefCell<Option<State>> = const { RefCell::new(None) }; }
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Trace {
    pub epoch_ns: u64,
    pub error_ns: u64,
    pub drift_ns: u64,
    pub timed: bool,
    pub overflow: bool,
    pub events: Vec<[u64; WORDS]>,
}
struct State {
    trace: Trace,
    anchor: Option<Instant>,
    site: u64,
    map: u64,
    active: bool,
    map_start: usize,
    seq: u64,
    previous_exit: u64,
    stage: [u64; 8], // ID, duration, points, bytes, finish
}
fn epoch() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos() as u64
}
#[inline]
pub fn enabled() -> bool {
    ENABLED.with(Cell::get)
}
/// Called at profile initialization, before the first map.
pub fn initialize() {
    let config = crate::switches::RunConfig::current();
    let on = config.on(crate::switches::Id::COLLECTIVE_TRACE);
    ENABLED.with(|v| v.set(on));
    if !on {
        return;
    }
    TRACE.with(|s| {
        if s.borrow().is_some() {
            return;
        }
        let timed = config.on(crate::switches::Id::STAGE_TIMERS);
        let (anchor, epoch_ns, error_ns) = if timed {
            let a = epoch();
            let i = Instant::now();
            let b = epoch();
            (
                Some(i),
                a + b.saturating_sub(a) / 2,
                b.saturating_sub(a) / 2 + 1,
            )
        } else {
            (None, 0, 0)
        };
        *s.borrow_mut() = Some(State {
            trace: Trace {
                epoch_ns,
                error_ns,
                timed,
                events: Vec::with_capacity(CAPACITY),
                ..Trace::default()
            },
            anchor,
            site: 0,
            map: 0,
            active: false,
            map_start: 0,
            seq: 0,
            previous_exit: 0,
            stage: [0; 8],
        });
    });
}
fn with(f: impl FnOnce(&mut State)) {
    if enabled() {
        TRACE.with(|s| {
            if let Some(s) = s.borrow_mut().as_mut() {
                f(s);
            }
        });
    }
}
impl State {
    fn now(&self) -> u64 {
        self.anchor.map_or(0, |i| i.elapsed().as_nanos() as u64)
    }
}
pub fn map_begin() {
    with(|s| {
        s.map += 1;
        s.active = true;
        s.map_start = s.trace.events.len();
        s.seq = 0;
        s.previous_exit = s.now();
        s.stage = [0; 8];
    });
}
pub fn map_commit(bucket: usize) {
    with(|s| {
        for e in &mut s.trace.events[s.map_start..] {
            e[1] = bucket as u64;
        }
        s.active = false;
        s.map_start = s.trace.events.len();
    });
}
/// Most recently completed stage, not an exclusive decomposition. IDs 1..64
/// are stage_clock slots; 101..118 FFT; 201..208 M2 operator; 301 projector GEMM.
pub fn stage(id: usize, ns: u64, points: usize, bytes: usize) {
    with(|s| {
        if s.active {
            s.stage = [id as u64, ns, points as u64, bytes as u64, s.now(), 0, 0, 0];
        }
    });
}
/// Mark one logical collective; off does not access TLS or read a clock.
pub fn begin() -> Option<u64> {
    if !enabled() {
        return None;
    }
    TRACE.with(|s| s.borrow().as_ref().filter(|s| s.active).map(State::now))
}
pub fn end(mark: Option<u64>, op: usize, bytes: usize) {
    let Some(enter) = mark else {
        return;
    };
    with(|s| {
        let exit = s.now();
        let [stage, ns, points, stage_bytes, finish, a, b, c] = s.stage;
        if s.trace.events.len() < CAPACITY {
            s.trace.events.push([
                s.map,
                0,
                s.seq,
                op as u64,
                bytes as u64,
                enter,
                exit,
                stage,
                ns,
                points,
                stage_bytes,
                finish,
                s.previous_exit,
                s.trace.events.len() as u64,
                a,
                b,
                c,
                s.site,
            ]);
        } else {
            s.trace.overflow = true;
        }
        s.site = 0;
        s.seq += 1;
        s.previous_exit = exit;
        // Never attribute the next interval to a stale pre-collective stage.
        s.stage = [0; 8];
    });
}
pub fn take() -> Trace {
    if !enabled() {
        return Trace::default();
    }
    TRACE.with(|s| {
        s.borrow_mut().take().map_or_else(Trace::default, |mut s| {
            if s.trace.timed {
                let a = epoch();
                let predicted = s.trace.epoch_ns + s.now();
                let b = epoch();
                s.trace.drift_ns = if predicted < a {
                    a - predicted
                } else {
                    predicted.saturating_sub(b)
                };
                s.trace.error_ns += b.saturating_sub(a);
            }
            s.trace
        })
    })
}
/// Lossless u64 wire representation; used only AFTER the SCF.
pub fn encode(t: &Trace) -> Vec<f64> {
    [
        t.epoch_ns,
        t.error_ns,
        t.drift_ns,
        t.timed as u64,
        t.overflow as u64,
    ]
    .into_iter()
    .chain(t.events.iter().flatten().copied())
    .flat_map(|v| [v as u32 as f64, (v >> 32) as f64])
    .collect()
}
pub fn decode(w: &[f64]) -> Trace {
    let v: Vec<u64> = w
        .as_chunks::<2>()
        .0
        .iter()
        .map(|p| p[0] as u64 | ((p[1] as u64) << 32))
        .collect();
    Trace {
        epoch_ns: v[0],
        error_ns: v[1],
        drift_ns: v[2],
        timed: v[3] != 0,
        overflow: v[4] != 0,
        events: v[5..].as_chunks::<WORDS>().0.to_vec(),
    }
}

/// Clock for a stage without existing timing: disabled trace reads no clock.
pub fn clock() -> Option<Instant> {
    if !enabled() {
        return None;
    }
    TRACE.with(|s| {
        s.borrow()
            .as_ref()
            .filter(|s| s.active && s.trace.timed)
            .map(|_| Instant::now())
    })
}
pub fn shape(a: usize, b: usize, c: usize) {
    with(|s| {
        s.stage[5..].copy_from_slice(&[a as u64, b as u64, c as u64]);
    });
}

/// Optional semantic site of the next collective; reset on collective exit.
pub fn site(id: usize) {
    with(|s| {
        if s.active {
            s.site = id as u64;
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::switches::RunConfig;
    #[test]
    fn arrival_bounded_count_only_and_disabled() {
        let off = RunConfig::from_pairs(&[("collective_trace", "off")]).unwrap();
        let guard = off.enter();
        initialize();
        assert!(begin().is_none());
        assert!(TRACE.with(|s| s.borrow().is_none()));
        drop(guard);
        let on =
            RunConfig::from_pairs(&[("collective_trace", "on"), ("stage_timers", "off")]).unwrap();
        let _guard = on.enter();
        initialize();
        map_begin();
        let pointer = TRACE.with(|s| s.borrow().as_ref().unwrap().trace.events.as_ptr());
        stage(101, 0, 12, 384);
        let mark = begin();
        end(mark, 3, 96);
        assert_eq!(
            pointer,
            TRACE.with(|s| s.borrow().as_ref().unwrap().trace.events.as_ptr())
        );
        TRACE.with(|s| {
            s.borrow_mut()
                .as_mut()
                .unwrap()
                .trace
                .events
                .resize(CAPACITY, [0; WORDS])
        });
        end(begin(), 3, 96);
        map_commit(2);
        let t = take();
        assert!(t.overflow);
        assert!(!t.timed);
        assert_eq!(t.events.len(), CAPACITY);
        assert_eq!(t.events[0][7], 101);
        assert_eq!(t.events[0][5], 0);
        assert_eq!(t.events[0][6], 0);
    }
}
