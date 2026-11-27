//! Full-screen monitor (ratatui + crossterm). Read-only: it tails the stream
//! and never talks to the job. Keys: q/Esc quit, Tab/Shift-Tab switch panel.
use super::{JobStatus, State, Tail, job_status, sci};
use ratatui::{
    Frame,
    crossterm::event::{self, Event, KeyCode, KeyEventKind},
    layout::{Constraint, Layout, Rect},
    style::{Color, Modifier, Style},
    symbols,
    text::{Line, Span},
    widgets::{Axis, Block, Cell, Chart, Dataset, GraphType, Paragraph, Row, Table, Tabs},
};
use std::{
    path::Path,
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};

/// Name, points and colour of one chart line.
type Series = (String, Vec<(f64, f64)>, Color);

const TABS: [&str; 4] = ["Overview", "SCF maps", "Relaxation", "Events"];

pub(super) fn run(dir: &Path, mut tail: Tail) -> std::io::Result<()> {
    let mut terminal = ratatui::init();
    let mut tab = 0usize;
    let mut status = job_status(dir, &tail.state);
    let mut last_poll = Instant::now();
    let result = loop {
        if let Err(e) = terminal.draw(|f| draw(f, dir, &tail.state, &status, tab)) {
            break Err(e);
        }
        match event::poll(Duration::from_millis(250)) {
            Ok(true) => match event::read() {
                Ok(Event::Key(k)) if k.kind == KeyEventKind::Press => match k.code {
                    KeyCode::Char('q') | KeyCode::Esc => break Ok(()),
                    KeyCode::Tab | KeyCode::Right => tab = (tab + 1) % TABS.len(),
                    KeyCode::BackTab | KeyCode::Left => tab = (tab + TABS.len() - 1) % TABS.len(),
                    KeyCode::Char(c @ '1'..='4') => tab = c as usize - '1' as usize,
                    _ => {}
                },
                Ok(_) => {}
                Err(e) => break Err(e),
            },
            Ok(false) => {}
            Err(e) => break Err(e),
        }
        if last_poll.elapsed() >= Duration::from_millis(500) {
            tail.poll();
            status = job_status(dir, &tail.state);
            last_poll = Instant::now();
        }
    };
    ratatui::restore();
    result
}

fn status_style(status: &JobStatus) -> Style {
    match status {
        JobStatus::Running { .. } | JobStatus::Remote { .. } => Style::new().black().on_green(),
        JobStatus::Finished { status, .. } if status == "completed" => {
            Style::new().black().on_cyan()
        }
        _ => Style::new().white().on_red(),
    }
}

fn draw(f: &mut Frame, dir: &Path, state: &State, status: &JobStatus, tab: usize) {
    let [header, tabs, body, footer] = Layout::vertical([
        Constraint::Length(4),
        Constraint::Length(1),
        Constraint::Min(8),
        Constraint::Length(1),
    ])
    .areas(f.area());
    let elapsed = match (state.start_unix, status.live()) {
        (Some(start), true) => {
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_or(0., |d| d.as_secs_f64())
                - start
        }
        _ => state.elapsed_s,
    };
    let title = Line::from(vec![
        Span::styled(" BAANDI ", Style::new().bold().black().on_yellow()),
        Span::raw(format!(" {} ", state.task())),
        Span::styled(format!(" {} ", status.label()), status_style(status)),
        Span::raw(format!("  {}  elapsed {}", dir.display(), clock(elapsed))),
    ]);
    let head = Paragraph::new(vec![
        Line::from(state.system_line()),
        Line::from(vec![
            Span::styled("Phase: ", Style::new().bold()),
            Span::raw(state.phase_line()),
        ]),
    ])
    .block(Block::bordered().title(title));
    f.render_widget(head, header);
    f.render_widget(
        Tabs::new(TABS)
            .select(tab)
            .highlight_style(Style::new().bold().reversed())
            .divider(symbols::line::VERTICAL),
        tabs,
    );
    match tab {
        0 => overview(f, state, body),
        1 => scf_table(f, state, body, usize::MAX),
        2 => relax_panel(f, state, body),
        _ => events(f, state, body, usize::MAX),
    }
    f.render_widget(
        Paragraph::new(
            " q quit   Tab/Shift-Tab or 1-4 switch panel   read-only view of progress.jsonl",
        )
        .style(Style::new().dark_gray()),
        footer,
    );
}

fn clock(seconds: f64) -> String {
    let s = seconds.max(0.) as u64;
    format!("{:02}:{:02}:{:02}", s / 3600, (s / 60) % 60, s % 60)
}

fn overview(f: &mut Frame, state: &State, area: Rect) {
    let [top, bottom] =
        Layout::vertical([Constraint::Percentage(55), Constraint::Percentage(45)]).areas(area);
    let [left, right] =
        Layout::horizontal([Constraint::Percentage(55), Constraint::Percentage(45)]).areas(top);
    scf_chart(f, state, left);
    if state.steps.is_empty() && state.relax_settings().is_none() {
        energy_chart(f, state, right);
    } else {
        relax_chart(f, state, right);
    }
    let [l, r] =
        Layout::horizontal([Constraint::Percentage(55), Constraint::Percentage(45)]).areas(bottom);
    scf_table(f, state, l, (l.height as usize).saturating_sub(3));
    events(f, state, r, (r.height as usize).saturating_sub(2));
}

fn bounds(values: impl Iterator<Item = f64>) -> [f64; 2] {
    let (lo, hi) = values
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| {
            (a.min(v), b.max(v))
        });
    if !lo.is_finite() {
        [0., 1.]
    } else if hi - lo < 1e-12 {
        [lo - 0.5, hi + 0.5]
    } else {
        let pad = 0.05 * (hi - lo);
        [lo - pad, hi + pad]
    }
}

/// Three axis labels with enough decimals to tell them apart (`digits` is
/// the minimum).
fn labels(b: [f64; 2], digits: usize) -> Vec<Line<'static>> {
    let mid = 0.5 * (b[0] + b[1]);
    let span = (b[1] - b[0]).abs();
    let digits = if span > 0. {
        digits.max((1. - (span / 2.).log10()).ceil().clamp(0., 10.) as usize)
    } else {
        digits
    };
    [b[0], mid, b[1]]
        .iter()
        .map(|v| Line::from(format!("{v:.digits$}")))
        .collect()
}

fn line_chart(
    f: &mut Frame,
    area: Rect,
    title: String,
    x_title: &str,
    y_title: &str,
    series: Vec<Series>,
    threshold: Option<f64>,
) {
    // x (map or step number) spans the data exactly; y is padded by 5 %.
    let (x0, x1) = series
        .iter()
        .flat_map(|s| s.1.iter().map(|p| p.0))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| {
            (a.min(v), b.max(v))
        });
    let xb = if x0.is_finite() && x1 > x0 {
        [x0, x1]
    } else {
        [0., 1.]
    };
    let x_digits = usize::from(xb[1] - xb[0] < 2.);
    let yb = bounds(
        series
            .iter()
            .flat_map(|s| s.1.iter().map(|p| p.1))
            .chain(threshold),
    );
    let limit: Vec<(f64, f64)> = threshold
        .map(|t| vec![(xb[0], t), (xb[1], t)])
        .unwrap_or_default();
    let mut datasets: Vec<Dataset> = series
        .iter()
        .map(|(name, data, color)| {
            Dataset::default()
                .name(name.clone())
                .marker(symbols::Marker::Braille)
                .graph_type(GraphType::Line)
                .style(Style::new().fg(*color))
                .data(data)
        })
        .collect();
    if threshold.is_some() {
        datasets.push(
            Dataset::default()
                .name("threshold")
                .marker(symbols::Marker::Dot)
                .graph_type(GraphType::Line)
                .style(Style::new().red())
                .data(&limit),
        );
    }
    let chart = Chart::new(datasets)
        .block(Block::bordered().title(title))
        .x_axis(
            Axis::default()
                .title(x_title.to_owned())
                .bounds(xb)
                .labels(labels(xb, x_digits)),
        )
        .y_axis(
            Axis::default()
                .title(y_title.to_owned())
                .bounds(yb)
                .labels(labels(yb, 0)),
        );
    f.render_widget(chart, area);
}

fn scf_chart(f: &mut Frame, state: &State, area: Rect) {
    let Some(c) = state.cycles.last() else {
        f.render_widget(
            Paragraph::new("waiting for the first SCF map")
                .block(Block::bordered().title("SCF convergence")),
            area,
        );
        return;
    };
    // x is the record number in this SCF (retry and verification maps
    // repeat an iteration number).
    let dr2: Vec<(f64, f64)> = (1..)
        .zip(&c.maps)
        .filter(|(_, m)| m.dr2_ry > 0.)
        .map(|(i, m)| (f64::from(i), m.dr2_ry.log10()))
        .collect();
    let de: Vec<(f64, f64)> = (1..)
        .zip(&c.maps)
        .filter_map(|(i, m)| {
            m.de_ha
                .filter(|d| *d != 0.)
                .map(|d| (f64::from(i), d.abs().log10()))
        })
        .collect();
    line_chart(
        f,
        area,
        format!("SCF {} convergence (log10)", c.index),
        "map #",
        "log10",
        vec![
            ("dr2 [Ry]".into(), dr2, Color::Cyan),
            ("|dE| [Ha]".into(), de, Color::Yellow),
        ],
        (c.threshold_ry > 0.).then(|| c.threshold_ry.log10()),
    );
}

fn energy_chart(f: &mut Frame, state: &State, area: Rect) {
    let Some(c) = state.cycles.last() else {
        f.render_widget(Block::bordered().title("Energy"), area);
        return;
    };
    let e: Vec<(f64, f64)> = (1..)
        .zip(&c.maps)
        .filter_map(|(i, m)| m.energy_ha.map(|e| (f64::from(i), e)))
        .collect();
    line_chart(
        f,
        area,
        "Total energy [Ha]".into(),
        "map #",
        "Ha",
        vec![("E".into(), e, Color::Green)],
        None,
    )
}

fn relax_chart(f: &mut Frame, state: &State, area: Rect) {
    let force: Vec<(f64, f64)> = state
        .steps
        .iter()
        .filter(|s| s.max_force_ha_bohr > 0.)
        .map(|s| (s.step as f64, s.max_force_ha_bohr.log10()))
        .collect();
    let tol = state
        .relax_settings()
        .and_then(|r| r["force_tol_ha_bohr"].as_f64())
        .filter(|t| *t > 0.)
        .map(f64::log10);
    let [a, b] =
        Layout::vertical([Constraint::Percentage(50), Constraint::Percentage(50)]).areas(area);
    line_chart(
        f,
        a,
        "Relaxation: log10 max|F| [Ha/Bohr]".into(),
        "step",
        "log10",
        vec![("max|F|".into(), force, Color::Magenta)],
        tol,
    );
    let energy: Vec<(f64, f64)> = state
        .steps
        .iter()
        .map(|s| (s.step as f64, s.energy_ha))
        .collect();
    line_chart(
        f,
        b,
        "Relaxation: energy [Ha]".into(),
        "step",
        "Ha",
        vec![("E".into(), energy, Color::Green)],
        None,
    );
}

fn scf_table(f: &mut Frame, state: &State, area: Rect, rows: usize) {
    let maps = state.cycles.last().map_or(&[][..], |c| &c.maps[..]);
    let skip = maps.len().saturating_sub(rows);
    let body: Vec<Row> = maps[skip..]
        .iter()
        .map(|m| {
            let style = match m.kind.as_str() {
                "verification" => Style::new().cyan(),
                "retry" => Style::new().red(),
                _ => Style::new(),
            };
            Row::new(vec![
                Cell::from(m.iteration.to_string()),
                Cell::from(m.kind.clone()),
                Cell::from(
                    m.energy_ha
                        .map_or_else(|| "-".into(), |e| format!("{e:.10}")),
                ),
                Cell::from(sci(m.de_ha)),
                Cell::from(sci(Some(m.dr2_ry))),
                Cell::from(format!("{:.3}", m.seconds)),
            ])
            .style(style)
        })
        .collect();
    let title = state.cycles.last().map_or_else(
        || "SCF maps".to_owned(),
        |c| {
            format!(
                "SCF {} of {} maps (stop dr2 < {} Ry)",
                c.index,
                state.cycles.len(),
                sci(Some(c.threshold_ry))
            )
        },
    );
    let table = Table::new(
        body,
        [
            Constraint::Length(5),
            Constraint::Length(13),
            Constraint::Length(18),
            Constraint::Length(10),
            Constraint::Length(10),
            Constraint::Length(8),
        ],
    )
    .header(
        Row::new([
            "map",
            "kind",
            "energy (Ha)",
            "dE (Ha)",
            "dr2 (Ry)",
            "time s",
        ])
        .style(Style::new().add_modifier(Modifier::BOLD)),
    )
    .block(Block::bordered().title(title));
    f.render_widget(table, area);
}

fn relax_panel(f: &mut Frame, state: &State, area: Rect) {
    let [top, bottom] =
        Layout::vertical([Constraint::Percentage(50), Constraint::Percentage(50)]).areas(area);
    relax_chart(f, state, top);
    let skip = state
        .steps
        .len()
        .saturating_sub((bottom.height as usize).saturating_sub(3));
    let rows: Vec<Row> = state.steps[skip..]
        .iter()
        .map(|s| {
            Row::new(vec![
                s.step.to_string(),
                format!("{:.10}", s.energy_ha),
                sci(Some(s.max_force_ha_bohr)),
                sci(Some(s.step_length_bohr)),
                s.pressure_ha_bohr3
                    .map_or_else(|| "-".into(), |p| format!("{p:.2e}")),
                s.scf_iterations.to_string(),
                format!("{:.1}", s.seconds),
                s.provenance.clone(),
            ])
        })
        .collect();
    let table = Table::new(
        rows,
        [
            Constraint::Length(5),
            Constraint::Length(18),
            Constraint::Length(11),
            Constraint::Length(10),
            Constraint::Length(10),
            Constraint::Length(5),
            Constraint::Length(8),
            Constraint::Min(8),
        ],
    )
    .header(
        Row::new([
            "step",
            "energy (Ha)",
            "max|F|",
            "step Bohr",
            "P Ha/B^3",
            "maps",
            "time s",
            "provenance",
        ])
        .style(Style::new().bold()),
    )
    .block(Block::bordered().title("Relaxation steps"));
    f.render_widget(table, bottom);
}

fn events(f: &mut Frame, state: &State, area: Rect, rows: usize) {
    let skip = state
        .events
        .len()
        .saturating_sub(rows.min(state.events.len()));
    let lines: Vec<Line> = state.events[skip..]
        .iter()
        .map(|(t, text)| {
            let style =
                if text.starts_with("warning") || text.contains("failed") || text.contains("NOT") {
                    Style::new().yellow()
                } else if text.starts_with("end") {
                    Style::new().bold()
                } else {
                    Style::new()
                };
            Line::from(vec![
                Span::styled(format!("{:>8.1}s ", t), Style::new().dark_gray()),
                Span::styled(text.clone(), style),
            ])
        })
        .collect();
    f.render_widget(
        Paragraph::new(lines).block(Block::bordered().title("Events and warnings")),
        area,
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use ratatui::{Terminal, backend::TestBackend};

    /// Every panel renders the recorded stream at a large and a tiny size
    /// without panicking; the overview carries the key labels.
    #[test]
    fn panels_render_the_recorded_stream() {
        let state = State::parse(include_str!("fixtures/relax.jsonl"));
        let status = JobStatus::Finished {
            status: "not_converged".into(),
            exit_code: Some(3),
        };
        for (w, h) in [(120, 40), (30, 10)] {
            let mut terminal = Terminal::new(TestBackend::new(w, h)).unwrap();
            for tab in 0..TABS.len() {
                terminal
                    .draw(|f| draw(f, Path::new("results"), &state, &status, tab))
                    .unwrap();
            }
        }
        let mut terminal = Terminal::new(TestBackend::new(120, 40)).unwrap();
        terminal
            .draw(|f| draw(f, Path::new("results"), &state, &status, 0))
            .unwrap();
        let text: String = terminal
            .backend()
            .buffer()
            .content()
            .iter()
            .map(|c| c.symbol())
            .collect::<Vec<_>>()
            .chunks(120)
            .map(|row| row.concat() + "\n")
            .collect();
        if std::env::var_os("PRINT_TUI").is_some() {
            println!("{text}");
        }
        for label in [
            "BAANDI",
            "NOT_CONVERGED",
            "SCF 4 convergence",
            "max|F|",
            "Events and warnings",
        ] {
            assert!(text.contains(label), "missing {label}\n{text}");
        }
    }
}
