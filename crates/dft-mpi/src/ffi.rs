//! Open MPI ABI adapter. Handles are pointers; predefined communicator,
//! datatype and operation objects are the `ompi_mpi_*` data symbols.
use crate::MpiError;
use std::ffi::{c_char, c_int, c_void};

type Handle = *mut c_void;

#[repr(C)]
struct Opaque {
    _private: [u8; 0],
}

unsafe extern "C" {
    static ompi_mpi_comm_world: Opaque;
    static ompi_mpi_double: Opaque;
    static ompi_mpi_op_sum: Opaque;
    fn MPI_Init(argc: *mut c_int, argv: *mut *mut *mut c_char) -> c_int;
    fn MPI_Initialized(flag: *mut c_int) -> c_int;
    fn MPI_Finalized(flag: *mut c_int) -> c_int;
    fn MPI_Finalize() -> c_int;
    fn MPI_Comm_rank(comm: Handle, rank: *mut c_int) -> c_int;
    fn MPI_Comm_size(comm: Handle, size: *mut c_int) -> c_int;
    fn MPI_Barrier(comm: Handle) -> c_int;
    fn MPI_Wtime() -> f64;
    fn MPI_Bcast(
        buffer: *mut c_void,
        count: c_int,
        datatype: Handle,
        root: c_int,
        comm: Handle,
    ) -> c_int;
    fn MPI_Allgatherv(
        sendbuf: *const c_void,
        sendcount: c_int,
        sendtype: Handle,
        recvbuf: *mut c_void,
        recvcounts: *const c_int,
        displs: *const c_int,
        recvtype: Handle,
        comm: Handle,
    ) -> c_int;
    fn MPI_Send(
        buf: *const c_void,
        count: c_int,
        datatype: Handle,
        dest: c_int,
        tag: c_int,
        comm: Handle,
    ) -> c_int;
    fn MPI_Recv(
        buf: *mut c_void,
        count: c_int,
        datatype: Handle,
        source: c_int,
        tag: c_int,
        comm: Handle,
        status: *mut c_void,
    ) -> c_int;
    fn MPI_Alltoall(
        sendbuf: *const c_void,
        sendcount: c_int,
        sendtype: Handle,
        recvbuf: *mut c_void,
        recvcount: c_int,
        recvtype: Handle,
        comm: Handle,
    ) -> c_int;
}

fn check(name: &'static str, code: c_int) -> Result<(), MpiError> {
    if code == 0 {
        Ok(())
    } else {
        Err(MpiError::Call(name, code))
    }
}
fn count(n: usize, what: &'static str) -> Result<c_int, MpiError> {
    c_int::try_from(n).map_err(|_| MpiError::Shape(what))
}

/// Initialises MPI (single-threaded semantics: the library is called from the
/// thread that created the Universe) and finalises it on drop.
pub struct Universe {
    owned: bool,
}
impl Universe {
    pub fn init() -> Result<Self, MpiError> {
        let (mut up, mut down): (c_int, c_int) = (0, 0);
        // SAFETY: both calls only write the flag cells passed by reference.
        unsafe {
            check("MPI_Initialized", MPI_Initialized(&mut up))?;
            check("MPI_Finalized", MPI_Finalized(&mut down))?;
        }
        if down != 0 {
            return Err(MpiError::State("MPI already finalised"));
        }
        if up != 0 {
            return Ok(Self { owned: false });
        }
        // SAFETY: MPI_Init accepts null argc/argv (MPI-2 and later).
        unsafe {
            check(
                "MPI_Init",
                MPI_Init(std::ptr::null_mut(), std::ptr::null_mut()),
            )?;
        }
        Ok(Self { owned: true })
    }
    pub fn world(&self) -> World {
        World { _private: () }
    }
}
impl Drop for Universe {
    fn drop(&mut self) {
        if self.owned {
            // SAFETY: finalise exactly once, after init by this Universe.
            let _ = unsafe { MPI_Finalize() };
        }
    }
}

/// MPI_COMM_WORLD.
#[derive(Clone, Copy)]
pub struct World {
    _private: (),
}
impl World {
    fn comm(&self) -> Handle {
        // Only the address of the library's predefined object is taken.
        std::ptr::addr_of!(ompi_mpi_comm_world) as Handle
    }
    fn double(&self) -> Handle {
        std::ptr::addr_of!(ompi_mpi_double) as Handle
    }
    #[allow(dead_code)]
    fn op_sum(&self) -> Handle {
        std::ptr::addr_of!(ompi_mpi_op_sum) as Handle
    }
    pub fn rank(&self) -> usize {
        let mut r: c_int = 0;
        // SAFETY: valid communicator, result cell passed by reference.
        let _ = unsafe { MPI_Comm_rank(self.comm(), &mut r) };
        r.max(0) as usize
    }
    pub fn size(&self) -> usize {
        let mut s: c_int = 1;
        // SAFETY: valid communicator, result cell passed by reference.
        let _ = unsafe { MPI_Comm_size(self.comm(), &mut s) };
        s.max(1) as usize
    }
    pub fn barrier(&self) -> Result<(), MpiError> {
        // SAFETY: valid communicator.
        check("MPI_Barrier", unsafe { MPI_Barrier(self.comm()) })
    }
    pub fn wtime(&self) -> f64 {
        // SAFETY: no arguments.
        unsafe { MPI_Wtime() }
    }
    pub fn bcast_f64(&self, data: &mut [f64], root: usize) -> Result<(), MpiError> {
        if root >= self.size() {
            return Err(MpiError::Shape("root out of range"));
        }
        // SAFETY: the buffer holds `data.len()` doubles and is uniquely borrowed.
        check("MPI_Bcast", unsafe {
            MPI_Bcast(
                data.as_mut_ptr().cast(),
                count(data.len(), "bcast count")?,
                self.double(),
                root as c_int,
                self.comm(),
            )
        })
    }
    pub fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), MpiError> {
        let size = self.size();
        if counts.len() != size
            || counts[self.rank()] != local.len()
            || counts.iter().sum::<usize>() != out.len()
        {
            return Err(MpiError::Shape("allgatherv counts"));
        }
        let rc: Vec<c_int> = counts
            .iter()
            .map(|&c| count(c, "allgatherv count"))
            .collect::<Result<_, _>>()?;
        let mut at = 0;
        let displs: Vec<c_int> = counts
            .iter()
            .map(|&c| {
                let d = at;
                at += c;
                count(d, "allgatherv displacement")
            })
            .collect::<Result<_, _>>()?;
        // SAFETY: counts/displacements describe disjoint slices of `out`, whose
        // length equals their sum; `local` holds counts[rank] doubles; all
        // buffers are borrowed for the call.
        check("MPI_Allgatherv", unsafe {
            MPI_Allgatherv(
                local.as_ptr().cast(),
                count(local.len(), "allgatherv send count")?,
                self.double(),
                out.as_mut_ptr().cast(),
                rc.as_ptr(),
                displs.as_ptr(),
                self.double(),
                self.comm(),
            )
        })
    }
    /// Elementwise sum over ranks added in rank order after an allgather
    /// (bitwise independent of the library's reduction tree).
    pub fn ordered_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), MpiError> {
        if out.len() != local.len() {
            return Err(MpiError::Shape("ordered_sum lengths"));
        }
        let size = self.size();
        let n = local.len();
        let counts = vec![n; size];
        let mut all = vec![0.0; n * size];
        self.allgatherv_f64(local, &counts, &mut all)?;
        for (i, o) in out.iter_mut().enumerate() {
            let mut s = 0.0;
            for r in 0..size {
                s += all[r * n + i];
            }
            *o = s;
        }
        Ok(())
    }
    /// Fixed-order binomial-tree sum: at step s = 1, 2, 4, ... a rank with bit s set
    /// sends its partial sum to rank - s, which adds the incoming partial AFTER its
    /// own (acc = acc + incoming); rank 0 then broadcasts. For a fixed rank count
    /// the association is a fixed function of the rank, so every rank returns
    /// the same bits on every run, with 2 log2(P) small messages instead of the
    /// P-fold data volume of `ordered_sum_f64`.
    pub fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), MpiError> {
        if out.len() != local.len() {
            return Err(MpiError::Shape("tree_sum lengths"));
        }
        let (rank, size) = (self.rank(), self.size());
        let n = count(local.len(), "tree_sum count")?;
        out.copy_from_slice(local);
        let mut incoming = vec![0.0; local.len()];
        let mut step = 1;
        while step < size {
            if rank & step != 0 {
                // SAFETY: `out` holds n doubles, borrowed for the blocking call.
                check("MPI_Send", unsafe {
                    MPI_Send(
                        out.as_ptr().cast(),
                        n,
                        self.double(),
                        (rank - step) as c_int,
                        7001,
                        self.comm(),
                    )
                })?;
                break;
            } else if rank + step < size {
                // SAFETY: `incoming` holds n doubles; the null status pointer is
                // MPI_STATUS_IGNORE in the Open MPI ABI.
                check("MPI_Recv", unsafe {
                    MPI_Recv(
                        incoming.as_mut_ptr().cast(),
                        n,
                        self.double(),
                        (rank + step) as c_int,
                        7001,
                        self.comm(),
                        std::ptr::null_mut(),
                    )
                })?;
                for (a, b) in out.iter_mut().zip(&incoming) {
                    *a += *b;
                }
            }
            step <<= 1;
        }
        self.bcast_f64(out, 0)
    }
    pub fn alltoall_f64(
        &self,
        send: &[f64],
        block: usize,
        recv: &mut [f64],
    ) -> Result<(), MpiError> {
        let size = self.size();
        if send.len() != block * size || recv.len() != block * size {
            return Err(MpiError::Shape("alltoall lengths"));
        }
        // SAFETY: both buffers hold block * size doubles and are borrowed for the call.
        check("MPI_Alltoall", unsafe {
            MPI_Alltoall(
                send.as_ptr().cast(),
                count(block, "alltoall block")?,
                self.double(),
                recv.as_mut_ptr().cast(),
                count(block, "alltoall block")?,
                self.double(),
                self.comm(),
            )
        })
    }
}
