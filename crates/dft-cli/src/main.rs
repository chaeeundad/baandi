#[path = "bands.rs"]
mod bands;
#[path = "dos.rs"]
mod dos;
#[path = "nscf.rs"]
mod nscf;
mod refscf_check;
mod refscf_map;
use dft_cli::checkpoint_io;
#[path = "pseudo.rs"]
mod pseudo;
use clap::{Parser, Subcommand};
use dft_cli::{config::Config, grid_report};
use std::{path::PathBuf, process::ExitCode};

#[derive(Parser)]
#[command(
    name = "baandi",
    version,
    about = "BAANDI plane-wave Kohn-Sham DFT engine"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}
#[derive(Subcommand)]
enum Command {
    /// Compare fixed-potential output densities from QE's own input trajectory.
    RefscfMap(refscf_map::Arguments),
    /// Replay independent Broyden mixing against a QE density trajectory.
    RefscfCheck {
        input: PathBuf,
        #[arg(
            long,
            required_unless_present = "stop_to_verified",
            conflicts_with = "stop_to_verified"
        )]
        qe: Option<PathBuf>,
        /// Read a completed reference-driver results.json; no SCF or QE replay.
        #[arg(long)]
        stop_to_verified: bool,
        #[arg(long, default_value_t = 6)]
        steps: usize,
        #[arg(long, default_value_t = 1e-12)]
        atol: f64,
        #[arg(long, default_value_t = 1e-9)]
        rtol: f64,
    },
    /// Relax atomic positions (and optionally the cell) with forces from verified SCF states.
    Relax {
        input: PathBuf,
        #[arg(long)]
        output_dir: PathBuf,
        #[arg(long)]
        resume: bool,
    },
    /// Validate a proposal and record an explicit review; never edits the input.
    PreflightReview {
        input: PathBuf,
        proposal: PathBuf,
        #[arg(long)]
        output_dir: Option<PathBuf>,
        #[arg(long, value_enum, requires_all = ["output_dir", "actor", "reason"])]
        decision: Option<dft_cli::config::preflight::proposal::Decision>,
        #[arg(long, requires = "decision")]
        actor: Option<String>,
        #[arg(long, requires = "decision")]
        reason: Option<String>,
    },
    /// Save an iteration boundary or an evaluated result in checkpoint v1 format.
    Checkpoint {
        input: PathBuf,
        output: PathBuf,
        #[arg(long, default_value_t = 3)]
        iteration: usize,
        #[arg(long, value_enum, default_value_t = checkpoint_io::Kind::Exact)]
        kind: checkpoint_io::Kind,
    },
    /// Resume an exact checkpoint, or explicitly consume a warm-start artifact.
    Restart {
        input: PathBuf,
        checkpoint: PathBuf,
        output: PathBuf,
        #[arg(long)]
        warm: bool,
    },

    /// SCF followed by a fixed-density band path and portable plotting files.
    Bands {
        input: PathBuf,
        #[arg(long)]
        output_dir: Option<PathBuf>,
    },
    /// SCF followed by a many-band fixed-density solve ([nscf] table; R1):
    /// dense, Chebyshev or Davidson bands with certified residuals.
    Nscf {
        input: PathBuf,
        #[arg(long)]
        output_dir: Option<PathBuf>,
        /// Per-k nscf_bands_v1 checkpoint units; existing matching units are
        /// resumed (re-certified), mismatching units reject.
        #[arg(long)]
        checkpoint_dir: Option<PathBuf>,
    },
    /// SCF followed by a weighted fixed-density Gaussian density of states.
    Dos {
        input: PathBuf,
        #[arg(long)]
        output_dir: Option<PathBuf>,
    },
    /// Evaluate a unitless expression (15 significant digits).
    Eval {
        #[arg(allow_hyphen_values = true)]
        expression: String,
    },
    /// Convert energy, length or pressure (15 significant digits).
    Convert {
        #[arg(allow_hyphen_values = true)]
        value: String,
        from: String,
        to: String,
    },
    /// Run SCF and write reproducible output files (Classic by default).
    Run {
        input: PathBuf,
        #[arg(long)]
        output_dir: Option<PathBuf>,
        #[arg(long, conflicts_with = "agent")]
        observe: bool,
        #[arg(long, conflicts_with = "observe")]
        agent: bool,
        /// JSONL mirror: '-' for stdout, or 'unix:/path' for a Unix datagram socket.
        #[arg(long)]
        events: Option<String>,
        /// Local proposal socket; defaults to OUTPUT/control/control.sock in agent mode.
        #[arg(long, requires = "agent")]
        control_socket: Option<PathBuf>,
    },
    /// Inspect the last durable event or the final results of a calculation.
    Status { calc: PathBuf },
    /// Explain a TOML configuration/preflight, or recorded calculation evidence.
    Explain { source: PathBuf },
    /// Inspect a UPF v2 NC pseudopotential and report strict support validation.
    InspectUpf {
        file: PathBuf,
        /// Hash-matched dataset provenance; never authorizes a calculation.
        #[arg(long)]
        manifest: Option<PathBuf>,
        #[arg(long)]
        json: bool,
    },
    /// Offline pseudopotential registry and verified cache.
    Pseudo {
        #[arg(long)]
        registry: Option<PathBuf>,
        #[command(subcommand)]
        command: pseudo::Command,
    },
    /// Full preflight including UPF identity and SCF controls; no SCF execution.
    #[command(alias = "check")]
    Validate {
        input: PathBuf,
        #[arg(long)]
        json: bool,
        #[arg(long)]
        explain: bool,
    },
    /// Print cell, cutoff, basis counts and FFT dimensions.
    Grid { input: PathBuf },
}
fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    match &cli.command {
        Command::RefscfMap(args) => return refscf_map::execute(args).map_err(Into::into),
        Command::RefscfCheck {
            input,
            qe,
            stop_to_verified,
            steps,
            atol,
            rtol,
        } => {
            if *stop_to_verified {
                return refscf_check::stop_to_verified(input).map_err(Into::into);
            }
            return refscf_check::execute(
                input,
                qe.as_deref().expect("clap requires --qe for replay"),
                *steps,
                refscf_check::Tolerance {
                    absolute: *atol,
                    relative: *rtol,
                },
            )
            .map_err(Into::into);
        }
        Command::Validate {
            input,
            json,
            explain,
        } => return dft_cli::config::preflight::cli_check(input, *json, *explain),
        Command::PreflightReview {
            input,
            proposal,
            output_dir,
            decision,
            actor,
            reason,
        } => {
            if let Some(decision) = decision {
                let path = dft_cli::config::preflight::proposal::record(
                    input,
                    proposal,
                    output_dir.as_deref().expect("clap requires output"),
                    *decision,
                    actor.as_deref().expect("clap requires actor"),
                    reason.as_deref().expect("clap requires reason"),
                )?;
                println!("Recorded {}; input unchanged", path.display());
                return Ok(());
            }
            return dft_cli::config::preflight::proposal::preview(input, proposal);
        }
        Command::Status { calc } => return dft_cli::inspect::status(calc),
        Command::Explain { source } => return dft_cli::inspect::explain(source),
        Command::Checkpoint {
            input,
            output,
            iteration,
            kind,
        } => return checkpoint_io::create(input, output, *iteration, *kind),
        Command::Restart {
            input,
            checkpoint,
            output,
            warm,
        } => return checkpoint_io::restart(input, checkpoint, output, *warm),
        _ => {}
    }
    if let Command::Pseudo { command, registry } = &cli.command {
        return pseudo::execute(command, registry.as_deref());
    }
    if let Command::InspectUpf {
        file,
        json,
        manifest,
    } = &cli.command
    {
        if *json {
            return pseudo::inspect_json(file, manifest.as_deref());
        }
        let report = pseudo::inspection_report(file, manifest.as_deref())?;
        println!(
            "SHA-256: {}\nHeader facts: {}\nLicense reference: {}\nQualification: {}",
            report["sha256"],
            report["header_facts"],
            report["license_reference"],
            report["qualification_status"]
        );
        for blocker in report["blockers"].as_array().into_iter().flatten() {
            println!(
                "Blocker: {}",
                blocker.as_str().unwrap_or("invalid metadata")
            );
        }
        let p = dft_upf::from_file(file)
            .map_err(|e| format!("{}: UNSUPPORTED or INVALID: {e}", file.display()))?;
        println!(
            "Element: {}\nz_valence: {}\nFunctional: {}\nType: NC",
            p.element(),
            p.z_valence(),
            p.functional()
        );
        println!(
            "Mesh: {} points, r = {} .. {} Bohr",
            p.mesh().r().len(),
            p.mesh().r()[0],
            p.mesh().r().last().unwrap()
        );
        println!("Projectors: {}", p.projectors().len());
        for (i, b) in p.projectors().iter().enumerate() {
            println!(
                "  {}: l={}, cutoff index={}",
                i + 1,
                b.l(),
                b.cutoff_index()
            );
        }
        println!(
            "D matrix: {} x {}, Ha, row-major {:?}",
            p.projectors().len(),
            p.projectors().len(),
            p.dij_hartree()
        );
        println!("Support: SUPPORTED (functional compatibility is decided by the engine)");
        return Ok(());
    }
    let input = match &cli.command {
        Command::Validate { input, .. } | Command::Grid { input } => input,
        Command::Status { .. }
        | Command::PreflightReview { .. }
        | Command::Explain { .. }
        | Command::Checkpoint { .. }
        | Command::Restart { .. }
        | Command::Eval { .. }
        | Command::Pseudo { .. }
        | Command::InspectUpf { .. }
        | Command::Bands { .. }
        | Command::Nscf { .. }
        | Command::Dos { .. }
        | Command::Run { .. }
        | Command::Relax { .. }
        | Command::RefscfCheck { .. }
        | Command::RefscfMap(_)
        | Command::Convert { .. } => {
            unreachable!()
        }
    };
    let text = std::fs::read_to_string(input)
        .map_err(|e| std::io::Error::new(e.kind(), format!("{}: {e}", input.display())))?;
    let config = Config::parse_at(&text, input)?;
    let prepared = config.prepare()?;
    match cli.command {
        Command::Validate { .. } => unreachable!(),
        Command::Status { .. }
        | Command::PreflightReview { .. }
        | Command::Explain { .. }
        | Command::Checkpoint { .. }
        | Command::Restart { .. }
        | Command::Eval { .. }
        | Command::Pseudo { .. }
        | Command::InspectUpf { .. }
        | Command::Bands { .. }
        | Command::Nscf { .. }
        | Command::Dos { .. }
        | Command::Run { .. }
        | Command::Relax { .. }
        | Command::RefscfCheck { .. }
        | Command::RefscfMap(_)
        | Command::Convert { .. } => {
            unreachable!()
        }
        Command::Grid { .. } => print!("{}", grid_report(&config, &prepared)),
    }
    Ok(())
}
fn main() -> ExitCode {
    use dft_cli::run_config::launch;
    let main_entry = std::time::Instant::now();
    // MPI lives for the whole process; without feature `mpi` this is a one-rank stand-in.
    let universe = match dft_mpi::Universe::init() {
        Ok(u) => u,
        Err(e) => {
            eprintln!("error: MPI initialisation: {e}");
            return ExitCode::from(5);
        }
    };
    let mpi_ready = std::time::Instant::now();
    let world = universe.world();
    let cli = Cli::parse();
    let args_ready = std::time::Instant::now();
    // RUN-CONFIG: resolve the input file's [performance]/[profile] tables with
    // the environment (environment > input > default) and install them as the
    // process configuration before the first switch is read.
    if let Err(e) = dft_cli::run_config::install_for_input(input_of(&cli.command)) {
        eprintln!("error: {e}");
        return ExitCode::from(1);
    }
    let config_ready = std::time::Instant::now();
    dft_cli::run::reference_launch_begin(main_entry, mpi_ready, args_ready, config_ready);
    launch::begin(main_entry, mpi_ready, args_ready, config_ready);
    match dft_numerics::switches::choice(dft_numerics::switches::Id::MPI_TREE) {
        "rd" => dft_mpi::set_tree_sum(dft_mpi::TreeSum::RecursiveDoubling),
        "rhd" => dft_mpi::set_tree_sum(dft_mpi::TreeSum::HalvingDoubling),
        _ => {}
    }
    dft_cli::parallel::set(world.rank(), world.size());
    if world.size() > 1 {
        // Every rank must resolve the same switch values before any optional
        // collective; a mismatch is the same error on every rank.
        if let Err(e) = launch::timed("switch_agreement", || {
            dft_cli::run_config::agree(&world, dft_numerics::switches::RunConfig::current())
        }) {
            eprintln!("error: {e}");
            return ExitCode::from(1);
        }
        dft_engine::parallel::set_communicator(Some(std::sync::Arc::new(world)));
    }
    let code = real_main(cli);
    launch::mark("main_dispatch_return");
    if dft_numerics::switches::present(dft_numerics::switches::Id::DGEMM_TIMERS) {
        for (site, calls, nanos, flops) in dft_numerics::dgemm_site_report() {
            eprintln!(
                "dgemm-timer {site} calls {calls} seconds {:.6} gflops_total {:.3}",
                nanos as f64 * 1e-9,
                flops as f64 * 1e-9
            );
        }
    }
    // Complete all engine work and console output before MPI_Finalize. Retain
    // ordinary Rust destruction; no process::exit bypass of MPI or I/O drops.
    launch::timed("communicator_drop", || {
        dft_engine::parallel::set_communicator(None)
    });
    // The feature-off Universe is a zero-sized stand-in without a Drop impl.
    #[allow(clippy::drop_non_drop)]
    launch::timed("mpi_finalize", || drop(universe));
    if let Err(e) = launch::finish() {
        eprintln!("error: launch profile output: {e}");
        return ExitCode::from(5);
    }
    code
}
/// The input file of a command that reads one (RUN-CONFIG).
fn input_of(command: &Command) -> Option<&std::path::Path> {
    match command {
        Command::RefscfMap(args) => Some(&args.input),
        Command::RefscfCheck {
            stop_to_verified: true,
            ..
        } => None,
        Command::Relax { input, .. }
        | Command::Checkpoint { input, .. }
        | Command::Restart { input, .. }
        | Command::Bands { input, .. }
        | Command::Nscf { input, .. }
        | Command::Dos { input, .. }
        | Command::Run { input, .. }
        | Command::Validate { input, .. }
        | Command::Grid { input }
        | Command::RefscfCheck { input, .. } => Some(input),
        Command::PreflightReview { .. }
        | Command::Eval { .. }
        | Command::Convert { .. }
        | Command::Status { .. }
        | Command::Explain { .. }
        | Command::InspectUpf { .. }
        | Command::Pseudo { .. } => None,
    }
}
fn real_main(cli: Cli) -> ExitCode {
    let post = match &cli.command {
        Command::RefscfMap(args) => Some(refscf_map::execute(args)),
        Command::Relax {
            input,
            output_dir,
            resume,
        } => Some(dft_cli::relax::execute(input, output_dir, *resume)),
        Command::Bands { input, output_dir } => Some(bands::execute(input, output_dir.as_deref())),
        Command::Dos { input, output_dir } => Some(dos::execute(input, output_dir.as_deref())),
        Command::Nscf {
            input,
            output_dir,
            checkpoint_dir,
        } => Some(nscf::execute(
            input,
            output_dir.as_deref(),
            checkpoint_dir.as_deref(),
        )),
        _ => None,
    };
    if let Some(result) = post {
        return match result {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("error: {e}");
                ExitCode::from(e.code)
            }
        };
    }
    let evaluation = match &cli.command {
        Command::Convert { value, from, to } => Some(convert(value, from, to)),
        Command::Eval { expression } => {
            Some(dft_units::eval_expr(expression).map_err(|e| e.to_string()))
        }
        _ => None,
    };
    if let Some(result) = evaluation {
        return match result {
            Ok(v) => {
                println!("{v:.14e}");
                ExitCode::SUCCESS
            }
            Err(e) => {
                eprintln!("error: {e}");
                ExitCode::from(2)
            }
        };
    }
    if let Command::Run {
        input,
        output_dir,
        observe,
        agent,
        events,
        control_socket,
    } = &cli.command
    {
        let mode = dft_cli::telemetry::RunMode {
            observe: *observe,
            agent: *agent,
            events: events.clone(),
            control_socket: control_socket.clone(),
        };
        return match dft_cli::run::execute_mode(input, output_dir.as_deref(), &mode) {
            Ok(report) => {
                use std::io::Write;
                let console_result = dft_cli::run_config::launch::timed("console_write", || {
                    if mode.events.as_deref() == Some("-") {
                        std::io::stderr()
                            .lock()
                            .write_all(report.console.as_bytes())
                    } else {
                        std::io::stdout()
                            .lock()
                            .write_all(report.console.as_bytes())
                    }
                });
                if let Err(e) = console_result {
                    eprintln!("error: console I/O: {e}");
                    return ExitCode::from(5);
                }
                if report.converged {
                    ExitCode::SUCCESS
                } else {
                    eprintln!("error: NOT CONVERGED: iteration budget exhausted");
                    ExitCode::from(3)
                }
            }
            Err(e) => {
                // Drain source-only diagnostics even when SCF fails, so the
                // failed iteration/root remains available to external oracles.
                let _ = dft_engine::work_accounting::take();
                eprintln!("error: {e}");
                ExitCode::from(e.code)
            }
        };
    }
    match run(cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

fn convert(value: &str, from: &str, to: &str) -> Result<f64, String> {
    use dft_units::*;
    let text = format!("{value} {from}");
    // Parsing the value separately prevents a unit embedded in the value argument.
    let v = eval_expr(value).map_err(|e| e.to_string())?;
    let result = if let Ok(unit) = from.parse::<EnergyUnit>() {
        convert_energy(v, unit, to.parse().map_err(|e: UnitError| e.to_string())?)
    } else if let Ok(unit) = from.parse::<LengthUnit>() {
        convert_length(v, unit, to.parse().map_err(|e: UnitError| e.to_string())?)
    } else if let Ok(unit) = from.parse::<PressureUnit>() {
        convert_pressure(v, unit, to.parse().map_err(|e: UnitError| e.to_string())?)
    } else {
        return Err(format!("unknown unit in {text:?}"));
    };
    if result.is_finite() {
        Ok(result)
    } else {
        Err(format!("nonfinite conversion of {text:?} to {to:?}"))
    }
}
