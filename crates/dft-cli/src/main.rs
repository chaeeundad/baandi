#[path = "bands.rs"]
mod bands;
#[path = "dos.rs"]
mod dos;
use dft_cli::checkpoint_io;
#[path = "pseudo.rs"]
mod pseudo;
use clap::{Parser, Subcommand};
use dft_cli::{config::Config, grid_report};
use std::{path::PathBuf, process::ExitCode};

#[derive(Parser)]
#[command(
    version,
    about = "Gamma-only plane-wave input validation and FFT grid inspection"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}
#[derive(Subcommand)]
enum Command {
    /// Relax atomic positions at fixed cell with strict certified forces.
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
        | Command::Dos { .. }
        | Command::Run { .. }
        | Command::Relax { .. }
        | Command::Convert { .. } => {
            unreachable!()
        }
    };
    let text = std::fs::read_to_string(input)
        .map_err(|e| std::io::Error::new(e.kind(), format!("{}: {e}", input.display())))?;
    let config = Config::parse(&text)?;
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
        | Command::Dos { .. }
        | Command::Run { .. }
        | Command::Relax { .. }
        | Command::Convert { .. } => {
            unreachable!()
        }
        Command::Grid { .. } => print!("{}", grid_report(&config, &prepared)),
    }
    Ok(())
}
fn main() -> ExitCode {
    // MPI lives for the whole process; without feature `mpi` this is a one-rank stand-in.
    let universe = match dft_mpi::Universe::init() {
        Ok(u) => u,
        Err(e) => {
            eprintln!("error: MPI initialisation: {e}");
            return ExitCode::from(5);
        }
    };
    let world = universe.world();
    dft_cli::parallel::set(world.rank(), world.size());
    if world.size() > 1 {
        dft_engine::parallel::set_communicator(Some(std::sync::Arc::new(world)));
    }
    let code = real_main();
    if dft_numerics::switches::present(dft_numerics::switches::Id::DGEMM_TIMERS) {
        for (site, calls, nanos, flops) in dft_numerics::dgemm_site_report() {
            eprintln!(
                "dgemm-timer {site} calls {calls} seconds {:.6} gflops_total {:.3}",
                nanos as f64 * 1e-9,
                flops as f64 * 1e-9
            );
        }
    }
    code
}
fn real_main() -> ExitCode {
    let cli = Cli::parse();
    let post = match &cli.command {
        Command::Relax {
            input,
            output_dir,
            resume,
        } => Some(dft_cli::relax::execute(input, output_dir, *resume)),
        Command::Bands { input, output_dir } => Some(bands::execute(input, output_dir.as_deref())),
        Command::Dos { input, output_dir } => Some(dos::execute(input, output_dir.as_deref())),
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
                let console_result = if mode.events.as_deref() == Some("-") {
                    std::io::stderr()
                        .lock()
                        .write_all(report.console.as_bytes())
                } else {
                    std::io::stdout()
                        .lock()
                        .write_all(report.console.as_bytes())
                };
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
