mod bands;
mod checkpoint_io;
mod dos;
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
    /// Run Gamma-only LDA-PZ SCF and write reproducible output files.
    Run {
        input: PathBuf,
        #[arg(long)]
        output_dir: Option<PathBuf>,
    },
    /// Inspect a UPF v2 NC pseudopotential and report strict support validation.
    InspectUpf {
        file: PathBuf,
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
    Validate {
        input: PathBuf,
        #[arg(long)]
        json: bool,
    },
    /// Print cell, cutoff, basis counts and FFT dimensions.
    Grid { input: PathBuf },
}
fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    match &cli.command {
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
    if let Command::InspectUpf { file, json } = &cli.command {
        if *json {
            return pseudo::inspect_json(file);
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
        Command::Checkpoint { .. }
        | Command::Restart { .. }
        | Command::Eval { .. }
        | Command::Pseudo { .. }
        | Command::InspectUpf { .. }
        | Command::Bands { .. }
        | Command::Dos { .. }
        | Command::Run { .. }
        | Command::Convert { .. } => {
            unreachable!()
        }
    };
    let text = std::fs::read_to_string(input)
        .map_err(|e| std::io::Error::new(e.kind(), format!("{}: {e}", input.display())))?;
    let config = Config::parse(&text)?;
    if let Command::Validate { json, .. } = &cli.command {
        let report = dft_cli::config::preflight::check(&config, input)?;
        if *json {
            println!("{}", serde_json::to_string_pretty(&report)?);
        } else {
            for record in &config.input_expressions {
                println!("{}: {:?} = {}", record.key_path, record.text, record.value);
            }
            report.print();
        }
        return Ok(());
    }
    let prepared = config.prepare()?;
    match cli.command {
        Command::Validate { .. } => unreachable!(),
        Command::Checkpoint { .. }
        | Command::Restart { .. }
        | Command::Eval { .. }
        | Command::Pseudo { .. }
        | Command::InspectUpf { .. }
        | Command::Bands { .. }
        | Command::Dos { .. }
        | Command::Run { .. }
        | Command::Convert { .. } => {
            unreachable!()
        }
        Command::Grid { .. } => print!("{}", grid_report(&config, &prepared)),
    }
    Ok(())
}
fn main() -> ExitCode {
    let cli = Cli::parse();
    let post = match &cli.command {
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
    if let Command::Run { input, output_dir } = &cli.command {
        return match dft_cli::run::execute(input, output_dir.as_deref()) {
            Ok(report) => {
                use std::io::Write;
                if let Err(e) = std::io::stdout()
                    .lock()
                    .write_all(report.console.as_bytes())
                {
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
