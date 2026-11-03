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
    InspectUpf { file: PathBuf },
    /// Validate input and basis/grid settings; does not read UPF files.
    Validate { input: PathBuf },
    /// Print cell, cutoff, basis counts and FFT dimensions.
    Grid { input: PathBuf },
}
fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    if let Command::InspectUpf { file } = &cli.command {
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
        Command::Validate { input } | Command::Grid { input } => input,
        Command::InspectUpf { .. } | Command::Run { .. } | Command::Convert { .. } => {
            unreachable!()
        }
    };
    let text = std::fs::read_to_string(input)
        .map_err(|e| std::io::Error::new(e.kind(), format!("{}: {e}", input.display())))?;
    let config = Config::parse(&text)?;
    let prepared = config.prepare()?;
    match cli.command {
        Command::Validate { .. } => println!(
            "Valid: {} (pseudopotential paths parsed; file existence not checked)",
            config.system.name
        ),
        Command::InspectUpf { .. } | Command::Run { .. } | Command::Convert { .. } => {
            unreachable!()
        }
        Command::Grid { .. } => print!("{}", grid_report(&config, &prepared)),
    }
    Ok(())
}
fn main() -> ExitCode {
    let cli = Cli::parse();
    if let Command::Convert { value, from, to } = &cli.command {
        return match convert(value, from, to) {
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
    let v: f64 = value
        .parse()
        .map_err(|_| format!("invalid value {value:?}"))?;
    if !v.is_finite() {
        return Err(format!("nonfinite value {value:?}"));
    }
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
