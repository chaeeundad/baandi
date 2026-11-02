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
        Command::InspectUpf { .. } => unreachable!(),
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
        Command::InspectUpf { .. } => unreachable!(),
        Command::Grid { .. } => print!("{}", grid_report(&config, &prepared)),
    }
    Ok(())
}
fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
