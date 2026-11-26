#!/usr/bin/env python3
"""Run the prepared Si tutorial with the native binding (Python 3.9+)."""
import argparse
from pathlib import Path
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

from baandi import Options, Session


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    args = parser.parse_args()
    config = tomllib.loads(args.input.read_text())
    system = config["system"]
    # This example deliberately uses Bohr, matching the native Session boundary.
    if system["lattice_unit"] != "bohr":
        parser.error("this script expects the supplied Bohr input")
    potentials = {
        label: (args.input.parent / path).read_text()
        for label, path in config["pseudopotentials"].items()
    }
    options = Options()
    options.reference_strict = False
    options.cutoff_hartree = config["basis"]["ecutwfc_ry"] / 2.0

    def session():
        return Session(
            system["lattice_vectors"],
            [atom["species"] for atom in system["atoms"]],
            [atom["fractional"] for atom in system["atoms"]],
            potentials,
            options,
        )

    batch = session()
    batch.run()
    if not batch.result().converged:
        raise RuntimeError(f"batch did not converge: {batch.status}")
    print("batch energy (Ha):", batch.result().energy_hartree)
    print("forces (Ha/Bohr):", batch.forces())

    stepped = session()
    stepped.initialize()
    while stepped.status == "Initialized":
        report = stepped.step()
        print("map:", report.kind, "iteration:", report.iteration,
              "estimator (Ha):", report.estimator_ha)
    if not stepped.result().converged:
        raise RuntimeError(f"stepwise did not converge: {stepped.status}")
    print("stepwise energy (Ha):", stepped.result().energy_hartree)


if __name__ == "__main__":
    main()
