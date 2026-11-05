"""Optional ASE energy and force calculator; no free-energy alias.

Pseudopotentials are UPF text keyed by ASE chemical symbol. Copying the text
at construction avoids a file changing underneath a cached calculation.
"""
from ase.calculators.calculator import Calculator, PropertyNotImplementedError, all_changes

from . import Options, Session


class ConvergenceError(RuntimeError):
    """SCF exhausted its budget; the session retains its uncertified result."""


class Baandi(Calculator):
    implemented_properties = ["energy", "forces", "stress"]
    # Recompute instead of reusing a cache with incomplete model identity.
    use_cache = False

    def __init__(self, *, pseudopotentials, options=None, **kwargs):
        if any(name in kwargs for name in ("hubbard", "hubbard_u")):
            raise NotImplementedError("BAANDI ASE Hubbard energies and forces are unavailable")
        super().__init__(**kwargs)
        self._pseudopotentials = dict(pseudopotentials)
        source = options or Options()
        self._options = Options()
        for name in (
            "cutoff_hartree", "bands", "extra_bands", "max_iterations",
            "density_tolerance", "energy_tolerance", "eigen_tolerance", "seed",
            "mixing_beta", "mixing_history",
        ):
            setattr(self._options, name, getattr(source, name))
        self.session = None

    def get_property(self, name, atoms=None, allow_calculation=True):
        if name not in self.implemented_properties:
            raise PropertyNotImplementedError(
                f"BAANDI does not implement {name}; free_energy "
                "are unavailable (no zero placeholders)."
            )
        # Deliberately disable ASE's results cache, including changes to settings.
        self.results = {}
        return super().get_property(name, atoms, allow_calculation)

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        self.results = {}
        for name in properties:
            if name not in self.implemented_properties:
                raise PropertyNotImplementedError(f"BAANDI does not implement {name}")
        super().calculate(atoms, properties, system_changes)
        if self.atoms is None or not self.atoms.pbc.all():
            raise ValueError("BAANDI requires a fully periodic, nonsingular 3D cell")
        if self.atoms.get_initial_magnetic_moments().any():
            raise NotImplementedError("BAANDI supports only spin-unpolarized calculations")
        if self.atoms.get_initial_charges().any():
            raise NotImplementedError("BAANDI does not support charged ASE atoms")
        from ._baandi import angstrom_to_bohr, hartree_to_ev
        rows = [[angstrom_to_bohr(float(x)) for x in row] for row in self.atoms.cell]
        self.session = Session(
            rows, self.atoms.get_chemical_symbols(),
            self.atoms.get_scaled_positions(wrap=False).tolist(),
            self._pseudopotentials, self._options,
        )
        self.session.run()
        result = self.session.result()
        if not result.converged:
            raise ConvergenceError(
                f"SCF did not converge after {result.iterations} iterations; "
                "uncertified diagnostics are in calculator.session.result()"
            )
        self.results = {"energy": hartree_to_ev(result.energy_hartree)}
        if "forces" in properties:
            import numpy as np
            # (Ha/Bohr) * (eV/Ha) * (Bohr/Angstrom).
            self.results["forces"] = np.asarray(self.session.forces()) * (
                hartree_to_ev(1.0) * angstrom_to_bohr(1.0)
            )

        if "stress" in properties:
            import numpy as np
            from ase.stress import full_3x3_to_voigt_6_stress
            # ASE: tensile-positive, (xx, yy, zz, yz, xz, xy), no shear factor.
            tensor = np.asarray(self.session.stress()) * (
                hartree_to_ev(1.0) * angstrom_to_bohr(1.0) ** 3
            )
            self.results["stress"] = full_3x3_to_voigt_6_stress(tensor)

    def get_eigenvalues(self, kpt=0, spin=0):
        if kpt != 0 or spin != 0:
            raise NotImplementedError("Only Gamma, spin-unpolarized eigenvalues are available")
        if self.session is None or not self.session.result().converged:
            raise ConvergenceError("Run a converged energy calculation first")
        from ._baandi import hartree_to_ev
        import numpy as np
        return np.array([hartree_to_ev(x) for x in self.session.result().eigenvalues_hartree])
