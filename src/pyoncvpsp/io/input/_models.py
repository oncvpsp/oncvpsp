"""ONVPSP input file handling."""

import logging
from os import PathLike
from typing import Annotated, Any, Iterator, TextIO

from pydantic import BaseModel, Field, model_validator
import toml

from .._utils import flatten, recursive_merge

LOGGER: logging.Logger = logging.getLogger(__name__)


class OncvpspAtomInput(BaseModel):  # pylint: disable=too-few-public-methods
    r"""Input parameters defining the atom being pseudized.

    Attributes:
        atsym (str): Atomic symbol of the element (e.g. 'H', 'He', ...).
        z (float): Nuclear charge.
        nc (int): Number of core states.
        nv (int): Number of valence states.
        iexc (int): Exchange-correlation functional. Negative values indicate Libxc; 1-4 are built-in functionals.
    """

    atsym: Annotated[str, Field(description="Atomic symbol of the element (e.g. 'H', 'He', ...)")]
    z: Annotated[float, Field(gt=0.0, le=120.0, description="Atomic charge (number of protons)")]
    nc: Annotated[int, Field(ge=0, description="Number of core states")]
    nv: Annotated[int, Field(ge=1, description="Number of valence states")]
    iexc: Annotated[
        int,
        Field(
            le=4,
            description=(
                "Exchange-correlation functional. Negative values indicate Libxc. 1-4 are built-in functionals."
            ),
        ),
    ]

    @model_validator(mode="after")
    def _check_atomic_symbol(self) -> "OncvpspAtomInput":
        if len(self.atsym) == 0:
            raise ValueError("Atomic symbol cannot be empty")
        if len(self.atsym) > 2:
            raise ValueError("Atomic symbol must be at most 2 characters")
        return self

    @model_validator(mode="after")
    def _check_exchange_correlation_type(self) -> "OncvpspAtomInput":
        if self.iexc == 0:
            raise ValueError("exchange_correlation_type can be negative (LibXC), 1-4 (built-in), but not 0")
        return self

    def as_dat_string(self, psfile: str, with_comment: bool = False) -> str:
        """String representation of the pseudopotential input in .dat file format."""
        s: str = f"{self.atsym:>10s} {self.z:>10.5f} {self.nc:>10d} " f"{self.nv:>10d} {self.iexc:>10d} {psfile:>10s}"
        if with_comment:
            s = f"# {'atsym':>8s} {'z':>10s} {'nc':>10s} " f"{'nv':>10s} {'iexc':>10s} {'psfile':>10s}\n" + s
        return s


class AtomicStateInput(BaseModel):  # pylint: disable=too-few-public-methods
    r"""Input parameters describing an atomic state.

    Attributes:
        n (int): Principal quantum number (1 <= n).
        l (int): Angular momentum quantum number (0 <= l <= n - 1).
        f (float): Occupation number (0.0 <= f <= 2(2l + 1)).
    """

    n: Annotated[int, Field(ge=1)]
    l: Annotated[int, Field(ge=0)]
    f: Annotated[float, Field(ge=0.0)]

    @model_validator(mode="after")
    def _check(self) -> "AtomicStateInput":
        if self.l > self.n - 1:
            raise ValueError(f"Invalid angular momentum: (n={self.n}), l={self.l}")
        if self.f > 2 * (2 * self.l + 1):
            raise ValueError(f"Invalid occupation number: (n={self.n}, l={self.l}), f={self.f}")
        return self

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, AtomicStateInput):
            raise NotImplementedError(f"Cannot compare AtomicState with {type(other)}")
        if self.n != other.n:
            LOGGER.debug("n differs: %d != %d", self.n, other.n)
            return False
        if self.l != other.l:
            LOGGER.debug("l differs: %d != %d", self.l, other.l)
            return False
        if self.f != other.f:
            LOGGER.debug("f differs: %f != %f", self.f, other.f)
            return False
        return True

    def as_dat_string(self, with_comment: bool = False) -> str:
        """String representation of the atomic state in .dat file format."""
        s: str = f"{self.n:>10d} {self.l:>10d} {self.f:>10.5f}"
        if with_comment:
            s = f"# {'n':>8s} {'l':>10s} {'f':>10s}\n" + s
        return s


class AtomicConfigurationInput(BaseModel):  # pylint: disable=too-few-public-methods
    r"""Input parameters describing a full atomic configuration (composed of multiple atomic states).

    Attributes:
        n (list[int]): List of principal quantum numbers for each state.
        l (list[int]): List of angular momentum quantum numbers for each state.
        f (list[float]): List of occupation numbers for each state.
    """

    n: Annotated[list[int], Field(min_length=1)]
    l: Annotated[list[int], Field(min_length=1)]
    f: Annotated[list[float], Field(min_length=1)]

    def __len__(self) -> int:
        return len(self.n)

    def __getitem__(self, index: int) -> AtomicStateInput:
        return AtomicStateInput(
            n=self.n[index],
            l=self.l[index],
            f=self.f[index],
        )

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    @property
    def charge(self) -> float:
        """Total occupation (charge) of the atomic configuration."""
        return sum(self.f)

    @model_validator(mode="after")
    def _check(self) -> "AtomicConfigurationInput":
        total_occupation = sum(self.f)
        if total_occupation != self.charge:
            raise ValueError(f"Total occupation {total_occupation} does not match charge {self.charge}")
        return self

    def as_dat_string(self, with_comment: bool = False) -> str:
        """String representation of the atomic state in .dat file format."""
        states = sorted(list(self), key=lambda s: (s.n, s.l))
        lines = [states[0].as_dat_string(with_comment)]
        if len(states) > 1:
            for state in states[1:]:
                lines.append(state.as_dat_string())
        return "\n".join(lines)


class PseudopotentialInput(BaseModel):  # pylint: disable=too-few-public-methods
    r"""Input parameters for a single pseudopotential channel.

    Attributes:
        l (int): Angular momentum quantum number.
        rc (float): [Bohr] Core radius beyond which the pseudopotential matches the all-electron potential.
            Correspondingly, the radius beyond which the pseudo-wavefunctions for this angular momentum
            match the all-electron wavefunctions.
        ep (float): [Hartree] Energy at which the pseudopotential is constructed.
        ncon (int): Number of constraints used in matching the pseudo-wavefunction to the all-electron
            wavefunction at rc. :math:`d^{i}\phi/dr^{i} = d^{i}\psi/dr^{i}` for :math:`i=0,...,n_{\mathrm{con}}-1`.
        nbas (int): Number of basis functions used in optimizing the pseudo-wavefunction (ncon + 2 <= nbas <= ncon + 5).
        qcut (float): [1/Bohr] Cutoff wavevector used for defining the residual kinetic energy.
    """

    l: Annotated[int, Field(ge=0)]
    rc: Annotated[float, Field(gt=0.0)]
    ep: float
    ncon: Annotated[int, Field(ge=3, le=5)]
    nbas: Annotated[int, Field(ge=5, le=10)]
    qcut: Annotated[float, Field(gt=0.0)]

    @model_validator(mode="after")
    def _check(self) -> "PseudopotentialInput":
        if not self.ncon + 2 <= self.nbas <= self.ncon + 5:
            raise ValueError(
                f"n_basis must be between n_constraints + 2 and n_constraints + 5: "
                f"n_basis={self.nbas}, n_constraints={self.ncon}"
            )
        return self

    def as_dat_string(self, with_comment: bool = False) -> str:
        """String representation of the pseudopotential input in .dat file format."""
        s: str = f"{self.l:>10d} {self.rc:>10.5f} {self.ep:10.5f} {self.ncon:>10d} {self.nbas:>10d} {self.qcut:>10.5f}"
        if with_comment:
            s = f"# {'l':>8s} {'rc':>10s} {'ep':>10s} {'ncon':>10s} {'nbas':>10s} {'qcut':>10s}\n" + s
        return s


class PseudoPotentialsInput(BaseModel):  # pylint: disable=too-few-public-methods
    r"""Input parameters for all pseudopotential channels.

    Attributes:
        lmax (int): Maximum angular momentum channel of the pseudopotential.
        l (list[int]): List of angular momentum quantum numbers for each channel (`all(li <= lmax for li in l)`).
        rc (list[float]): List of core radii for each channel.
        ep (list[float]): List of energies at which each pseudopotential channel is constructed.
        ncon (list[int]): List of number of constraints for each channel.
        nbas (list[int]): List of number of basis functions for each channel.
        qcut (list[float]): List of cutoff wavevectors for each channel.
    """

    lmax: Annotated[
        int,
        Field(ge=0, description="Maximum angular momentum channel of the pseudopotential."),
    ]
    l: Annotated[list[int], Field(min_length=1)]
    rc: Annotated[list[float], Field(min_length=1)]
    ep: Annotated[list[float], Field(min_length=1)]
    ncon: Annotated[list[int], Field(min_length=1)]
    nbas: Annotated[list[int], Field(min_length=1)]
    qcut: Annotated[list[float], Field(min_length=1)]

    def __len__(self) -> int:
        return self.lmax + 1

    def __getitem__(self, index: int) -> PseudopotentialInput:
        return PseudopotentialInput(
            l=self.l[index],
            rc=self.rc[index],
            ep=self.ep[index],
            ncon=self.ncon[index],
            nbas=self.nbas[index],
            qcut=self.qcut[index],
        )

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    @model_validator(mode="after")
    def _check(self) -> "PseudoPotentialsInput":
        if sorted(self.l) != list(range(self.lmax + 1)):
            raise ValueError(f"Pseudo potentials must cover all l from 0 to {self.lmax}")
        if not len(self.l) == len(self.rc) == len(self.ep) == len(self.ncon) == len(self.nbas) == len(self.qcut):
            raise ValueError("All pseudopotential parameter lists must have the same length")
        for l in self.l:
            if l > self.lmax:
                raise ValueError(f"Pseudopotential with l={l} exceeds l_max={self.lmax}")
        return self

    def as_dat_string(self, with_comment: bool = False) -> str:
        """String representation of the pseudopotential input in .dat file format."""
        psps: list[PseudopotentialInput] = [self[l] for l in range(self.lmax + 1)]
        lines: list[str] = [
            f"# {'lmax':>8s}\n{self.lmax:>10d}\n#",
            psps[0].as_dat_string(with_comment),
        ]
        if len(psps) > 1:
            for psp in psps[1:]:
                lines.append(psp.as_dat_string())
        return "\n".join(lines)


class LocalPotentialInput(BaseModel):  # pylint: disable=too-few-public-methods
    """Input parameters for the local part of the KB separable form of the pseudopotential

    Attributes:
        lloc (int): Angular momentum channel used to construct the local potential. If `lloc=4`,
            then a polynomial continuation of the all-electron potential from `r=rcloc` to `r=0` is used.
        lpopt (None | int): If `lloc=4`, specifies the type of polynomial used in the continuation (1 <= lpopt <= 5).
            1) Match the 0th, 2nd, and 4th radial derivatives at `rcloc`.
            2) Match the 0th, 4th, and 6th radial derivatives at `rcloc`.
            3) Match the 0th, 4th, 5th, and 6th radial derivatives at `rcloc`.
            4) Match the 0th, 4th, 6th, and 8th radial derivatives at `rcloc`.
            5) Match the 0th, 2nd, 4th, and 6th radial derivatives at `rcloc`.
        rcloc (None | float): [Bohr] The radius beyond which the local potential matches the all-electron potential.
            If `lloc != 4`, `all(rcloc < rc[l] for l in range(lmax + 1))` must hold.
            If `lloc == 4`, `all(rcloc <= rc[l] for l in range(lmax + 1))` must hold.
        dvloc0 (None | float): [Hartree] If `lloc=4`, shifts the local potential at `r = 0` depending on `lpopt`:
            1) Shift by `dvloc0 * (1-(r/rcloc)**2)**3`
            2) Shift by `dvloc0 * (1-(r/rcloc)**4)**3`
            3) Shift by `dvloc0 * (1-(r/rcloc)**4)**4`
            4) Shift by `dvloc0 * (1-(r/rcloc)**4)**4`
            5) Shift by `dvloc0 * (1-(r/rcloc)**4)**4`
    """

    lloc: Annotated[int, Field(ge=0)]
    lpopt: None | Annotated[int, Field(ge=1, le=5)] = None
    rcloc: None | Annotated[float, Field(gt=0.0)] = None
    dvloc0: None | Annotated[float, Field(ge=0.0)] = None

    @model_validator(mode="after")
    def _check(self) -> "LocalPotentialInput":
        if self.lloc == 4:
            exceptions = []
            if self.lpopt is None:
                exceptions.append(ValueError("lpopt must be set when lloc=4"))
            if self.rcloc is None:
                exceptions.append(ValueError("rcloc must be set when lloc=4"))
            if self.dvloc0 is None:
                exceptions.append(ValueError("dvloc0 must be set when lloc=4"))
            if exceptions:
                raise ExceptionGroup("Invalid LocalPotentialInput for lloc=4", exceptions)
        return self

    def as_dat_string(self, with_comment: bool = False) -> str:
        """String representation of the pseudopotential input in .dat file format."""
        lpopt: int = self.lpopt if self.lpopt else 0
        rc: float = self.rcloc if self.rcloc else 0.0
        dvloc0: float = self.dvloc0 if self.dvloc0 else 0.0
        s: str = f"{self.lloc:>10d} {lpopt:>10d} {rc:>10.5f} {dvloc0:>10.5f}"
        if with_comment:
            s = f"# {'lloc':>8s} {'lpopt':>10s} {'rc(5)':>10s} {'dvloc0':>10s}\n" + s
        return s


class VKBProjectorInput(BaseModel):  # pylint: disable=too-few-public-methods
    """Input parameters for the Vanderbilt / Kleinman-Bylander projector(s) of one angular momentum channel.

    Attributes:
        l (int): Angular momentum quantum number.
        nproj (int): Number of projectors for this angular momentum channel.
        debl (float): [Hartree] Energy added to that of the previous pseudo-wavefunction for constructing
            the 2nd (and subsequent) projector(s) if no corresponding occupied bound state exists.
    """

    l: Annotated[int, Field(ge=0)]
    nproj: Annotated[int, Field(ge=0, le=5)]
    debl: Annotated[float, Field(ge=0.0)] = 0.0

    def as_dat_string(self, with_comment: bool = False) -> str:
        """String representation of the pseudopotential input in .dat file format."""
        s: str = f"{self.l:>10d} {self.nproj:>10d} {self.debl:>10.5f}"
        if with_comment:
            s = f"# {'l':>8s} {'nproj':>10s} {'debl':>10s}\n" + s
        return s


class VKBProjectorsInput(BaseModel):  # pylint: disable=too-few-public-methods
    """Input parameters for all Vanderbilt/ Kleinman-Bylander projectors.

    Attributes:
        l (list[int]): List of angular momentum quantum numbers for each channel (`all(li <= lmax for li in l)`).
        nproj (list[int]): List of number of projectors for each channel.
        debl (list[float]): List of energy increments for constructing additional projectors for each channel.
    """

    l: Annotated[list[int], Field(min_length=1)]
    nproj: Annotated[list[int], Field(min_length=1)]
    debl: Annotated[list[float], Field(min_length=1)]

    def __len__(self) -> int:
        return len(self.l)

    def __getitem__(self, index: int) -> VKBProjectorInput:
        return VKBProjectorInput(
            l=self.l[index],
            nproj=self.nproj[index],
            debl=self.debl[index],
        )

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    @model_validator(mode="after")
    def _check(self) -> "VKBProjectorsInput":
        if sorted(self.l) != list(range(max(self.l) + 1)):
            raise ValueError(f"VKB projectors must cover all l from 0 to {max(self.l)}")
        return self

    def as_dat_string(self, with_comment: bool = False) -> str:
        """String representation of the pseudopotential input in .dat file format."""
        projectors: list[VKBProjectorInput] = [self[l] for l in range(max(self.l) + 1)]
        lines: list[str] = [projectors[0].as_dat_string(with_comment)]
        if len(projectors) > 1:
            for p in projectors[1:]:
                lines.append(p.as_dat_string(with_comment=False))
        return "\n".join(lines)


class ModelCoreChargeInput(BaseModel):  # pylint: disable=too-few-public-methods
    r"""Input parameters controlling the construction of a model core charge density.

    Attributes:
        icmod (int): Model core charge type:
            0) No model core charge.
            1) Polynomial continuation of the all-electron core charge density. Matches the
                all-electron core charge density and its first 4 radial derivatives at
                the radius where :math:`\rho_{c}(r) = f_{c}^{\mathrm{fact}} \rho_v(r)`.
            2) Teter model core charge where the Teter parameters are adjusted such that
                the model core charge matches the all-electron core charge and its first
                derivative at the radius where :math:`\rho_{c}(r) = f_{c}^{\mathrm{fact}} \rho_v(r)`.
            3) Teter model core charge with fixed amplitude and scale factor (if `teter_relative=True`)
                or fixed amplitude and scale (if `teter_relative=False`).
            4) Optimized Teter model core charge where the Teter parameters are selected to minimize
                the RMSE in the 2nd derivative of the exchange-correlation energy w.r.t. the all-electron
                atom.
        fcfact (float): Used to determine the matching radius for `icmod in (1, 2)`. Unused for `icmod in (0, 3, 4)`.
        rcfact (float): Unused.
        teter_relative (bool): If True, `teter_amp*` and `teter_scale*` specify relative
            ranges for optimizing the Teter parameters. If False, they specify absolute ranges.
        teter_amp_min (None | float): Minimum Teter amplitude for grid search (relative or absolute).
        teter_amp_max (None | float): Maximum Teter amplitude for grid search (relative or absolute).
        teter_amp_step (None | float): Teter amplitude step size for grid search (relative or absolute).
        teter_scale_min (None | float): Minimum Teter scale factor for grid search (relative or absolute).
        teter_scale_max (None | float): Maximum Teter scale factor for grid search (relative or absolute).
        teter_scale_step (None | float): Teter scale factor step size for grid search (relative or absolute).
    """

    icmod: Annotated[int, Field(ge=0, le=4)]
    fcfact: Annotated[float, Field(ge=0.0)] = 0.0
    rcfact: Annotated[float, Field(ge=0.0)] = 0.0
    teter_relative: Annotated[bool, Field()] = True
    teter_amp_min: None | Annotated[float, Field()] = None
    teter_amp_max: None | Annotated[float, Field()] = None
    teter_amp_step: None | Annotated[float, Field()] = None
    teter_scale_min: None | Annotated[float, Field()] = None
    teter_scale_max: None | Annotated[float, Field()] = None
    teter_scale_step: None | Annotated[float, Field()] = None

    @model_validator(mode="after")
    def _check(self) -> "ModelCoreChargeInput":
        if self.icmod in (1, 3):
            if self.fcfact == 0.0:
                raise ValueError("fcfact must be greater than 0 for icmod 1 and 3")
        if self.icmod == 3:
            if self.rcfact == 0.0:
                raise ValueError("rcfact must be greater than 0 for icmod 3")
        return self

    def as_dat_string(self, with_comment: bool = False) -> str:
        """String representation of the pseudopotential input in .dat file format."""
        if not self.teter_relative:
            raise ValueError("Absolute Teter parameters not supported by text input")
        if self.icmod in (0, 1, 2):
            s = f"{self.icmod:>10d} {self.fcfact:>10.5f}"
            if with_comment:
                s = f"# {'icmod':>8s} {'fcfact':>10s}\n" + s
        elif self.icmod == 3:
            s = f" {self.icmod:>10d} {self.fcfact:>10.5f} {self.rcfact:>10.5f}"
            if with_comment:
                s = f"# {'icmod':>8s} {'fcfact':>10s} {'rcfact':>10s}\n" + s
        else:
            teter_grid_params = (
                self.teter_amp_min,
                self.teter_amp_max,
                self.teter_amp_step,
                self.teter_scale_min,
                self.teter_scale_max,
                self.teter_scale_step,
            )
            if all(v is None for v in teter_grid_params):
                s = f"{self.icmod:>10d} {self.fcfact:>10.5f}"
                if with_comment:
                    s = f"# {'icmod':>8s} {'fcfact':>10s}\n" + s
            else:
                s = (
                    f" {self.icmod:>10d} {self.fcfact:>10.5f} {self.rcfact:>10.5f} "
                    f"{self.teter_amp_min:>10.5f} {self.teter_amp_max:>10.5f} {self.teter_amp_step:>10.5f} "
                    f"{self.teter_scale_min:>10.5f} {self.teter_scale_max:>10.5f} {self.teter_scale_step:>10.5f}"
                )
                if with_comment:
                    s = (
                        f"# {'icmod':>8s} {'fcfact':>10s} {'rcfact':>10s} "
                        f"{'fcfact_min':>10s} {'fcfact_max':>10s} {'fcfact_step':>10s} "
                        f"{'rcfact_min':>10s} {'rcfact_max':>10s} {'rcfact_step':>10s}\n" + s
                    )
        return s


class LogDerivativeInput(BaseModel):  # pylint: disable=too-few-public-methods
    """Input parameters controlling the calculation of logarithmic derivatives.

    Attributes:
        epsh1 (float): [Hartree] Minimum energy for logarithmic derivative calculation.
        epsh2 (float): [Hartree] Maximum energy for logarithmic derivative calculation.
        depsh (float): [Hartree] Energy step size for logarithmic derivative calculation.
        rxpsh (None | float): [Bohr] Radius at which logarithmic derivatives are calculated.
            If `None`, the radius is chosen automatically at a few grid points beyond the largest
            core radius used in the pseudopotential construction.
    """

    epsh1: float
    epsh2: float
    depsh: Annotated[float, Field(gt=0.0)]
    rxpsh: None | Annotated[float, Field(gt=0.0)] = None

    @model_validator(mode="after")
    def _check(self) -> "LogDerivativeInput":
        if self.epsh2 <= self.epsh1:
            raise ValueError("e_max must be greater than e_min")
        return self

    def as_dat_string(self, with_comment: bool = False) -> str:
        """String representation of the pseudopotential input in .dat file format."""
        if self.rxpsh is None:
            s = f"{self.epsh1:>10.5f} {self.epsh2:>10.5f} {self.depsh:>10.5f}"
            if with_comment:
                s = f"# {'epsh1':>8s} {'epsh2':>10s} {'depsh':>10}\n" + s
        else:
            s = f"{self.epsh1:>10.5f} {self.epsh2:>10.5f} {self.depsh:>10.5f} {self.rxpsh:>10.5f}"
            if with_comment:
                s = f"# {'epsh1':>8s} {'epsh2':>10s} {'depsh':>10s} {'rxpsh':>10s}\n" + s
        return s


class OutputGridInput(BaseModel):  # pylint: disable=too-few-public-methods
    """Input parameters controlling the linear mesh (starting from `r=0`)
    used for writing the output pseudopotential file(s).

    Attributes:
        rlmax (float): [Bohr] Maximum radius of the output grid (exclusive).
        drl (float): [Bohr] Spacing of the output grid.
    """

    rlmax: Annotated[float, Field(gt=0.0)]
    drl: Annotated[float, Field(gt=0.0)]

    def as_dat_string(self, with_comment: bool = False) -> str:
        """String representation of the pseudopotential input in .dat file format."""
        s = f"{self.rlmax:>10.5f} {self.drl:>10.5f}"
        if with_comment:
            s = f"# {'rlmax':>8s} {'drl':>10s}\n" + s
        return s


class OncvpspPPOutputInput(BaseModel):  # pylint: disable=too-few-public-methods
    """Input parameter controlling the pseudopotential file formats written as output.

    Attributes:
        psfile (str): Pseudopotential output format: 'psp8', 'upf', 'psml', 'both', 'all', or 'none'.
            `both` writes 'psp8' and 'upf' formats.
            `all` writes 'psp8', 'upf', and 'psml' formats.
            `none` writes no pseudopotential files (useful for testing/ optimization runs).
    """

    psfile: Annotated[
        str,
        Field(description="Pseudopotential output format: 'psp8', 'upf', 'psml', 'both', 'all', or 'none'."),
    ]

    @model_validator(mode="after")
    def _check_psp_output_format(self) -> "OncvpspPPOutputInput":
        valid_formats = {"psp8", "upf", "psml", "both", "all", "none"}
        if self.psfile not in valid_formats:
            raise ValueError(
                f"Invalid pseudopotential output format: {self.psfile}. "
                f"Valid formats are: {', '.join(valid_formats)}"
            )
        return self


class OncvpspInput(BaseModel):
    """ONCVPSP Input parameters"""

    oncvpsp: OncvpspAtomInput
    reference_configuration: AtomicConfigurationInput
    pseudopotentials: PseudoPotentialsInput
    local_potential: LocalPotentialInput
    vkb_projectors: VKBProjectorsInput
    model_core_charge: ModelCoreChargeInput
    log_derivative_analysis: LogDerivativeInput
    test_configurations: list[AtomicConfigurationInput] = []
    linear_mesh: OutputGridInput
    pp_output: OncvpspPPOutputInput

    def merge(self, other: dict[str, Any] | None = None, **kwargs) -> "OncvpspInput":
        """Merge this set of input parameters with another set of parameters or keyword arguments.
        Keys in `other` and `kwargs` will override those in `self`.
        Keys in `kwargs` take precedence over both `other` and `self`.

        Returns:
            OncvpspInput: Merged ONCVPSP input parameters.
        """
        data = self.model_dump()
        merged = recursive_merge(data, other or {})
        merged = recursive_merge(merged, kwargs)
        return OncvpspInput(**merged)

    def as_dat_string(self, with_comment=True) -> str:
        """Write the input to a string in the format expected by ONCVPSP for stdin input.

        Args:
            with_comment (bool, optional): Add explanatory comments to the input string. Defaults to True.

        Returns:
            str: ONCVPSP stdin string
        """
        if with_comment:
            return "\n".join(
                [
                    "# ATOM AND REFERENCE CONFIGURATION",
                    self.oncvpsp.as_dat_string(psfile=self.pp_output.psfile, with_comment=with_comment),
                    "#",
                    self.reference_configuration.as_dat_string(with_comment=with_comment),
                    "#",
                    "# PSEUDOPOTENTIAL AND OPTIMIZATION",
                    self.pseudopotentials.as_dat_string(with_comment=with_comment),
                    "#",
                    "# LOCAL POTENTIAL",
                    self.local_potential.as_dat_string(with_comment=with_comment),
                    "#",
                    "# VANDERBILT-KLEINMAN-BYLANDER PROJECTORS",
                    self.vkb_projectors.as_dat_string(with_comment=with_comment),
                    "#",
                    "# MODEL CORE CHARGE",
                    self.model_core_charge.as_dat_string(with_comment=with_comment),
                    "#",
                    "# LOG DERIVATIVE ANALYSIS",
                    self.log_derivative_analysis.as_dat_string(with_comment=with_comment),
                    "#",
                    "# OUTPUT GRID",
                    self.linear_mesh.as_dat_string(with_comment=with_comment),
                    "#",
                    "# TEST CONFIGURATIONS",
                    f"# {'ncnf':>8s}",
                    f"{len(self.test_configurations):>10d}",
                ]
                + [
                    f"#\n# {'nvcnf':>8s}\n{len(config):>10d}\n" + config.as_dat_string(with_comment=with_comment)
                    for config in self.test_configurations
                ]
            )
        return "\n".join(
            [
                self.oncvpsp.as_dat_string(psfile=self.pp_output.psfile, with_comment=with_comment),
                self.reference_configuration.as_dat_string(with_comment=with_comment),
                self.pseudopotentials.as_dat_string(with_comment=with_comment),
                self.local_potential.as_dat_string(with_comment=with_comment),
                self.vkb_projectors.as_dat_string(with_comment=with_comment),
                self.model_core_charge.as_dat_string(with_comment=with_comment),
                self.log_derivative_analysis.as_dat_string(with_comment=with_comment),
                self.linear_mesh.as_dat_string(with_comment=with_comment),
                f"{len(self.test_configurations):>10d}",
            ]
            + [
                f"{len(config):>10d}\n" + config.as_dat_string(with_comment=with_comment)
                for config in self.test_configurations
            ]
        )

    @classmethod
    def from_dat_file(cls, file: str | PathLike | TextIO) -> "OncvpspInput":
        """Load the input from an ONCVPSP dat file.

        Args:
            file (PathLike | TextIO): Path to ONCVPSP dat file or file-like object.
        Returns:
            OncvpspInput: ONCVPSP input parameters.
        """
        if isinstance(file, (str, PathLike)):
            with open(file, "r", encoding="utf-8") as f:
                return cls.from_dat_string(f.read())
        else:
            return cls.from_dat_string(file.read())

    def to_dat_file(self, file: str | PathLike | TextIO, with_comment: bool = True) -> None:
        """Write the input to an ONCVPSP dat file.

        Args:
            file (PathLike | TextIO): Path to ONCVPSP dat file or file-like object.
            with_comment (bool, optional): Add explanatory comments to the input string. Defaults to True.
        """
        s = self.as_dat_string(with_comment=with_comment)
        if isinstance(file, (str, PathLike)):
            with open(file, "w", encoding="utf-8") as f:
                f.write(s)
        else:
            file.write(s)

    @classmethod
    def from_dat_string(cls, s: str) -> "OncvpspInput":  # pylint: disable=too-many-locals,too-many-statements
        """Load the input from an ONCVPSP stdin string.

        Args:
            s (str): Input in ONCVPSP stdin format

        Returns:
            OncvpspInput: ONCVPSP input parameters.
        """
        lines: list[str] = s.strip().splitlines()
        # Strip whitespace
        lines = [line.strip() for line in lines]
        # Remove empty lines
        lines = [line for line in lines if line]
        # Remove empty and comment lines
        lines = [line for line in lines if not line.startswith("#")]
        # Remove inline comments
        lines = [line.split("#", 1)[0].strip() for line in lines]
        iter_: Iterator[str] = iter(lines)
        # ATOM AND REFERENCE CONFIGURATION
        # atsym, z, nc, nv, iexc, psfile
        words: list[str] = next(iter_).split()
        atsym = words[0]
        z = float(words[1])
        nc = int(words[2])
        nv = int(words[3])
        iexc = int(words[4])
        psfile = words[5]
        # n, l, f (nc + nv lines)
        na: list[int] = [0] * (nc + nv)
        la: list[int] = [0] * (nc + nv)
        fa: list[float] = [0.0] * (nc + nv)
        for i in range(nc + nv):
            words = next(iter_).split()
            na[i] = int(words[0])
            la[i] = int(words[1])
            fa[i] = float(words[2])
        # PSEUDOPOTENTIAL AND OPTIMIZATION
        # lmax
        lmax = int(next(iter_).split()[0])
        # l, rc, ep, ncon, nbas, qcut (lmax + 1 lines, l's must be in order)
        rc: list[float] = [0.0] * (lmax + 1)
        ep: list[float] = [0.0] * (lmax + 1)
        ncon: list[int] = [0] * (lmax + 1)
        nbas: list[int] = [0] * (lmax + 1)
        qcut: list[float] = [0.0] * (lmax + 1)
        for l in range(lmax + 1):
            words = next(iter_).split()
            assert l == int(words[0])
            rc[l] = float(words[1])
            ep[l] = float(words[2])
            ncon[l] = int(words[3])
            nbas[l] = int(words[4])
            qcut[l] = float(words[5])
        # LOCAL POTENTIAL
        # lloc, lpopt, rc(5), dvloc0
        words = next(iter_).split()
        lloc = int(words[0])
        lpopt = int(words[1])
        rc5 = float(words[2])
        dvloc0 = float(words[3])
        # VANDERBILT-KLEINMAN-BYLANDER PROJECTORs
        # nproj, debl (lmax + 1 lines, l's must be in order)
        nproj: list[int] = [0] * (lmax + 1)
        debl: list[float] = [0.0] * (lmax + 1)
        for l in range(lmax + 1):
            words = next(iter_).split()
            assert l == int(words[0])
            nproj[l] = int(words[1])
            debl[l] = float(words[2])
        # MODEL CORE CHARGE
        # icmod, fcfact, (rcfact)
        words = next(iter_).split()
        icmod = int(words[0])
        fcfact = float(words[1])
        rcfact: float = 0.0 if len(words) < 3 else float(words[2])
        # LOG DERIVATIVE ANALYSIS
        # epsh1, epsh2, depsh, (rxpsh)
        words = next(iter_).split()
        epsh1 = float(words[0])
        epsh2 = float(words[1])
        depsh = float(words[2])
        rxpsh = None if len(words) < 4 else float(words[3])
        # OUTPUT GRID
        # rlmax, drl
        words = next(iter_).split()
        rlmax = float(words[0])
        drl = float(words[1])
        # TEST CONFIGURATIONS
        # ncnf
        ncnf = int(next(iter_).strip())
        # nvcnf (repeated ncnf times)
        # n, l, f (nvcnf lines, repeated following ncnf's ncnf times)
        test_configurations: list[dict[str, list[int | float]]] = []
        nvcnf: list[int] = [0] * ncnf
        for i in range(ncnf):
            nvcnf[i] = int(next(iter_).strip())
            test_configurations.append({"n": [0] * nvcnf[i], "l": [0] * nvcnf[i], "f": [0.0] * nvcnf[i]})
            for j in range(nvcnf[i]):
                words = next(iter_).split()
                test_configurations[i]["n"][j] = int(words[0])
                test_configurations[i]["l"][j] = int(words[1])
                test_configurations[i]["f"][j] = float(words[2])
        data = {
            "oncvpsp": {
                "atsym": atsym,
                "z": z,
                "nc": nc,
                "nv": nv,
                "iexc": iexc,
            },
            "reference_configuration": {
                "n": na,
                "l": la,
                "f": fa,
            },
            "pseudopotentials": {
                "lmax": lmax,
                "l": list(range(lmax + 1)),
                "rc": rc,
                "ep": ep,
                "ncon": ncon,
                "nbas": nbas,
                "qcut": qcut,
            },
            "local_potential": {
                "lloc": lloc,
                "lpopt": lpopt,
                "rcloc": rc5,
                "dvloc0": dvloc0,
            },
            "vkb_projectors": {
                "l": list(range(lmax + 1)),
                "nproj": nproj,
                "debl": debl,
            },
            "model_core_charge": {
                "icmod": icmod,
                "fcfact": fcfact,
                "rcfact": rcfact,
            },
            "log_derivative_analysis": {
                "epsh1": epsh1,
                "epsh2": epsh2,
                "depsh": depsh,
                "rxpsh": rxpsh,
            },
            "linear_mesh": {
                "rlmax": rlmax,
                "drl": drl,
            },
            "pp_output": {
                "psfile": psfile,
            },
            "test_configurations": test_configurations,
        }
        return cls(**data)  # type: ignore

    def as_dict(self, flat=False) -> dict[str, Any]:
        """Return the ONCVPSP input parameters as a nested dictionary.

        Returns:
            dict[str, Any]: ONCVPSP input parameters.
        """
        dict_ = self.model_dump()
        if flat:
            dict_ = flatten(dict_)
        return dict_

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "OncvpspInput":
        """Create an OncvpspInput instance from a nested dictionary."""
        return cls(**data)

    def as_toml_string(self) -> str:
        """Return the ONCVPSP input parameters as a TOML string.

        Returns:
            str: ONCVPSP input parameters in TOML format.
        """
        data = self.model_dump()
        return toml.dumps(data)

    @classmethod
    def from_toml_string(cls, s: str) -> "OncvpspInput":
        """Create an OncvpspInput instance from a TOML string.

        Args:
            s (str): ONCVPSP input parameters in TOML format.
        Returns:
            OncvpspInput: ONCVPSP input parameters.
        """
        data = toml.loads(s)
        return cls(**data)

    def to_toml_file(self, file: str | PathLike | TextIO) -> None:
        """Write the ONCVPSP input parameters to a TOML file.

        Args:
            file (PathLike | TextIO): Path to TOML file or file-like object.
        """
        if isinstance(file, (str, PathLike)):
            with open(file, "w", encoding="utf-8") as fp:
                toml.dump(self.model_dump(), fp)
        else:
            toml.dump(self.model_dump(), file)

    @classmethod
    def from_toml_file(cls, file: str | PathLike | TextIO) -> "OncvpspInput":
        """Load the ONCVPSP input parameters from a TOML file.

        Args:
            file (PathLike | TextIO): Path to TOML file or file-like object.
        Returns:
            OncvpspInput: ONCVPSP input parameters.
        """
        if isinstance(file, (str, PathLike)):
            with open(file, "r", encoding="utf-8") as fp:
                data = toml.load(fp)
        else:
            data = toml.load(file)
        return cls(**data)

    @model_validator(mode="after")
    def _check_number_of_states(self) -> "OncvpspInput":
        if len(self.reference_configuration) != self.oncvpsp.nc + self.oncvpsp.nv:
            raise ValueError(
                f"Number of states in reference configuration ({len(self.reference_configuration)}) "
                f"does not match n_core + n_valence ({self.oncvpsp.nc} + {self.oncvpsp.nv} = "
                f"{self.oncvpsp.nc + self.oncvpsp.nv})"
            )
        if not 1 <= (self.oncvpsp.nc + self.oncvpsp.nv) <= 30:
            raise ValueError(
                f"Total number of states (n_core + n_valence = {self.oncvpsp.nc + self.oncvpsp.nv}) "
                f"must be between 1 and 30"
            )
        return self

    @model_validator(mode="after")
    def _check_reference_configuration(self) -> "OncvpspInput":
        for i, state in enumerate(self.reference_configuration):
            if state.l > self.pseudopotentials.lmax:
                raise ValueError(
                    f"State at index {i} with l={state.l} exceeds"
                    f" l_max={self.pseudopotentials.lmax} in reference configuration"
                )
        total_occupation = sum(state.f for state in self.reference_configuration)
        if total_occupation > self.oncvpsp.z:
            raise ValueError(
                "The reference configuration is a negative ion: "
                f"total occupation {total_occupation} > atomic charge {self.oncvpsp.z}"
            )
        return self

    @model_validator(mode="after")
    def _check_pseudo_potentials(self) -> "OncvpspInput":
        if self.local_potential.lloc == 4:
            for psp in self.pseudopotentials:
                if psp.rc < self.local_potential.rcloc:  # type: ignore
                    raise ValueError(
                        f"Pseudopotential with l={psp.l} has rc={psp.rc} "
                        f"smaller than local potential rc={self.local_potential.rcloc}"
                    )
        return self

    @model_validator(mode="after")
    def _check_output_grid(self) -> "OncvpspInput":
        rc_max: float = max(self.pseudopotentials.rc)
        if self.local_potential.lloc == 4:
            rc_max = max(rc_max, self.local_potential.rcloc)  # type: ignore
        if self.linear_mesh.rlmax < rc_max:
            raise ValueError(
                f"Output grid rlmax={self.linear_mesh.rlmax} is smaller than "
                f"maximum rc={rc_max} of pseudopotentials"
            )
        return self
