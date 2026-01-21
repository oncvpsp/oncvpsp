from collections import defaultdict
from functools import cached_property
from os import PathLike
from typing import TextIO, Any
import pathlib as pl
import re

import numpy as np

from pyoncvpsp.io._utils import fort_float
from pyoncvpsp.io.input._models import OncvpspInput

RE_FLOAT: str = r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][+-]?\d+)?"

WARNINGS = [
    {
        "name": "icmod1_not_converged",
        "pattern": r"WARNING - modcore not converged",
        "match_string": "WARNING - modcore not converged",
        "line_count": 1,
        "subroutine": "modcore",
        "description": "Optimization of polynomial model core charge did not converge",
    },
    {
        "name": "icmod4_nm_not_converged",
        "pattern": r"WARNING: not fully converged in 100 steps",
        "match_string": "WARNING: not fully converged in 100 steps",
        "line_count": 1,
        "subroutine": "modcore3",
        "description": "Nelder-Mead optimization of model core charge Teter parameters did not converge",
    },
    {
        "name": "atom_not_converged",
        "pattern": r"(?P<rel>sr|rel)atom: WARNING failed to converge",
        "match_string": "atom: WARNING failed to converge",
        "line_count": 1,
        "subroutine": "(sr|rel)atom",
        "description": "All-electron calculation did not converge",
    },
    {
        "name": "ldiracfb_not_converged",
        "pattern": r"runconfig: WARNING ldiracfb convergence ERROR n,l,kap,iter=",
        "match_string": "runconfig: WARNING ldiracfb convergence ERROR n,l,kap,iter=",
        "line_count": 1,
        "subroutine": "run_config_r",
        "description": "",
    },
    {
        "name": "fr_no_ae_solution",
        "pattern": r"run_config_r: WARNING  for AE atom,",
        "match_string": "run_config_r: WARNING  for AE atom,",
        "line_count": 1,
        "subroutine": "run_config_r",
        "description": "",
    },
    {
        "name": "fr_fully_non_local_ps_atom",
        "pattern": r"run_config_r: WARNING for fully non-local PS atom",
        "match_string": "run_config_r: WARNING for fully non-local PS atom",
        "line_count": 1,
        "subroutine": "run_config_r",
        "description": "",
    },
    {
        "name": "localized_bound_state_for_projector",
        "pattern": r"WARNING wellstate: localized bound state found for n=",
        "match_string": "WARNING wellstate: localized bound state found for n=",
        "line_count": 1,
        "subroutine": "wellstate[_r]",
        "description": "",
    },
    {
        "name": "moderately_localized_bound_state_for_projector",
        "pattern": r"WARNING wellstate: moderately localized bound state found for n=",
        "match_string": "WARNING wellstate: moderately localized bound state found for n=",
        "line_count": 1,
        "subroutine": "wellstate[_r]",
        "description": "",
    },
    {
        "name": "scattering_first_projector",
        "pattern": r"WARNING wellstate: negative energy specified for n=\s*(?P<enn>\d+) l=\s*(?P<ell>\d+)",
        "match_string": "WARNING wellstate: negative energy specified for n=",
        "line_count": 3,
        "subroutine": "wellstate",
        "description": "",
    },
    {
        "name": "lschvkbbe_not_converged",
        "pattern": rf"(?P<subroutine>\w+): lschvkbbe ERROR\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+({RE_FLOAT})\s+({RE_FLOAT})",
        "match_string": "lschvkbbe ERROR",
        "line_count": 1,
        "subroutine": "lschvkbbe",
        "description": "Unable to converge fixed-energy solution in local+non-local potential.",
    },
]

ERRORS = [
    {
        "name": "lschvkbb_not_converged",
        "pattern": r"psatom(?:_r)?: WARNING lschvkbb convergence error n,l,(?:kap,)?iter=",
        "match_string": "WARNING lschvkbb convergence error n,l,",
        "line_count": 1,
        "subroutine": "psatom(_r)",
        "description": "Unable to converge bound-state solution in local+non-local potential.",
    },
    {
        "name": "wellstate_not_converged",
        "pattern": r"ERROR wellstate(?:_r)?: well potential iteration failed to converge, n=",
        "match_string": "ERROR wellstate: well potential iteration failed to converge, n=",
        "line_count": 1,
        "subroutine": "wellstate(_r)",
        "description": "Unable to converge well-state solution in full potential.",
    },
    {
        "name": "bad_input",
        "pattern": r"ERROR: test_data found\s*(?P<nerr>\d+)\s*ERROR; stopping",
        "match_string": "ERROR: test_data found",
        "line_count": 1,
        "subroutine": "check_data",
        "description": "Invalid input data.",
    },
    {
        "name": "input_file_read_error",
        "pattern": r"Error reading input data file",
        "match_string": "Error reading input data file",
        "line_count": 1,
        "subroutine": "input",
        "description": "Unable to read input data file.",
    },
    {
        "name": "no_classical_turning_point",
        "pattern": r"ERROR no classical turning point",
        "match_string": "ERROR no classical turning point",
        "line_count": 1,
        "subroutine": "lsch*",
        "description": "No classical turning point found.",
    },
    {
        "name": "no_core_valence_crossing",
        "pattern": r"ERROR ircc (cor-valence charge crossover",
        "match_string": "ERROR ircc (cor-valence charge crossover",
        "line_count": 1,
        "subroutine": "modcore*",
        "description": "No core-valence charge density crossover radius found.",
    },
    {
        "name": "first_pswf_has_node",
        "pattern": r"ERROR pspot:  first pseudo wave function has node,",
        "match_string": "ERROR pspot:  first pseudo wave function has node,",
        "line_count": 2,
        "subroutine": "pspot",
        "description": "First pseudo wavefunction has a node.",
    },
    {
        "name": "ps0norm_gt_uunorm",
        "pattern": r"optimize: ERROR ps0norm > uunorm, code will stop",
        "match_string": "optimize: ERROR ps0norm > uunorm, code will stop",
        "line_count": 1,
        "subroutine": "optimize",
        "description": "Norm of unoptimized pseudo wavefunction exceeds that of all-electron wavefunction.",
    },
    {
        "name": "rxpsh_too_small",
        "pattern": rf"ERROR: run_phsft: rxpsh= (?P<rxpsh>{RE_FLOAT}) < rc\((?P<l>\d)\)= (?P<rc>{RE_FLOAT})",
        "match_string": "ERROR: run_phsft: rxpsh=",
        "line_count": 1,
        "subroutine": "run_phsft",
        "description": "Matching radius is smaller than pseudopotential core radius for given angular momentum.",
    }
]


def _kappa(ell: int, jay: float | None = None, ess: float | None = None, spin_sign: int | None = None) -> int:
    sources = iter(x is not None for x in [jay, ess, spin_sign])
    assert any(sources) and not any(sources), "Exactly one input source must be provided."
    if spin_sign is not None:
        ess: float = spin_sign * 1 / 2
    if ess is not None:
        jay: float = ell + ess
    return int((ell - jay) * (2 * jay + 1))


class OncvpspOutputError(Exception):

    line_numbers: list[int]
    context_lines: list[str]
    context_line_numbers: list[int]
    name: str
    subroutine: str
    description: str

    def __init__(
        self,
        line_numbers: list[int],
        context_lines: list[str],
        context_line_numbers: list[int],
        name: str,
        subroutine: str,
        description: str,
    ) -> None:
        self.line_numbers = line_numbers
        self.context_lines = context_lines
        self.context_line_numbers = context_line_numbers
        self.name = name
        self.subroutine = subroutine
        self.description = description
        message = (
            f"ONCVPSP error encountered at line {line_numbers[0]}:\n"
            f"  Subroutine: {subroutine}\n"
            f"  Description: {description}\n"
            "  Context:\n"
            + "\n".join(
                f"\t{ln:>6d} {'>' if ln in line_numbers else ' '}|{line}"
                for ln, line in zip(context_line_numbers, context_lines)
            )
        )
        super().__init__(message)

    @property
    def error_lines(self) -> list[str]:
        """Lines which contain the error.

        Returns:
            list[str]: Lines containing the error.
        """
        return [line for ln, line in zip(self.context_line_numbers, self.context_lines) if ln in self.line_numbers]


class OncvpspTextParser:

    _lines: list[str]

    def __init__(
        self,
        path: str | PathLike | None = None,
        io: TextIO | None = None,
        text: str | None = None,
        lines: list[str] | None = None,
        check: bool = True,
        werror: bool = False,
    ) -> None:
        sources = iter(x is not None for x in [path, io, text, lines])
        # First `any`: consume the iteratble until the first True or it is exhausted
        # Second `any`: check if there are any remaining True values in the iterator
        assert any(sources) and not any(sources), "Exactly one input source must be provided."
        if text is not None:
            lines = text.splitlines(keepends=False)
        elif path is not None:
            path = pl.Path(path)
            with path.open("r", encoding="latin-1") as f:
                lines = f.readlines()
        elif io is not None:
            lines = io.readlines()
        else:
            raise RuntimeError("unreachable")
        self._lines = [line.rstrip() for line in lines]
        if check:
            self.check(werror=werror)

    def _parse_errors(self, definitions: list[dict], context_window: int | tuple[int, int] = 3) -> list[dict]:
        _context_window = (context_window, context_window) if isinstance(context_window, int) else context_window
        errors = []
        for i, line in enumerate(self._lines):
            if line.strip().startswith("DATA FOR PLOTTING"):
                break
            line_number = i + 1
            for definition in definitions:
                if definition["match_string"] in line:
                    line_numbers = range(line_number, line_number + definition.get("line_count", 1))
                    context_start = max(0, i - _context_window[0])
                    context_end = min(len(self._lines), i + _context_window[1] + 1)
                    context_lines = self._lines[context_start:context_end]
                    context_line_numbers = list(range(context_start + 1, context_end + 1))
                    warning = {
                        "line_numbers": list(line_numbers),
                        "context_lines": context_lines,
                        "context_line_numbers": context_line_numbers,
                        "name": definition.get("name", "unknown_error"),
                        "subroutine": definition.get("subroutine", "Unknown"),
                        "description": definition.get("description", "No description available."),
                    }
                    errors.append(warning)
        return errors

    def _get_lines_between(
        self,
        start_pattern: str | re.Pattern[str],
        stop_pattern: str | re.Pattern[str],
        include_start: bool = False,
        include_stop: bool = False,
    ) -> list[str]:
        start_match = self.search(start_pattern)
        if not start_match:
            return []
        start_idx = start_match[0] + (0 if include_start else 1)
        stop_match = self.search(stop_pattern, start=start_idx)
        if not stop_match:
            return []
        end_idx = stop_match[0] + (1 if include_stop else 0)
        return self.lines[start_idx:end_idx]

    def search(self, pattern: str | re.Pattern[str], start: int = 0) -> tuple[int, re.Match[str]] | None:
        """Find the line in `self.lines` which matches the regex `pattern` starting from `self.lines[start]`.

        Args:
            pattern (str | re.Pattern[str]): Pattern to search for.
            start (int, optional): Index of line to start searching from. Defaults to 0.

        Returns:
            tuple[int, re.Match[str]] | None: Tuple of (line index, regex match object) if found, else None.
        """
        regex = re.compile(pattern) if isinstance(pattern, str) else pattern
        for i in range(start, len(self.lines)):
            if match := regex.search(self.lines[i]):
                return i, match
        return None

    def findall(
        self, pattern: str | re.Pattern[str], start: int = 0, stop: int = -1
    ) -> list[tuple[int, re.Match[str]]]:
        """Find all lines in `self.lines[start:stop]` which match the regex `pattern`.

        Args:
            pattern (str | re.Pattern[str]): Pattern to serch for.
            start (int, optional): Lower bound of search window (inclusive). Defaults to 0.
            stop (int, optional): Upper bound of search window (exclusive). Defaults to -1.

        Returns:
            list[tuple[int, re.Match[str]]]: List of tuples of (line index, regex match object).
        """
        regex = re.compile(pattern) if isinstance(pattern, str) else pattern
        stop = stop if stop >= 0 else len(self.lines)
        matches = []
        for i in range(start, stop):
            if match := regex.search(self.lines[i]):
                matches.append((i, match))
        return matches

    @cached_property
    def warnings(self) -> list[dict]:
        """List of warning summaries.

        Returns:
            list[dict]: List of warnings encountered in the output.
        """
        return self._parse_errors(WARNINGS)

    @cached_property
    def errors(self) -> list[dict]:
        """List of error summaries.

        Returns:
            list[dict]: List of errors encountered in the output.
        """
        return self._parse_errors(ERRORS)

    def check(self, doraise: bool = True, werror: bool = False) -> list[OncvpspOutputError]:
        """Check the output for errors and warnings (if `werror`).
        Raise an `ExceptionGroup` containing `OncvpspOutputError`s, if found.

        Args:
            werror (bool, optional): Raise errors for any warnings encountered. Defaults to False.

        Raises:
            ExceptionGroup: Group of `OncvpspOutputError`s describing errors and warnings encountered.
        """
        exceptions = []
        # Only give 'incomplete output' error when the program did not report any other errors
        # but is nonetheless unfinished (e.g. still running, killed, etc.)
        if not self.errors and not self.completed:
            exceptions.append(
                OncvpspOutputError(
                    line_numbers=[len(self._lines)],
                    context_lines=self._lines[-5:],
                    context_line_numbers=list(range(len(self._lines) - 4, len(self._lines) + 1)),
                    name="incomplete_output",
                    subroutine="main",
                    description="'DATA FOR PLOTTING' section not found.",
                )
            )
        exceptions.extend([OncvpspOutputError(**error) for error in self.errors])
        if werror:
            exceptions.extend([OncvpspOutputError(**warning) for warning in self.warnings])
        if exceptions and doraise:
            raise ExceptionGroup(
                f"ONCVPSP reported {len(self.errors)} errors during execution.",
                exceptions,
            )
        return exceptions

    @cached_property
    def completed(self) -> bool:
        """Whether ONCVPSP/METAPSP successfully pseudized the atom."""
        return self.search("DATA FOR PLOTTING") is not None

    @cached_property
    def lines(self) -> list[str]:
        """Lines of the output file, excluding those corresponding to warnings and errors."""
        # Remove lines corresponding to warnings and errors prior to parsing output quantities
        # because they may be inserted within formatted data blocks.
        warning_line_indices = [i - 1 for warning in self.warnings for i in warning["line_numbers"]]
        error_line_indices = [i - 1 for error in self.errors for i in error["line_numbers"]]
        skip_line_indices = set(warning_line_indices + error_line_indices)
        return [line for i, line in enumerate(self._lines) if i not in skip_line_indices]

    @cached_property
    def program_information(self) -> dict:
        """Description of the program (name, version, relativistic treatment, etc.)

        Raises:
            ValueError: If the program name is unknown.

        Returns:
            dict:
                program (str): Name of the program (ONCVPSP or METAPSP).
                version (str): Version of the program.
                date (str): Date of the program run.
                relativistic_treatment (str): Relativistic treatment used.
        """
        program: str = self.lines[0].strip().split()[0]
        if program == "ONCVPSP":
            line = self.lines[1]
            relativistic_treatment = line.strip().split()[0]
        elif program == "METAPSP":
            line = self.lines[2]
            relativistic_treatment = "non-relativistic"
        else:
            raise ValueError(f"Unknown program: {program}")
        version_match: re.Match[str] | None = re.search(r"(\d+\.\d+\.\d+)", line)
        version: str = version_match.group(1) if version_match else "0.0.0"
        date_match: re.Match[str] | None = re.search(r"(\d{2}/\d{2}/\d{4})", line)
        date: str = date_match.group(1) if date_match else "00/00/0000"
        return {
            "program": program,
            "version": version,
            "date": date,
            "relativistic_treatment": relativistic_treatment,
        }

    @cached_property
    def is_relativistic(self) -> bool:
        """Whether the fully-relativistic Dirac equation was solved."""
        return self.program_information["relativistic_treatment"] in ("relativistic", "fully-relativistic")

    @cached_property
    def input(self) -> dict:
        """Parsed input data used for the pseudopotential generation."""
        if line_match := self.search(r"^# ATOM AND REFERENCE CONFIGURATION$"):
            i, _ = line_match
            return OncvpspInput.from_dat_string("\n".join(self.lines[i:])).model_dump()
        return {}

    def get_reference_quantum_numbers(self, index) -> list[tuple[int, int, int]]:
        """Get the quantum numbers (n, l, sign(s)) describing the states of the reference configuration.
        `index` acts upon the reference configuration _as it is written in the input file_, i.e. before
        expanding for spin-orbit splitting in the relativistic case.

        Args:
            index (int | slice): State ind(ex|ices)

        Returns:
            list[tuple[int, int, int]]: Quantum number tuples (n, l, sign(s)).
        """
        config = self.input["reference_configuration"]
        quantum_numbers = []
        for enn, ell in zip(config["n"][index], config["l"][index]):
            spin_signs = [0] if not self.is_relativistic else ([1] if ell == 0 else [-1, 1])
            for spin_sign in spin_signs:
                quantum_numbers.append((enn, ell, spin_sign))
        return quantum_numbers

    @property
    def core_quantum_numbers(self) -> list[tuple[int, int, int]]:
        """Quantum number tuples (n, l, sign(s)) corresponding to the core states.

        Returns:
            list[tuple[int, int, int]]: Quantum number tuples (n, l, sign(s)).
        """
        n_core = self.input["oncvpsp"]["nc"]
        return self.get_reference_quantum_numbers(index=slice(0, n_core))

    @property
    def valence_quantum_numbers(self) -> list[tuple[int, int, int]]:
        """Quantum number tuples (n, l, sign(s)) corresponding to the valence states.

        Returns:
            list[tuple[int, int, int]]: Quantum number tuples (n, l, sign(s))
        """
        n_core = self.input["oncvpsp"]["nc"]
        return self.get_reference_quantum_numbers(index=slice(n_core, None))

    @cached_property
    def reference_configuration(self) -> list[dict]:
        """Information about the reference configuration states, including all-electron eigenvalues.

        Raises:
            RuntimeError: If the "ATOM AND REFERENCE CONFIGURATION" block cannot be found.

        Returns:
            dict:
                (l, n, spin_sign) (tuple[int, int, int]): State identifier tuple.
                    n (int): Principal quantum number.
                    l (int): Angular momentum quantum number.
                    f (float): Occupation number.
                    spin_sign (int): Spin sign (+1 for j=l+1/2, -1 for j=l-1/2, 0 for non-relativistic).
                    eig_ae (float): All-electron eigenvalue (Ha).
                    eig_pbe (float, optional): PBE eigenvalue (Ha), if available.
        """
        nlf_pattern = rf"\s+(?P<n>\d+)\s+(?P<l>\d+)\s+(?P<f>{RE_FLOAT})"
        pattern = re.compile(nlf_pattern + rf"\s+(?P<eig>{RE_FLOAT})(?:\s+(?P<eig_dw>{RE_FLOAT}))?")
        n_states = self.input["oncvpsp"]["nc"] + self.input["oncvpsp"]["nv"]
        if self.program_information["program"] == "METAPSP":
            header_match = self.search(r"#   n    l    f      MGGA eval \(Ha\)         PBE         delta")
            pattern = re.compile(nlf_pattern + rf"\s+(?P<eig>{RE_FLOAT})\s+(?P<eig_pbe>{RE_FLOAT})")
        elif self.program_information["version"].startswith("3") and self.is_relativistic:
            header_match = self.search(r"#   n    l    f        l\+1/2             l-1/2")
        else:
            header_match = self.search(r"#   n    l    f        energy \(Ha\)")
        if not header_match:
            raise RuntimeError("Failed to find reference configuration section in output.")
        data: list[dict] = []
        for _, match in self.findall(pattern, start=header_match[0] + 1, stop=header_match[0] + 1 + n_states):
            enn = int(match.group("n"))
            ell = int(match.group("l"))
            occ = float(match.group("f"))
            if not self.is_relativistic:
                eigs = {0: fort_float(match.group("eig"))}
            else:
                eigs = {1: fort_float(match.group("eig"))}
                if ell > 0:
                    eigs[-1] = fort_float(match.group("eig_dw"))
            for spin_sign, eig in eigs.items():
                data.append(
                    {
                        "n": enn,
                        "l": ell,
                        "f": occ,
                        "spin_sign": spin_sign,
                        "eig_ae": eig,
                    }
                )
                if self.program_information["program"] == "METAPSP":
                    data[-1]["eig_pbe"] = fort_float(match.group("eig_pbe"))
        return data

    @cached_property
    def wellstate_metadata(self) -> list[dict]:
        """Information about any wellstate solutions found in the output.

        Returns:
            dict:
                n (int): Principal quantum number.
                l (int): Angular momentum quantum number.
                iproj (int): Projector index (n - l for ONCVPSP v3, else 0).
                k (int): Dirac quantum number kappa (0 for non-relativistic).
                spin_sign (int): Spin sign (+1 for j=l+1/2, -1 for j=l-1/2, 0 for non-relativistic).
                eig_ae (float): All-electron eigenvalue (Ha).
                asymptotic_potential (float): Asymptotic potential (Ha).
                half_point_radius (float): Half-point radius (a.u.).
        """
        wellstate_metadata: list[dict] = []
        matches = self.findall(
            r"\s+Wellstate for l =\s+(?P<ell>\d)  n =\s+(?P<enn>\d+)(?:\s+kap(?:pa)?=\s+(?P<kappa>-?\d+))?"
        )
        eig_pattern = re.compile(rf"\s+eigenvalue =\s*(?P<eig>{RE_FLOAT})")
        ap_pattern = re.compile(rf"\s+asymptotic potential =\s*(?P<ap>{RE_FLOAT})")
        hpr_pattern = re.compile(rf"\s+half-point radius =\s*(?P<hpr>{RE_FLOAT})")
        for i, match in matches:
            eig_match = re.match(eig_pattern, self.lines[i + 1])
            ap_match = re.match(ap_pattern, self.lines[i + 2])
            hpr_match = re.match(hpr_pattern, self.lines[i + 3])
            enn = int(match.group("enn"))
            ell = int(match.group("ell"))
            kappa = int(match.group("kappa")) if match.group("kappa") is not None else 0
            iproj = enn - ell if self.program_information["version"].startswith("3") else 0
            spin_sign = 0 if kappa == 0 else (-1 if kappa > 0 else 1)
            wellstate_metadata.append(
                {
                    "n": enn,
                    "l": ell,
                    "iproj": iproj,
                    "k": kappa,
                    "spin_sign": spin_sign,
                    "eig_ae": float(eig_match.group("eig")) if eig_match else None,
                    "asymptotic_potential": (float(ap_match.group("ap")) if ap_match else None),
                    "half_point_radius": (float(hpr_match.group("hpr")) if hpr_match else None),
                }
            )
        return wellstate_metadata

    @cached_property
    def _optimize_line_indices(self) -> list[tuple[int, int, int, int]]:
        indices: dict = {}
        # ONCVPSP v3
        for start, match in self.findall(
            r"Calculating (?P<proj_name>first|second) optimized projector for l=\s+(?P<ell>\d+)"
        ):
            ell = int(match.group("ell"))
            iproj = 1 if match.group("proj_name") == "first" else 2
            if (ell, iproj, -1) in indices:
                spin_sign = +1
            elif self.is_relativistic and ell == 0:
                spin_sign = 1
            elif self.is_relativistic and ell > 0:
                spin_sign = -1
            else:
                spin_sign = 0
            indices[(ell, iproj, spin_sign)] = start
        # ONCVPSP v4 / METAPSP
        for start, match in self.findall(r"Calculating optimized projector #\s+(?P<iproj>\d+)"):
            iproj = int(match.group("iproj"))
            if ell_match := self.search(r" for l=\s+(?P<ell>\d+)", start=start):
                ell = int(ell_match[1].group("ell"))
            else:
                raise RuntimeError(f"Failed to find ang. mom. in `optimize` block starting at line {start+1}")
            if (ell, iproj, -1) in indices:
                spin_sign = +1
            elif self.is_relativistic and ell == 0:
                spin_sign = 1
            elif self.is_relativistic and ell > 0:
                spin_sign = -1
            else:
                spin_sign = 0
            indices[(ell, iproj, spin_sign)] = start
        return [(ell, iproj, spin_sign, start) for (ell, iproj, spin_sign), start in indices.items()]

    @property
    def _optimize_convergence_profiles(self) -> list[dict]:
        convergence_header_pattern = re.compile(r"\s+Ha\s+eV\s+Ha")
        convergence_body_pattern = re.compile(rf"\s+(?P<eresid>{RE_FLOAT})\s+{RE_FLOAT}\s+(?P<ecut>{RE_FLOAT})")
        data: dict[tuple, dict[str, list]] = defaultdict(lambda: defaultdict(list))
        for ell, iproj, spin_sign, start in self._optimize_line_indices:
            convergence_match = self.search(convergence_header_pattern, start=start)
            if convergence_match:
                for _, match in self.findall(
                    convergence_body_pattern,
                    start=convergence_match[0],
                    stop=convergence_match[0] + 5,
                ):
                    data[(ell, iproj, spin_sign)]["eresid"].append(float(match.group("eresid")))
                    data[(ell, iproj, spin_sign)]["ecut"].append(float(match.group("ecut")))
        return [
            {**{"l": key[0], "iproj": key[1], "spin_sign": key[2]}, **{k: np.array(v) for k, v in value.items()}}
            for key, value in data.items()
        ]

    @cached_property
    def _run_vkb_indices(self) -> list[tuple[int, int, int, int]]:
        b_mat_pattern = re.compile(
            r"B matrix(?: Hermiticity error)?, ll=\s+(?P<ell>\d+)(?:,)?(?:\s+kap=\s+(?P<kappa>-?\d+))?"
        )
        coeff_pattern = re.compile(rf"\s+Orthonormal projector coefficients(?:\s+{RE_FLOAT})+")
        blocks: list[tuple[int, int, int, int]] = []
        for start, start_match in self.findall(b_mat_pattern):
            if (_ := self.search(coeff_pattern, start=start + 1)) is None:
                raise RuntimeError(f"Failed to find end of VKB block (starting at line {start+1})")
            else:
                stop, _ = _
            ell = int(start_match.group("ell"))
            kappa = int(start_match.group("kappa")) if start_match.group("kappa") is not None else None
            spin_sign = 0 if kappa is None else (-1 if kappa > 0 else 1)
            blocks.append((ell, spin_sign, start, stop))
        return blocks

    @cached_property
    def vkb_projector_coefficients(self) -> list[dict]:
        """Orthonormal Vanderbilt-Kleinman-Bylander coefficients `Dii`.

        Returns:
            list[dict]:
                ell (int): Angular momentum quantum number.
                iproj (int): Projector index.
                spin_sign (int): Spin sign (+1 for j=l+1/2, -1 for j=l-1/2, 0 for non-relativistic).
                coefficient (float): VKB coefficient Dii.
        """
        coefficients: list[dict] = []
        for ell, spin_sign, _, stop in self._run_vkb_indices:
            for iproj, value in enumerate(self.lines[stop].split()[3:], start=1):
                coefficients.append({"l": ell, "iproj": iproj, "spin_sign": spin_sign, "coefficient": value})
        return coefficients

    @cached_property
    def vkb_scalar_projector_coefficients(self) -> list[dict]:
        coefficients: list[dict] = []
        for i, match in self.findall(r" Orthonormal scalar projector coefficients, l =\s+(?P<ell>\d)"):
            ell = int(match.group("ell"))
            coefficients.append(
                {
                    "l": ell,
                    "coefficients": np.array([float(c) for c in self.lines[i + 1].split()]),
                }
            )
        return coefficients

    @cached_property
    def vkb_spin_orbit_projector_coefficients(self) -> list[dict]:
        coefficients: list[dict] = []
        for i, match in self.findall(r" Orthonormal spin-orbit projector coefficients, l =\s+(?P<ell>\d)"):
            ell = int(match.group("ell"))
            coefficients.append(
                {
                    "l": ell,
                    "coefficients": np.array([float(c) for c in self.lines[i + 1].split()]),
                }
            )
        return coefficients

    @cached_property
    def vkb_hermiticity_errors(self) -> list[dict]:
        v3_herm_err_pattern = re.compile(rf"Hermiticity error\s+(?P<herm_err>{RE_FLOAT})")
        # TODO: v3_bij_pattern: ^B11, B12, B22   -3.8877D-03    1.0547D-02   -2.8717D-02$
        v4_herm_err_pattern = re.compile(rf"\s+(?P<i>\d)\s+(?P<j>\d)\s+(?P<herm_err>{RE_FLOAT})")
        hermiticity_errors: list[dict] = []
        for ell, spin_sign, start, stop in self._run_vkb_indices:
            # Parse hermiticity errors
            if herm_err_match := re.match(v3_herm_err_pattern, self.lines[start + 1]):
                hermiticity_errors.append(
                    {
                        "l": ell,
                        "spin_sign": spin_sign,
                        "i": 1,
                        "j": 2,
                        "hermiticity_error": fort_float(herm_err_match.group("herm_err")),
                    }
                )
            else:
                herm_err_matches = self.findall(v4_herm_err_pattern, start=start, stop=stop)
                for _, herm_err_match in herm_err_matches:
                    hermiticity_errors.append(
                        {
                            "l": ell,
                            "spin_sign": spin_sign,
                            "i": int(herm_err_match.group("i")),
                            "j": int(herm_err_match.group("j")),
                            "hermiticity_error": fort_float(herm_err_match.group("herm_err")),
                        }
                    )
        return hermiticity_errors

    @cached_property
    def d2exc(self) -> dict:
        """Second derivatives of the exchange-correlation energy w.r.t. valence state occupations.

        Returns:
            dict: AE/PS d2Exc/dfi with/without non-linear core corrections.
        """
        pattern = re.compile(
            r"d2exc(?P<ae_ps>ae|ps) - (?:all-electron|pseudofunction) derivatives with(?P<no> no)? core correction"
        )
        matches = self.findall(pattern)
        data = defaultdict(list)
        for i, match in matches:
            ae_ps = match.group("ae_ps")
            with_nlcc = match.group("no") is None
            key = f"{ae_ps}_with_nlcc" if with_nlcc else ae_ps
            for j in range(self.input["oncvpsp"]["nv"]):
                line = self.lines[i + 2 + j]
                data[key].append([fort_float(x) for x in line.split()])
        return {k: np.array(v) for k, v in data.items()}

    @cached_property
    def d2exc_rmse(self) -> float | None:
        """Root-mean-square error of PS d2Exc/dfi with non-linear core corrections w.r.t. AE d2Exc/dfi."""
        pattern = re.compile(r"rms 2nd derivative error\s+(?P<rmse>{RE_FLOAT})")
        match = self.search(pattern)
        return float(match[1].group("rmse")) if match else None

    @cached_property
    def teter_parameters(self) -> dict:
        """Parameters of the Teter-function non-linear core correction model core charge density.

        Teter function is defined as: A(sin(2πx) / [(2π(x))(1 - 4(x)^2)(1 - (x)^2)])^2
        where x = r / (rcfact * rmatch) and A = fcfact * rhocmatch.

        N.B.: if icmod=2, all of the quantities have _different_ meanings!

        Returns:
            dict:
                fcfact (float): Amplitude prefactor.
                rcfact (float): Scale prefactor.
                rmatch (float): Core-valence charge density crossover radius (a.u.).
                rhocmatch (float): Core charge density at crossover radius (a.u.^-3).
        """
        parameters: dict[str, float] = {}
        if self.input["model_core_charge"]["icmod"] == 3:
            parameters["fcfact"] = self.input["model_core_charge"]["fcfact"]
            parameters["rcfact"] = self.input["model_core_charge"]["rcfact"]
        if self.input["model_core_charge"]["icmod"] == 4:
            param_match = self.search(r"amplitude prefactor, scale prefactor")
            if param_match:
                parameters["fcfact"] = float(self.lines[param_match[0] + 1].split()[0])
                parameters["rcfact"] = float(self.lines[param_match[0] + 1].split()[1])
        match_match = self.search(r"rmatch, rhocmatch\s+(?P<rmatch>{RE_FLOAT})\s+(?P<rhocmatch>{RE_FLOAT})")
        if match_match:
            parameters["rmatch"] = float(match_match[1].group("rmatch"))
            parameters["rhocmatch"] = float(match_match[1].group("rhocmatch"))
        return parameters

    @cached_property
    def teter_coarse_grid(self) -> dict:
        """Coarse grid evaulation of d2Exc/dfi RMSE over different values of Teter prefactors.
        Present only when icmod=4 in ONCVPSP.

        Returns:
            dict:
                fcfact (np.ndarray): Array of fcfact values.
                rcfact (np.ndarray): Array of rcfact values.
                d2exc_rmse (np.ndarray): 2D array of d2Exc/dfi RMSE values.
        """
        match = self.search(r"Coarse scan for minimum error")
        if not match:
            return {}
        grid: dict[str, list] = defaultdict(list)
        grid["fcfact"] = [float(x) for x in self.lines[match[0] + 5].split()]
        i = match[0] + 7
        while words := self.lines[i].strip().split():
            grid["rcfact"].append(float(words[0]))
            grid["d2exc_rmse"].append([float(x) for x in words[1:]])
            i += 1
        return {k: np.array(v) for k, v in grid.items()}

    @cached_property
    def diagnostic_tests_semilocal(self) -> list[dict]:
        """Diagnostic test information when solving for states in the semilocal pseudopotential.

        Returns:
            list[dict]:
                l (int): Angular momentum quantum number.
                k (int): Dirac quantum number kappa (0 for non-relativistic).
                spin_sign (int): Spin sign (+1 for j=l+1/2, -1 for j=l-1/2, 0 for non-relativistic).
                rcore (float): Core radius (a.u.).
                rmatch (float): Matching radius (a.u.).
                e_in (float): Input all-electron eigenvalue (Ha).
                e_test (float): Output pseudoeigenvalue (Ha).
                norm_test (float): Norm w/in rc of pseudo wavefunction / AE wavefunction.
                slope_test (float): Slope at rc of pseudo wavefuction / AE wavefunction.
        """
        start_match = self.search(r"Diagnostic tests using semi-local pseudopotentials")
        stop_match = self.search(
            r"Diagnostic tests using Vanderbilt-Kleinman-Bylander pseudopotentials",
            start=start_match[0] if start_match else 0,
        )
        if not start_match or not stop_match:
            return []
        pattern = re.compile(
            rf"\s+(?P<ell>\d)(?:\s+(?P<kap>-?\d))?\s+(?P<rcore>{RE_FLOAT})\s+(?P<rmatch>{RE_FLOAT})\s+(?P<e_in>{RE_FLOAT})\s+(?P<e_test>{RE_FLOAT})\s+(?P<norm_test>{RE_FLOAT})\s+(?P<slope_test>{RE_FLOAT})"
        )
        return [
            {
                "l": int(match.group("ell")),
                "k": int(match.group("kap")) if match.group("kap") is not None else 0,
                "spin_sign": 0 if match.group("kap") is None else (-1 if int(match.group("kap")) > 0 else 1),
                "rc": float(match.group("rcore")),
                "rmatch": float(match.group("rmatch")),
                "e_in": float(match.group("e_in")),
                "e_test": float(match.group("e_test")),
                "norm_test": float(match.group("norm_test")),
                "slope_test": float(match.group("slope_test")),
            }
            for _, match in self.findall(pattern, start=start_match[0], stop=stop_match[0])
        ]

    @cached_property
    def diagnostic_tests_vkb(self) -> list[dict]:
        """Diagnostic test information when solving for states in the VKB pseudopotential.

        Returns:
            list[dict]:
                l (int): Angular momentum quantum number.
                k (int): Dirac quantum number kappa (0 for non-relativistic).
                spin_sign (int): Spin sign (+1 for j=l+1/2, -1 for j=l-1/2, 0 for non-relativistic).
                rcore (float): Core radius (a.u.).
                rmatch (float): Matching radius (a.u.).
                e_in (float): Input all-electron eigenvalue (Ha).
                e_test (float | None): Output pseudoeigenvalue (Ha), if ONCVPSP v3.
                delta_e (float | None): Eigenvalue error (Ha), if ONCVPSP v4 or METAPSP v1.
                norm_test (float): Norm w/in rc of pseudo wavefunction / AE wavefunction.
                slope_test (float): Slope at rc of pseudo wavefuction / AE wavefunction.
        """
        start_match = self.search(r"Diagnostic tests using Vanderbilt-Kleinman-Bylander pseudopotentials")
        stop_match = self.search(r"\s*Test", start=start_match[0] if start_match else 0)
        if not start_match:
            raise RuntimeError("Failed to find start of VKB diagnostic test block")
        if not stop_match:
            raise RuntimeError("Failed to find end of VKB diagnostic test block")
        is_v3 = self.program_information["version"].startswith("3")
        pattern = re.compile(
            rf"\s+(?P<ell>\d)(?:\s+(?P<kap>-?\d))?\s+(?P<rcore>{RE_FLOAT})\s+(?P<rmatch>{RE_FLOAT})\s+(?P<e_in>{RE_FLOAT})\s+(?P<e_or_delta_e_test>{RE_FLOAT})\s+(?P<norm_test>{RE_FLOAT})\s+(?P<slope_test>{RE_FLOAT})"
        )
        return [
            {
                "l": int(match.group("ell")),
                "k": int(match.group("kap")) if match.group("kap") is not None else 0,
                "spin_sign": 0 if match.group("kap") is None else (-1 if int(match.group("kap")) > 0 else 1),
                "rc": float(match.group("rcore")),
                "rmatch": float(match.group("rmatch")),
                "e_in": float(match.group("e_in")),
                "e_test": float(match.group("e_or_delta_e_test")) if is_v3 else None,
                "delta_e": float(match.group("e_or_delta_e_test")) if not is_v3 else None,
                "norm_test": float(match.group("norm_test")),
                "slope_test": float(match.group("slope_test")),
            }
            for _, match in self.findall(pattern, start=start_match[0], stop=stop_match[0])
        ]

    @cached_property
    def ghosts(self) -> list[dict]:
        """Information about potential ghost states.

        Returns:
            list[dict]:
                l (int): Angular momentum quantum number.
                k (int): Dirac quantum number kappa (0 for non-relativistic).
                spin_sign (int): Spin sign (+1 for j=l+1/2, -1 for j=l-1/2, 0 for non-relativistic).
                eig (float): Eigenvalue of the ghost state (Ha).
                ecut (float): Cutoff energy where the ghost appears (Ha).
                sign (str): Sign of the ghost state ('+' or '-').
                r_rc (float | None): r / rc, if unbound (sign == '+').
        """
        ghosts = []
        pattern = re.compile(
            rf"\s+(?P<ell>\d)(?:\s+(?P<kap>-?\d))?\s+(?:(?P<r_rc>{RE_FLOAT})\s+)?(?P<eig>{RE_FLOAT})\s+(?P<ecut>{RE_FLOAT})\s+WARNING - GHOST\((?P<sign>[+-])\)"
        )
        matches = self.findall(pattern)
        for _, match in matches:
            ghosts.append(
                {
                    "l": int(match.group("ell")),
                    "k": int(match.group("kap")) if match.group("kap") is not None else 0,
                    "spin_sign": 0 if match.group("kap") is None else (-1 if int(match.group("kap")) > 0 else 1),
                    "eig": float(match.group("eig")),
                    "ecut": float(match.group("ecut")),
                    "sign": match.group("sign"),
                    "r_rc": float(match.group("r_rc")) if match.group("r_rc") is not None else None,
                }
            )
        return ghosts

    @cached_property
    def test_configurations(self) -> list[list[dict]]:
        """Test configuration results.

        Returns:
            list[list[dict]]:
                n (int): Principal quantum number.
                l (int): Angular momentum quantum number.
                k (int): Dirac quantum number kappa (0 for non-relativistic).
                spin_sign (int): Spin sign (+1 for j=l+1/2, -1 for j=l-1/2, 0 for non-relativistic).
                f (float): Occupation number.
                eig_ae (float): All-electron eigenvalue (Ha).
                eig_ps (float | None): Pseudoeigenvalue (Ha), if available.
        """
        test_configurations: list[list[dict]] = []
        pattern = re.compile(
            rf"\s+(?P<enn>\d)\s+(?P<ell>\d)\s+(?:(?P<kappa>-?\d)\s+)?(?P<occ>{RE_FLOAT})\s+(?P<eig_ae>{RE_FLOAT})(?:\s+(?P<eig_ps>{RE_FLOAT}))?(?:\s+(?P<diff>{RE_FLOAT}))?"
        )
        test_matches = self.findall(r"Test configuration\s+(?P<itest>\d+)")
        for i, match in test_matches:
            test_lines = self._get_lines_between(
                start_pattern=self.lines[i],
                stop_pattern=rf"PSP excitation error=\s+{RE_FLOAT}",
                include_start=True,
                include_stop=True,
            )
            test_configuration: list[dict] = []
            for line in test_lines:
                if config_match := re.match(pattern, line):
                    enn = int(config_match.group("enn"))
                    ell = int(config_match.group("ell"))
                    kappa = int(config_match.group("kappa")) if config_match.group("kappa") is not None else 0
                    spin_sign = 0 if kappa == 0 else (-1 if kappa > 0 else 1)
                    test_configuration.append(
                        {
                            "n": enn,
                            "l": ell,
                            "k": kappa,
                            "spin_sign": spin_sign,
                            "f": float(config_match.group("occ")),
                            "eig_ae": float(config_match.group("eig_ae")),
                            "eig_ps": (
                                float(config_match.group("eig_ps"))
                                if config_match.group("eig_ps") is not None
                                else None
                            ),
                        }
                    )
            assert int(match.group("itest")) == len(test_configurations)
            test_configurations.append(test_configuration)
        return test_configurations

    def _get_plot_blocks(self, pattern: re.Pattern[str]) -> list[list[tuple[int, re.Match[str]]]]:
        matches = self.findall(pattern)
        blocks = []
        block = []
        i_prev = matches[0][0] - 1 if matches else 0
        for i, match in matches:
            if i - i_prev == 1:
                block.append((i, match))
                i_prev = i
            else:
                blocks.append(block)
                block = [(i, match)]
                i_prev = i
        if block:
            blocks.append(block)
        return blocks

    @cached_property
    def unscreened_semilocal_potentials(self) -> list[dict]:
        r"""Unscreened semilocal pseudopotential radial functions.

        Returns:
            list[dict]:
                l (int): Angular momentum quantum number.
                r (np.ndarray): Radial grid points (a.u.).
                v_sl (np.ndarray): Unscreened semilocal potential values (Ha).
        """
        pattern = re.compile(
            rf"!p\s+(?P<r>{RE_FLOAT})\s+(?P<rho_val>{RE_FLOAT})"
            + "".join(rf"\s+(?P<v_sl_{l}>{RE_FLOAT})" for l in range(self.input["pseudopotentials"]["lmax"] + 1))
        )
        block = self.findall(pattern)
        data: dict[int, dict[str, list]] = defaultdict(lambda: defaultdict(list))
        for _, match in block:
            r = float(match.group("r"))
            for ell in range(self.input["pseudopotentials"]["lmax"] + 1):
                data[ell]["r"].append(r)
                data[ell]["v_sl"].append(float(match.group(f"v_sl_{ell}")))
        return [{"l": ell, **{k: np.array(v) for k, v in values.items()}} for ell, values in data.items()]

    @cached_property
    def local_potential(self) -> dict:
        """Local potential radial function, if `lloc = 4`.

        Returns:
            dict:
                r (np.ndarray): Radial grid points (a.u.).
                v_loc (np.ndarray): Local potential values (Ha).
        """
        pattern = re.compile(rf" !L\s+(?P<r>{RE_FLOAT})\s+(?P<v_loc>{RE_FLOAT})")
        block = self.findall(pattern)
        data = defaultdict(list)
        for _, match in block:
            data["r"].append(float(match.group("r")))
            data["v_loc"].append(float(match.group("v_loc")))
        return {**{k: np.array(v) for k, v in data.items()}, **self.input.get("local_potential", {})}

    @cached_property
    def charge_densities(self) -> dict:
        r"""Radial charge densities.

        N.B.: defined as int[0,∞] r²ρi(r) dr = 1 for a single-particle density ρi(r).

        Returns:
            dict:
                r (np.ndarray): Radial grid points (a.u.).
                rho_val (np.ndarray): Valence charge density (e/a.u.³).
                rho_core (np.ndarray): Core charge density (e/a.u.³).
                rho_model_core (np.ndarray): Model core charge density (e/a.u.³).
        """
        pattern = re.compile(
            rf"!r\s+(?P<r>{RE_FLOAT})\s+(?P<rho_val>{RE_FLOAT})\s+(?P<rho_core>{RE_FLOAT})\s+(?P<rho_model_core>{RE_FLOAT})"
        )
        block = self.findall(pattern)
        data = defaultdict(list)
        for _, match in block:
            data["r"].append(float(match.group("r")))
            data["rho_val"].append(float(match.group("rho_val")))
            data["rho_core"].append(float(match.group("rho_core")))
            data["rho_model_core"].append(float(match.group("rho_model_core")))
        return {k: np.array(v) for k, v in data.items()}

    @cached_property
    def kinetic_energy_densities(self) -> dict:
        r"""Radial kinetic energy densities.

        Returns:
            dict:
                r (np.ndarray): Radial grid points (a.u.).
                tau_ps (np.ndarray): Pseudodensity kinetic energy density.
                tau_mod (np.ndarray): Model core kinetic energy density.
        """
        pattern = re.compile(rf"!t\s+(?P<r>{RE_FLOAT})\s+(?P<tau_ps>{RE_FLOAT})\s+(?P<tau_mod>{RE_FLOAT})")
        block = self.findall(pattern)
        data = defaultdict(list)
        for _, match in block:
            data["r"].append(float(match.group("r")))
            data["tau_ps"].append(float(match.group("tau_ps")))
            data["tau_mod"].append(float(match.group("tau_mod")))
        return {k: np.array(v) for k, v in data.items()}

    @cached_property
    def kinetic_energy_potentials(self) -> dict:
        r"""Radial kinetic energy potentials.

        Returns:
            dict:
                r (np.ndarray): Radial grid points (a.u.).
                vtau_ps (np.ndarray): Pseudodensity kinetic energy potential.
                vtau_mod (np.ndarray): Model core kinetic energy potential.
        """
        pattern = re.compile(rf"!vt\s+(?P<r>{RE_FLOAT})\s+(?P<vtau_ps>{RE_FLOAT})\s+(?P<vtau_mod>{RE_FLOAT})")
        block = self.findall(pattern)
        data = defaultdict(list)
        for _, match in block:
            data["r"].append(float(match.group("r")))
            data["vtau_ps"].append(float(match.group("vtau_ps")))
            data["vtau_mod"].append(float(match.group("vtau_mod")))
        return {k: np.array(v) for k, v in data.items()}

    @cached_property
    def wavefunctions(self) -> list[dict]:
        r"""Radial wavefunctions.

        Returns:
            list[dict]:
                r (np.ndarray): Radial grid points (a.u.).
                wfn_ae (np.ndarray): All-electron wavefunction values.
                wfn_ps (np.ndarray): Pseudowavefunction values.
                is_bound (bool | None): Whether the state is a bound state (True), a wellstate (False), or unknown (None).
                n (int | None): Principal quantum number, None if unknown.
                l (int): Angular momentum quantum number.
                iproj (int): Projector index (1-based).
                k (int): Dirac quantum number kappa (0 for non-relativistic).
                spin_sign (int): Spin sign (+1 for j=l+1/2, -1 for j=l-1/2, 0 for non-relativistic).
                f (float, optional): Occupation number, if available.
                eig_ae (float, optional): All-electron eigenvalue (Ha), if available.
                eig_ps (float, optional): Pseudoeigenvalue (Ha), if available.
                asymptotic_potential (float, optional): Asymptotic potential (Ha), if available.
                half_point_radius (float, optional): Half-point radius (a.u.), if available.
        """
        data_pattern = re.compile(
            rf"&\s+(?P<spin_sign>[\s-])(?P<iproj>\d)?(?P<ell>\d)\s+(?P<r>{RE_FLOAT})\s+(?P<wfn_ae>{RE_FLOAT})\s+(?P<wfn_ps>{RE_FLOAT})"
        )
        header_pattern = re.compile(r"n=\s*(?P<enn>\d+),\s+l=\s*(?P<ell>\d)")
        v4_scattering_header_pattern = re.compile(r"scattering, iprj=\s*(?P<iproj>\d+),\s+l=\s*(?P<ell>\d)")
        blocks = self._get_plot_blocks(data_pattern)
        # Used to augment state data with metadata from  reference configuration
        test_config = (
            {(cfg["n"], cfg["l"], cfg["spin_sign"]): cfg for cfg in self.test_configurations[0]}
            if self.test_configurations
            else {}
        )
        ref_config = {(cfg["n"], cfg["l"], cfg["spin_sign"]): cfg for cfg in self.reference_configuration}
        wellstates = {(well["n"], well["l"], well["spin_sign"]): well for well in self.wellstate_metadata}
        data: list[dict] = []
        for block in blocks:
            block_data: dict[str, Any] = defaultdict(list)
            # Parse the block data
            for _, match in block:
                block_data["r"].append(float(match.group("r")))
                block_data["wfn_ae"].append(float(match.group("wfn_ae")))
                block_data["wfn_ps"].append(float(match.group("wfn_ps")))
            block_data = {k: np.array(v) for k, v in block_data.items()}
            # Parse the header line and first block line for enn, ell, bound vs. scattering state, iproj
            first_ln, first_match = block[0]
            header_line: str = self.lines[first_ln - 2]
            block_data["is_bound"] = None
            if header_match := re.match(header_pattern, header_line):
                if self.program_information["version"].startswith("4"):
                    block_data["is_bound"] = True
                block_data["n"] = int(header_match.group("enn"))
                ell = int(header_match.group("ell"))
                if first_match.group("iproj"):
                    iproj = int(first_match.group("iproj"))
                elif int(first_match.group("ell")) == ell:
                    iproj = 1
                elif int(first_match.group("ell")) == ell + 4:
                    iproj = 2
                else:
                    raise RuntimeError(
                        f"Failed to infer projector index for ONCVPSP v3 from ell={ell} in header line: '{header_line}'"
                    )
            elif header_match := re.match(v4_scattering_header_pattern, header_line):
                block_data["is_bound"] = False
                ell = int(header_match.group("ell"))
                iproj = int(header_match.group("iproj"))
                block_data["n"] = None
            else:
                raise RuntimeError(f"Failed to parse wavefunction header line: '{header_line}'")
            # Determine the sign of the spin quantum number
            if self.is_relativistic:
                spin_sign = 1 if first_match.group("spin_sign") == "-" or ell == 0 else -1
                jay = ell + spin_sign * 1 / 2
                block_data["k"] = int((ell - jay) * (2 * jay + 1))
            else:
                spin_sign = 0
                block_data["k"] = 0
            # Add additional information about the state from the wellstate metadata or reference configuration
            # Order the state keys at this ell and kappa: first boundstates (sorted by enn), then wellstates (ditto)
            state_keys = sorted(
                [k for k in self.valence_quantum_numbers if k[1] == ell and k[2] == spin_sign]
            ) + sorted([k for k in wellstates if k[1] == ell and k[2] == spin_sign])
            state_key = state_keys[iproj - 1] if (iproj - 1 < len(state_keys)) else None
            state_metadata = {
                "is_bound": None,
                "n": 0,
                "f": np.nan,
                "eig_ae": np.nan,
                "eig_ps": np.nan,
                "asymptotic_potential": np.nan,
                "half_point_radius": np.nan,
                "r_rc": np.nan,
            }  # Default empty metadata
            if not block_data["is_bound"] is False and state_key in self.valence_quantum_numbers:
                # If we're sure this is not a wellstate, get a boundstate
                state_metadata.update(test_config.get(state_key, {}))
                state_metadata.update(ref_config[state_key])
                state_metadata["is_bound"] = True
            if not block_data["is_bound"] is True and state_key in wellstates:
                # If we're sure this is not a boundstate, get a wellstate
                state_metadata.update(wellstates[state_key])
                state_metadata["is_bound"] = False
                state_metadata["f"] = 0.0
            block_data.update(state_metadata)
            # Finalize block data
            block_data["l"] = ell
            block_data["iproj"] = iproj
            block_data["spin_sign"] = spin_sign
            data.append(block_data)
        return data

    @cached_property
    def vkb_projectors(self) -> list[dict]:
        """Vanderbilt-Kleinman-Bylander projectors.

        Returns:
            list[dict]:
                r (np.ndarray): Radial grid points (a.u.).
                proj (np.ndarray): VKB projector values.
                coeff (float | None): VKB coefficient Dii, if available.
                l (int): Angular momentum quantum number.
                iproj (int): Projector index (1-based).
                spin_sign (int): Spin sign (+1 for j=l+1/2, -1 for j=l-1/2, 0 for non-relativistic).
        """
        max_nproj: int = max(self.input["vkb_projectors"]["nproj"])
        pattern = re.compile(
            rf"(?:@|!J)\s+(?P<spin_sign>[\s-])(?P<ell>\d)\s+(?P<r>{RE_FLOAT})\s+(?P<proj_1>{RE_FLOAT})"
            + "".join(rf"(?:\s+(?P<proj_{i}>{RE_FLOAT}))?" for i in range(2, max_nproj + 1))
        )
        # Used to look up VKB coefficients
        coeffs = {(c["l"], c["iproj"], c["spin_sign"]): c["coefficient"] for c in self.vkb_projector_coefficients}
        blocks = self._get_plot_blocks(pattern)
        data: list[dict] = []
        for block in blocks:
            _, first_match = block[0]
            ell: int = int(first_match.group("ell"))
            spin_sign: int = 0
            if self.is_relativistic:
                spin_sign = 1 if first_match.group("spin_sign") == "-" or ell == 0 else -1
            block_data: dict[str, Any] = defaultdict(list)
            for _, match in block:
                block_data["r"].append(float(match.group("r")))
                for i in range(1, max_nproj + 1):
                    proj_key = f"proj_{i}"
                    if match.group(proj_key) is not None:
                        block_data[proj_key].append(float(match.group(proj_key)))
            block_data = {k: np.array(v) for k, v in block_data.items()}
            for iproj in range(1, max_nproj + 1):
                proj_key = f"proj_{iproj}"
                if proj_key in block_data:
                    jay = ell + spin_sign * 1 / 2
                    kappa = int((ell - jay) * (2 * jay + 1))
                    data.append(
                        {
                            "l": ell,
                            "iproj": iproj,
                            "k": kappa,
                            "spin_sign": spin_sign,
                            "r": block_data["r"],
                            "proj": block_data[proj_key],
                            "coeff": coeffs.get((ell, iproj, spin_sign), None),
                        }
                    )
        return data

    @property
    def _plot_convergence_profiles(self) -> list[dict]:
        pattern = re.compile(rf"!C\s+(?P<ell>\d)\s+(?P<ecut>{RE_FLOAT})\s+(?P<eresid>{RE_FLOAT})")
        block = self.findall(pattern)
        data: dict[tuple[int, int, int], dict[str, list]] = defaultdict(lambda: defaultdict(list))
        spin_sign: int = 1 if self.is_relativistic else 0
        for _, match in block:
            ell = int(match.group("ell"))
            data[(ell, 1, spin_sign)]["ecut"].append(float(match.group("ecut")))
            data[(ell, 1, spin_sign)]["eresid"].append(float(match.group("eresid")))
        return [
            {**{"l": key[0], "iproj": key[1], "spin_sign": key[2]}, **{k: np.array(v) for k, v in value.items()}}
            for key, value in data.items()
        ]

    @cached_property
    def convergence_profiles(self) -> list[dict]:
        """Convergence behavior of the VKB projectors w.r.t. Ecut.

        Returns:
            list[dict]:
                l (int): Angular momentum quantum number.
                iproj (int): Projector index (1-based).
                spin_sign (int): Spin sign (+1 for j=l+1/2, -1 for j=l-1/2, 0 for non-relativistic).
                ecut (np.ndarray): Cutoff energies (Ha).
                eresid (np.ndarray): Kinetic energy residual L1 errors of the VKB projectors.
        """
        return self._optimize_convergence_profiles + self._plot_convergence_profiles

    @cached_property
    def ecut_recommended(self) -> float:
        ecut = -np.inf
        for prof in self.convergence_profiles:
            if prof["iproj"] != 1:
                continue
            mask = np.isclose(prof["eresid"], 1e-5, atol=1e-10)
            if np.sum(mask) == 0 or not np.any(mask):
                continue
            ecut_proj = prof["ecut"][mask][0]
            if ecut_proj > ecut:
                ecut = ecut_proj
        return ecut

    @cached_property
    def log_derivatives(self) -> list[dict]:
        """Phase-shift-like AE/PS logarithmic derivatives.

        Returns:
            list[dict]:
                l (int): Angular momentum quantum number.
                k (int): Dirac quantum number kappa (0 for non-relativistic).
                spin_sign (int): Spin sign (+1 for j=l+1/2, -1 for j=l-1/2, 0 for non-relativistic).
                rxpsh (float): Radius for logarithmic derivative evaluation (a.u.).
                e (np.ndarray): Energies (Ha).
                log_deriv_ae (np.ndarray): All-electron logarithmic derivative.
                log_deriv_ps (np.ndarray): Pseudopotential logarithmic derivative.
        """
        pattern = re.compile(
            rf"!\s+(?P<spin_sign>[\s-])(?P<ell>\d)\s+(?P<e>{RE_FLOAT})\s+(?P<log_deriv_ae>{RE_FLOAT})\s+(?P<log_deriv_ps>{RE_FLOAT})"
        )
        blocks = self._get_plot_blocks(pattern)
        data: dict[tuple[int, int], dict[str, Any]] = defaultdict(lambda: defaultdict(list))
        for block in blocks:
            first_ln, first_match = block[0]
            header_line: str = self.lines[first_ln - 3]
            ell = int(first_match.group("ell"))
            spin_sign: int = 0
            if self.is_relativistic:
                spin_sign = 1 if first_match.group("spin_sign") == "-" or ell == 0 else -1
            data[(ell, spin_sign)]["rxpsh"] = float(header_line.split("=")[-1])
            for _, match in block:
                ell = int(match.group("ell"))
                if self.is_relativistic:
                    spin_sign = 1 if first_match.group("spin_sign") == "-" or ell == 0 else -1
                data[(ell, spin_sign)]["e"].append(float(match.group("e")))
                data[(ell, spin_sign)]["log_deriv_ae"].append(float(match.group("log_deriv_ae")))
                data[(ell, spin_sign)]["log_deriv_ps"].append(float(match.group("log_deriv_ps")))
        return [
            {
                **{"l": key[0], "spin_sign": key[1], "k": _kappa(ell=key[0], spin_sign=key[1])},
                **{k: np.array(v) for k, v in value.items()},
            }
            for key, value in data.items()
        ]

    @cached_property
    def plot_data(self) -> dict:
        """Dictionary containing plot data (radial functions).

        Returns:
            dict:
                local_potential (dict): See `OncvpspTextParser.local_potential`.
                charge_densities (dict): See `OncvpspTextParser.charge_densities`.
                kinetic_energy_densities (dict): See `OncvpspTextParser.kinetic_energy_densities`.
                kinetic_energy_potentials (dict): See `OncvpspTextParser.kinetic_energy_potentials`.
                wavefunctions (dict): See `OncvpspTextParser.wavefunctions`.
                vkb_projectors (dict): See `OncvpspTextParser.vkb_projectors`.
                convergence_profiles (dict): See `OncvpspTextParser.convergence_profiles`.
                log_derivatives (dict): See `OncvpspTextParser.log_derivatives`.
        """
        return {
            "local_potential": self.local_potential,
            "charge_densities": self.charge_densities,
            "kinetic_energy_densities": self.kinetic_energy_densities,
            "kinetic_energy_potentials": self.kinetic_energy_potentials,
            "wavefunctions": self.wavefunctions,
            "vkb_projectors": self.vkb_projectors,
            "convergence_profiles": self.convergence_profiles,
            "log_derivatives": self.log_derivatives,
        }

    @cached_property
    def plot_string(self) -> str:
        """String containing data for plotting."""
        return "\n".join(self._get_lines_between("DATA FOR PLOTTING", "GNUSCRIPT"))

    @cached_property
    def gnuplot_script(self) -> str:
        """String containing gnuplot script for plotting."""
        return "\n".join(self._get_lines_between("GNUSCRIPT", "END_GNU"))

    @cached_property
    def upf_string(self) -> str:
        """String containing the pseudopotential in UPF format."""
        return "\n".join(
            self._get_lines_between('<UPF version="2.0.1">', "</UPF>", include_start=True, include_stop=True)
        )

    @cached_property
    def psp8_string(self) -> str:
        """String containing the pseudopotential in Abinit PSP8 format."""
        return "\n".join(self._get_lines_between("Begin PSPCODE8", "</INPUT>"))

    @cached_property
    def output(self) -> dict:
        """Parsed output as a dictionary."""
        output: dict[str, Any] = {}
        output["completed"] = self.completed
        output["errors"] = self.errors
        output["warnings"] = self.warnings
        output["program_information"] = self.program_information
        output["input"] = self.input
        output["wellstate_metadata"] = self.wellstate_metadata
        output["vkb_hermiticity_errors"] = self.vkb_hermiticity_errors
        output["vkb_scalar_projector_coefficients"] = self.vkb_scalar_projector_coefficients
        output["vkb_spin_orbit_projector_coefficients"] = self.vkb_spin_orbit_projector_coefficients
        output["d2exc"] = self.d2exc
        output["d2exc_rmse"] = self.d2exc_rmse
        output["teter_parameters"] = self.teter_parameters
        output["teter_coarse_grid"] = self.teter_coarse_grid
        output["diagnostic_tests_semilocal"] = self.diagnostic_tests_semilocal
        output["diagnostic_tests_vkb"] = self.diagnostic_tests_vkb
        output["ghosts"] = self.ghosts
        output["test_configurations"] = self.test_configurations
        output["plot_data"] = self.plot_data
        output["gnuplot_script"] = self.gnuplot_script
        output["upf"] = self.upf_string
        output["psp8"] = self.psp8_string
        return output

    def to_dict(self) -> dict:
        """Parsed output as a serializable dictionary."""

        def serialize(value: Any) -> Any:
            if isinstance(value, np.ndarray):
                return value.tolist()
            if isinstance(value, dict):
                return {k: serialize(v) for k, v in value.items()}
            if isinstance(value, list):
                return [serialize(v) for v in value]
            return value

        return serialize(self.output)

    def parse(self, check: bool = True, werror: bool = False) -> dict:
        """Parse the output file.

        Args:
            check (bool, optional): Check and raise exceptions for errors. Defaults to True.
            werror (bool, optional): Consider warnings as errors. Defaults to False.

        Returns:
            dict: `self.output` dictionary
        """
        if check:
            self.check(werror=werror)
        return self.output
