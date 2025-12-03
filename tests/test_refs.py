"""Run tests comparing program output against reference outputs for various elements."""

import io
import pathlib
import subprocess

import pytest

from pyoncvpsp.io.output import OncvpspTextParser

from .nested_approx import nested_approx

TEST_NAMES: list[str] = [
    "03_Li",
    "07_N",
    "07_N_dl",
    "08_O",
    "14_Si",
    "14_Si_st",
    "14_Si_UPF",
    "17_Cl",
    "19_K",
    "19_K_st",
    "20_Ca",
    "22_Ti",
    "29_Cu",
    "32_Ge",
    "34_Se",
    "38_Sr_lxc",
    "40_Zr",
    "56_Ba",
    "57_La",
    "60_Nd_GHOST",
    "73_Ta",
    "79_Au_lxc",
    "83_Bi",
    "52_Te",
    "74_W",
    "80_Hg",
]


def _determine_executable(output_file: pathlib.Path, filepaths_executables: dict[str, pathlib.Path]) -> pathlib.Path:
    with open(output_file, "r", encoding="latin-1") as fp:
        output_first_line: str = fp.readline().strip()
        output_second_line: str = fp.readline().strip()
    program: str = output_first_line.split()[0]
    if program == "ONCVPSP":
        relativistic: str = output_second_line.split()[0]
        if relativistic == "non-relativistic":
            exe_name = "oncvpspnr.x"
        elif relativistic == "scalar-relativistic":
            exe_name = "oncvpsp.x"
        elif relativistic == "relativistic":
            exe_name = "oncvpspr.x"
        else:
            raise ValueError(f"Unknown relativistic option: {relativistic}")
    elif program == "METAPSP":
        exe_name = "oncvpspm.x"
    else:
        raise ValueError(f"Unknown program: {program}")

    try:
        exe_path: pathlib.Path = filepaths_executables[exe_name]
    except KeyError as exc:
        raise RuntimeError(f"Executable for {exe_name} not found.") from exc

    return exe_path


def _inequalities_message(inequalities: dict) -> str:
    messages: list[str] = []
    for path, ((val1, val2), reason) in inequalities.items():
        messages.append(f"Difference at {path}: {val1} vs {val2} ({reason})")
    return "\n".join(messages)


@pytest.mark.parametrize("name", TEST_NAMES)
def test_against_reference(
    name: str,
    filepaths_executables: dict[str, pathlib.Path],
    filepath_data: pathlib.Path,
    filepath_refs: pathlib.Path,
    tmp_path: pathlib.Path,
) -> None:
    """Integration test comparing program output against reference output."""
    # Find input and output files
    input_file: pathlib.Path = filepath_data / f"{name}.dat"
    if not input_file.exists():
        raise FileNotFoundError(f"Input file {input_file} does not exist.")
    # Output file could be either name.out (scalar-relativistic) or name_r.out (fully-relativistic)
    output_file: pathlib.Path = filepath_refs / f"{name}.out"
    if not output_file.exists():
        output_file = filepath_refs / f"{name}_r.out"
        if not output_file.exists():
            raise FileNotFoundError(f"Reference output file {filepath_refs}/{name}[_r].out does not exist.")
    exe_path: pathlib.Path = _determine_executable(output_file, filepaths_executables)
    workdir: pathlib.Path = tmp_path / name
    workdir.mkdir(exist_ok=True)
    with open(input_file, "r", encoding="latin-1") as fp:
        proc: subprocess.CompletedProcess[str] = subprocess.run(
            args=[str(exe_path)], stdin=fp, capture_output=True, text=True, check=False, cwd=str(workdir)
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"Execution of {exe_path.name} failed with return code {proc.returncode}."
                f"\nstdout:\n{proc.stdout}"
                f"\nstderr:\n{proc.stderr}"
            )
    # Parse reference and test outputs
    with open(output_file, "r", encoding="latin-1") as fp:
        ref_output = OncvpspTextParser(io=fp).parse()
    with io.StringIO(proc.stdout) as fp:
        tst_output = OncvpspTextParser(io=fp).parse()
    # Compare outputs recursively
    inequalities = nested_approx(
        ref_output,
        tst_output,
        abs=1e-5,
        rel=1e-12,
        exclude=[
            "program_information/date",
            "gnuplot_script",
            "psp8",
            "upf"
        ],
    )
    if inequalities:
        pytest.fail(f"Outputs differ\n{_inequalities_message(inequalities)}")
