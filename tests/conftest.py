"""Test configuration for ONCVPSP / METAPSP tests."""

import os
import pathlib
import subprocess
import typing

import pytest

EXECUTABLE_NAMES: typing.Final[list[str]] = [
    "oncvpsp.x",
    "oncvpspnr.x",
    "oncvpspr.x",
]


def _filepath_tests() -> pathlib.Path:
    return pathlib.Path(__file__).parent.resolve()


@pytest.fixture
def filepath_tests() -> pathlib.Path:
    """Absolute path of "oncvpsp/tests/".

    Returns:
        pathlib.Path: Absolute path of "oncvpsp/tests/".
    """
    return _filepath_tests()


@pytest.fixture
def filepath_data(
    filepath_tests: pathlib.Path,  # pylint: disable=redefined-outer-name
) -> pathlib.Path:
    """Absolute path of "oncvpsp/tests/data/".

    Args:
        filepath_tests (pathlib.Path): Absolute path of "oncvpsp/tests/".

    Returns:
        pathlib.Path: Absolute path of "oncvpsp/tests/data/".
    """
    return filepath_tests / "data"


@pytest.fixture
def filepath_refs(
    filepath_tests: pathlib.Path,  # pylint: disable=redefined-outer-name
) -> pathlib.Path:
    """Absolute path of "oncvpsp/tests/refs/".

    Args:
        filepath_tests (pathlib.Path): Absolute path of "oncvpsp/tests/".

    Returns:
        pathlib.Path: _description_
    """
    return filepath_tests / "refs"


def _filepaths_executables(
    filepath_tests: pathlib.Path,
) -> dict[str, pathlib.Path]:
    executable_paths: dict[str, pathlib.Path] = {}
    # Search in PATH
    for exe_name in EXECUTABLE_NAMES:
        proc: subprocess.CompletedProcess[str] = subprocess.run(
            ["which", exe_name], capture_output=True, text=True, check=False
        )
        which_result: str = proc.stdout.strip()
        if which_result:
            executable_paths[exe_name] = pathlib.Path(which_result).resolve()
    # Search the directory in which `make` places the binaries (prioritize over PATH)
    filepath_make_bin: pathlib.Path = filepath_tests.parent / "src"
    for exe_name in EXECUTABLE_NAMES:
        exe_path: pathlib.Path = filepath_make_bin / exe_name
        if exe_path.exists():
            executable_paths[exe_name] = exe_path.resolve()
    # Search the suggested CMake build directory (prioritize over make)
    filepath_cmake_bin: pathlib.Path = filepath_tests.parent / "build" / "bin"
    for exe_name in EXECUTABLE_NAMES:
        exe_path: pathlib.Path = filepath_cmake_bin / exe_name
        if exe_path.exists():
            executable_paths[exe_name] = exe_path.resolve()
    # Search ONCVPSP_BIN_DIR environment variable (prioritize over CMake)
    env_bin_dir: str | None = os.getenv("ONCVPSP_BIN_DIR")
    if env_bin_dir is not None:
        filepath_env_bin: pathlib.Path = pathlib.Path(env_bin_dir)
        for exe_name in EXECUTABLE_NAMES:
            exe_path: pathlib.Path = filepath_env_bin / exe_name
            if exe_path.exists():
                executable_paths[exe_name] = exe_path.resolve()
    return executable_paths


@pytest.fixture
def filepaths_executables(
    filepath_tests: pathlib.Path,  # pylint: disable=redefined-outer-name
) -> dict[str, pathlib.Path]:
    """Absolute paths to the ONCVPSP executables.

    Searches the following sources (later sources have priority over earlier ones):
        1. The system PATH via `which`.
        2. "oncvpsp/src/" (the directory in which `make` places the binaries).
        3. "oncvpsp/build/src/" (the suggested CMake build directory).
        4. ONCVPSP_BIN_DIR environment variable (if set).

    If an executable is not found, it will not be included in the returned dictionary.

    Args:
        filepath_tests (pathlib.Path): Absolute path of "oncvpsp/tests/".

    Returns:
        dict[str, pathlib.Path]: Dictionary mapping executable names to their absolute paths.
    """
    return _filepaths_executables(filepath_tests)


@pytest.hookimpl()
def pytest_sessionstart(session) -> None:  # pylint: disable=unused-argument
    """Hook to print a message at the start of the test session."""
    executable_paths: dict[str, pathlib.Path] = _filepaths_executables(_filepath_tests())
    print("ONCVPSP executables:")
    for exe_name in EXECUTABLE_NAMES:
        if exe_name in executable_paths:
            print(f"  {exe_name:>12s}: {executable_paths[exe_name]}")
        else:
            print(f"  {exe_name:>12s}: NOT FOUND")
