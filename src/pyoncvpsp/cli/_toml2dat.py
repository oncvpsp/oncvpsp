"""CLI command to convert ONCVPSP TOML input files to dat format."""

from os import PathLike

import click

from pyoncvpsp.io import OncvpspInput


@click.command(name="toml2dat")
@click.option(
    "--input",
    "-i",
    "toml_filepath",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the ONCVPSP TOML file to convert.",
)
@click.option(
    "--output",
    "-o",
    "dat_filepath",
    type=click.Path(dir_okay=False),
    required=False,
    help="Path to save the converted dat file.",
)
def toml2dat(toml_filepath: PathLike, dat_filepath: PathLike | None) -> None:
    """Convert an ONCVPSP TOML input file to a dat input file.

    Args:
        toml_filepath (PathLike): Path to the ONCVPSP TOML input file.
        dat_filepath (PathLike): Path to save the converted dat input file.
    """
    oncvpsp_input: OncvpspInput = OncvpspInput.from_toml_file(toml_filepath)
    if dat_filepath is None:
        print(oncvpsp_input.as_dat_string())
    else:
        oncvpsp_input.to_dat_file(dat_filepath)
