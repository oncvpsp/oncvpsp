"""CLI command to convert ONCVPSP dat files to TOML format."""

from os import PathLike

import click

from pyoncvpsp.io import OncvpspInput


@click.command(name="dat2toml")
@click.option(
    "--input",
    "-i",
    "dat_filepath",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the ONCVPSP dat file to convert.",
)
@click.option(
    "--output",
    "-o",
    "toml_filepath",
    type=click.Path(dir_okay=False),
    required=False,
    help="Path to save the converted TOML file.",
)
def dat2toml(dat_filepath: PathLike, toml_filepath: PathLike | None) -> None:
    """Convert an ONCVPSP dat input file to a TOML input file.

    Args:
        dat_filepath (PathLike): Path to the ONCVPSP dat input file.
        toml_filepath (PathLike): Path to save the converted TOML input file.
    """
    oncvpsp_input: OncvpspInput = OncvpspInput.from_dat_file(dat_filepath)
    if toml_filepath is None:
        print(oncvpsp_input.as_toml_string())
    else:
        oncvpsp_input.to_toml_file(toml_filepath)
