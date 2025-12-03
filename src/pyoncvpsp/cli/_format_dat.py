"""CLI command to format ONCVPSP dat files."""

from os import PathLike

import click

from pyoncvpsp.io import OncvpspInput


@click.command(name="format-dat")
@click.option(
    "--input",
    "-i",
    "input_filepath",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the ONCVPSP dat file to convert.",
)
@click.option(
    "--output",
    "-o",
    "output_filepath",
    type=click.Path(dir_okay=False),
    required=False,
    help="Path to save the converted TOML file.",
)
def format_dat(input_filepath: PathLike, output_filepath: PathLike | None) -> None:
    """Format an ONCVPSP dat input file.

    Args:
        input_filepath (PathLike): Path to the ONCVPSP dat input file.
        output_filepath (PathLike): Path to save the formatted dat input file.
    """
    oncvpsp_input: OncvpspInput = OncvpspInput.from_dat_file(input_filepath)
    if output_filepath is None:
        print(oncvpsp_input.as_dat_string())
    else:
        oncvpsp_input.to_dat_file(output_filepath)
