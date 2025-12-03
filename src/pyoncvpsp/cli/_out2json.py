"""CLI command to convert ONCVPSP stdout files to JSON format."""

from os import PathLike
import json
from gzip import open as gzip_open

import click

from pyoncvpsp.io import OncvpspTextParser

@click.command(name="out2json")
@click.option(
    "--input",
    "-i",
    "out_filepath",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the ONCVPSP stdout file to convert.",
)
@click.option(
    "--output",
    "-o",
    "json_filepath",
    type=click.Path(dir_okay=False),
    required=False,
    help="Path to save the converted JSON file.",
)
@click.option(
    "--pretty",
    is_flag=True,
    default=False,
    help="Pretty-print the JSON output.",
)
@click.option(
    "--gzip",
    is_flag=True,
    default=False,
    help="gz-compress the JSON output.",
)
def out2json(out_filepath: PathLike, json_filepath: PathLike | None, pretty: bool, gzip: bool) -> None:
    """Convert an ONCVPSP stdout file to a JSON file.

    Args:
        out_filepath (PathLike): Path to the ONCVPSP stdout file.
        json_filepath (PathLike): Path to save the converted JSON file.
    """
    parser: OncvpspTextParser = OncvpspTextParser(path=out_filepath)
    output: dict = parser.to_dict()
    indent = 2 if pretty and not gzip else None
    if json_filepath is None:
        print(json.dumps(output, indent=indent))
    else:
        if gzip:
            with gzip_open(json_filepath, "wt", encoding="utf-8") as fp:
                json.dump(output, fp, indent=indent)
        else:
            with open(json_filepath, "w", encoding="utf-8") as fp:
                json.dump(output, fp, indent=indent)
