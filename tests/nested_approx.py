"""Perform approximate equality checks on the leaves of a nested dictionary / list structure."""

from typing import Mapping, Sequence
import numpy as np

import pytest


def nested_approx(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    data1,
    data2,
    abs: float | None = None,  # pylint: disable=redefined-builtin
    rel: float | None = None,
    path: list | None = None,
    inequalities: dict | None = None,
    exclude: Sequence[str] | None = None,
) -> dict[str, tuple[tuple, str]]:
    """Recursively compare two nested dictionaries with approximate equality for numerical values.

    Args:
        data1: First object to compare.
        data2: Second object to compare.
        atol (float): Absolute tolerance.
        rtol (float): Relative tolerance.

    Returns:
        bool: True if the dictionaries are approximately equal, False otherwise.
    """
    if path is None:
        path = []

    if inequalities is None:
        inequalities = {}

    if exclude is None:
        exclude = []

    if isinstance(data1, Mapping) and isinstance(data2, Mapping):
        if data1.keys() != data2.keys():
            path_str = "/".join(path)
            if path_str not in exclude:
                inequalities["/".join(path)] = (
                    (data1.keys(), data2.keys()),
                    "keys mismatch",
                )
            return inequalities
        for key in data1:
            inequalities = nested_approx(
                data1[key],
                data2[key],
                abs=abs,
                rel=rel,
                path=path + [str(key)],
                inequalities=inequalities,
                exclude=exclude,
            )
        return inequalities

    if isinstance(data1, np.ndarray) and isinstance(data2, np.ndarray):
        if data1.shape != data2.shape:
            path_str = "/".join(path)
            if path_str not in exclude:
                inequalities["/".join(path)] = ((data1.shape, data2.shape), "shape mismatch")
            return inequalities
        if not np.allclose(data1, data2, atol=abs, rtol=rel, equal_nan=True):
            path_str = "/".join(path)
            if path_str not in exclude:
                inequalities["/".join(path)] = ((data1, data2), "ndarray mismatch")
        return inequalities

    if isinstance(data1, (list, tuple)) and isinstance(data2, (list, tuple)):
        if len(data1) != len(data2):
            path_str = "/".join(path)
            if path_str not in exclude:
                inequalities["/".join(path)] = ((len(data1), len(data2)), "length mismatch")
            return inequalities
        for key, (item1, item2) in enumerate(zip(data1, data2)):
            inequalities = nested_approx(
                item1,
                item2,
                abs=abs,
                rel=rel,
                path=path + [str(key)],
                inequalities=inequalities,
                exclude=exclude,
            )
        return inequalities

    if isinstance(data1, float) and isinstance(data2, float):
        if not data1 == pytest.approx(data2, abs=abs, rel=rel):
            path_str = "/".join(path)
            if path_str not in exclude:
                inequalities["/".join(path)] = ((data1, data2), "float mismatch")
        return inequalities

    if not data1 == data2:
        path_str = "/".join(path)
        if path_str not in exclude:
            inequalities["/".join(path)] = ((data1, data2), "value mismatch")

    return inequalities
