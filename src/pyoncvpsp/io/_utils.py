"""I/O utilities"""

from collections import deque
from copy import deepcopy
from typing import Any, Mapping, Sequence

__all__ = ("fort_float",)


def fort_float(s: str) -> float:
    """Convert Fortran-style floating-point string to Python float.
    Replaces 'D' or 'd' with 'e' before conversion.

    Args:
        s: String to convert.

    Returns:
        Converted float.
    """
    return float(s.lower().replace("d", "e"))


def flatten(dictionary: Mapping[str, Any], separator="__") -> dict[str, Any]:
    """Flatten an arbitrarily nested dictionary by joining keys at nesting boundaries with `separator`.
    If lists are present at any level in the dictionary, they will be treated as if they were dictionaries
    with their list indices as keys.

    Args:
        dictionary (Mapping[str, Any]): Nested dictionary, potentially containing lists with arbitrary-typed elements.
        separator (str, optional): Separator used to join keys at nesting boundaries. Defaults to "__".

    Returns:
        dict[str, Any]: Flattened dictionary.
    """
    flat = {}
    queue = deque([(dictionary, "")])
    while queue:
        current_obj, parent_key = queue.popleft()
        if isinstance(current_obj, Mapping):
            generator = current_obj.items()
        elif isinstance(current_obj, Sequence):
            generator = enumerate(current_obj)
        else:
            continue
        for key, value in generator:
            new_key = f"{parent_key}{separator}{key}" if parent_key else str(key)
            if isinstance(value, (dict, list)):
                queue.append((value, new_key))  # type: ignore
            else:
                flat[new_key] = value
    return flat


def unflatten(
    dictionary: Mapping[str, Any], separator="__"
):  # pylint: disable=too-many-nested-blocks,too-many-branches
    """Unflatten a dictionary whose nestings have been reduced by joining keys with `separator`.
    Integer keys at any level will lead to the creation of a list rather than a dictionary at that level.

    Args:
        dictionary (Mapping[str, Any]): Flattened dictionary
        separator (str, optional): Separator used to join keys at nesting boundaries. Defaults to "__".

    Returns:
        dict[str, Any]: Reconstructed nested dictionary.
    """
    if all(key.split(separator)[0].isdigit() for key in dictionary.keys()):
        nested = []
    else:
        nested = {}
    for key, value in dictionary.items():
        parts = key.split(separator)
        d = nested
        for i, part in enumerate(parts[:-1]):
            if part.isdigit():
                part = int(part)
                if i == len(parts) - 1:
                    while len(d) <= part:
                        d.append(None)
                else:
                    if parts[i + 1].isdigit():
                        while len(d) <= part:
                            d.append([])
                    else:
                        while len(d) <= part:
                            d.append({})
            elif part not in d:
                if parts[i + 1].isdigit():
                    d[part] = []
                else:
                    d[part] = {}
            d = d[part]
        if isinstance(d, list):
            index = int(parts[-1])
            while len(d) <= index:
                d.append(None)
            d[index] = value
        else:
            d[parts[-1]] = value
    return nested


def recursive_merge(dict1: dict, dict2: dict) -> dict:
    """Recursively merge two dictionaries. Keys from dict2 take precedence.

    Args:
        dict1 (dict): First dictionary (overwritten)
        dict2 (dict): Second dictionary (overwriting)

    Returns:
        dict: Merged dictionary
    """
    dict1 = deepcopy(dict1)
    for key, value in dict2.items():
        if key in dict1 and isinstance(dict1[key], dict) and isinstance(value, dict):
            # Recursively merge nested dictionaries
            dict1[key] = recursive_merge(dict1[key], value)
        else:
            # Merge non-dictionary values
            dict1[key] = value
    return dict1
