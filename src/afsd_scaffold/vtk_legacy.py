"""Minimal deterministic ASCII legacy VTK point-data I/O.

It covers only the subset used by the independent geometry fixture.  Physics
outputs may later use VTU/HDF5, but must preserve the same logical field names.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence


Number = int | float
Vector3 = tuple[float, float, float]


class LegacyVtkError(ValueError):
    """Raised when fixture VTK data violates the supported subset."""


def write_point_cloud(
    path: str | Path,
    points: Sequence[Vector3],
    scalar_fields: Mapping[str, Sequence[Number]],
    vector_fields: Mapping[str, Sequence[Vector3]],
    *,
    title: str,
) -> None:
    target = Path(path)
    count = len(points)
    for name, values in scalar_fields.items():
        if len(values) != count:
            raise LegacyVtkError(f"scalar {name!r} length does not match points")
    for name, values in vector_fields.items():
        if len(values) != count:
            raise LegacyVtkError(f"vector {name!r} length does not match points")

    lines = [
        "# vtk DataFile Version 3.0",
        title,
        "ASCII",
        "DATASET POLYDATA",
        f"POINTS {count} double",
    ]
    lines.extend("{:.17g} {:.17g} {:.17g}".format(*point) for point in points)
    lines.append(f"VERTICES {count} {2 * count}")
    lines.extend(f"1 {index}" for index in range(count))
    lines.append(f"POINT_DATA {count}")
    for name, values in scalar_fields.items():
        vtk_type = "long" if all(isinstance(value, int) for value in values) else "double"
        lines.append(f"SCALARS {name} {vtk_type} 1")
        lines.append("LOOKUP_TABLE default")
        lines.extend(str(value) if isinstance(value, int) else f"{value:.17g}" for value in values)
    for name, values in vector_fields.items():
        lines.append(f"VECTORS {name} double")
        lines.extend("{:.17g} {:.17g} {:.17g}".format(*value) for value in values)
    target.write_text("\n".join(lines) + "\n", encoding="ascii", newline="\n")


def read_point_cloud(path: str | Path) -> dict[str, object]:
    """Read files produced by :func:`write_point_cloud`."""

    lines = Path(path).read_text(encoding="ascii").splitlines()
    if len(lines) < 6 or lines[0] != "# vtk DataFile Version 3.0":
        raise LegacyVtkError("not a supported legacy VTK file")
    if lines[2] != "ASCII" or lines[3] != "DATASET POLYDATA":
        raise LegacyVtkError("only ASCII POLYDATA is supported")
    point_header = lines[4].split()
    if len(point_header) != 3 or point_header[0] != "POINTS":
        raise LegacyVtkError("POINTS section missing")
    count = int(point_header[1])
    cursor = 5
    points: list[Vector3] = []
    for _ in range(count):
        values = tuple(float(item) for item in lines[cursor].split())
        if len(values) != 3:
            raise LegacyVtkError("point must have three components")
        points.append(values)  # type: ignore[arg-type]
        cursor += 1

    vertex_header = lines[cursor].split()
    if vertex_header[:1] != ["VERTICES"] or int(vertex_header[1]) != count:
        raise LegacyVtkError("VERTICES section missing or inconsistent")
    cursor += 1 + count
    point_data_header = lines[cursor].split()
    if point_data_header != ["POINT_DATA", str(count)]:
        raise LegacyVtkError("POINT_DATA section missing or inconsistent")
    cursor += 1

    scalars: dict[str, list[Number]] = {}
    vectors: dict[str, list[Vector3]] = {}
    while cursor < len(lines):
        header = lines[cursor].split()
        if header[0] == "SCALARS":
            if len(header) != 4 or header[3] != "1":
                raise LegacyVtkError("only one-component SCALARS are supported")
            name, vtk_type = header[1], header[2]
            cursor += 1
            if lines[cursor] != "LOOKUP_TABLE default":
                raise LegacyVtkError("SCALARS lookup table is missing")
            cursor += 1
            if vtk_type == "long":
                scalars[name] = [int(lines[cursor + offset]) for offset in range(count)]
            elif vtk_type == "double":
                scalars[name] = [float(lines[cursor + offset]) for offset in range(count)]
            else:
                raise LegacyVtkError(f"unsupported scalar type {vtk_type!r}")
            cursor += count
        elif header[0] == "VECTORS":
            if len(header) != 3 or header[2] != "double":
                raise LegacyVtkError("only double VECTORS are supported")
            name = header[1]
            cursor += 1
            values: list[Vector3] = []
            for _ in range(count):
                vector = tuple(float(item) for item in lines[cursor].split())
                if len(vector) != 3:
                    raise LegacyVtkError("vector must have three components")
                values.append(vector)  # type: ignore[arg-type]
                cursor += 1
            vectors[name] = values
        else:
            raise LegacyVtkError(f"unsupported section {header[0]!r}")
    return {
        "title": lines[1],
        "points": points,
        "scalars": scalars,
        "vectors": vectors,
    }
