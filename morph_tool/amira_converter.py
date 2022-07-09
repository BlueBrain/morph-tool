"""Amira morphology format converter."""
from collections import defaultdict
import numpy as np
import pandas as pd

from morphio import PointLevel, SectionType
from morphio.mut import Morphology


def get_section_data(lines):
    """Read file and extract data per section into a dict of lists."""
    _converter = {
        "@1": "node_points",
        "@2": "node_label",
        "@3": "edge",
        "@4": "n_points",
        "@5": "label",
        "@6": "point",
        "@7": "radius",
    }
    sections = defaultdict(list)
    current_section = ""
    for line in lines:
        if "@" in line:
            current_section = line
        elif not line:
            current_section = ""
        elif current_section:
            sections[_converter[current_section]].append(line)
    return sections


def get_labels(lines, level=0):
    """Get label mapping."""
    label = None
    line_id = 0
    retry = False
    _id = None
    labels = {}
    while line_id <= len(lines):
        line = lines[line_id]

        if "{" in line and not retry:
            label = line.strip().split(" ")[0]

        if "{" in line and "{" in lines[line_id + 1] and not retry:
            _labels, n_lines = get_labels(lines[line_id + 1 :], level + 1)
            labels.update(_labels)
            line_id += n_lines
            retry = True
        elif retry and "{" in line:
            _labels, n_lines = get_labels(lines[line_id:], level + 1)
            labels.update(_labels)
            line_id += n_lines - 1
            retry = True
        if "Id" in line:
            if line.endswith(","):
                line = line[:-1]
            _id = int(line.split(" ")[-1])

        if "}" in line:
            if _id is not None:
                labels[_id] = label
            return labels, line_id + 1
        line_id += 1


def create_dfs(sections, labels):
    """Create dataframe with morphology data."""
    point_id = 0
    df = pd.DataFrame(columns=["u", "v", "label", "points", "diameters"])

    us = []
    vs = []
    _labels = []
    points = []
    diameters = []
    for edge, n_points, label_id in zip(sections["edge"], sections["n_points"], sections["label"]):
        u, v = np.fromstring(edge, dtype=int, sep=" ")
        point = []
        diameter = []
        for _ in range(int(n_points)):
            point.append(np.fromstring(sections["point"][point_id], dtype=float, sep=" "))
            diameter.append(2 * float(sections["radius"][point_id]))
            point_id += 1
        us.append(u)
        vs.append(v)
        _labels.append(labels[int(label_id)])
        points.append([list(p) for p in point])
        diameters.append(diameter)
    df["u"] = us
    df["v"] = vs
    df["label"] = _labels
    df["points"] = points
    df["diameters"] = diameters
    dfs = {}
    for label, _df in df.groupby("label"):
        dfs[label] = _df.sort_values(by="u")

    return dfs


def make_soma(dfs, morph):
    """Make a soma."""
    if "Soma" in dfs:
        soma_points = []
        soma_diameters = []
        for gid in dfs["Soma"].index:
            soma_points += [list(p) for p in dfs["Soma"].loc[gid, "points"]]
            soma_diameters += list(dfs["Soma"].loc[gid, "diameters"])
        morph.soma.points = soma_points
        morph.soma.diameters = soma_diameters
        soma_ids = dfs["Soma"]["v"].to_list()
        del dfs["Soma"]
    else:
        raise Exception("No Soma found")
    return soma_ids


def make_morph(dfs):
    """Make morphology."""
    _convert = {
        "Axon": SectionType.axon,
        "BasalDendrite": SectionType.basal_dendrite,
        "ApicalDendrite": SectionType.apical_dendrite,
    }

    morph = Morphology()
    soma_ids = make_soma(dfs, morph)
    for label, _df in dfs.items():
        section_type = _convert[label]
        _df["id"] = -1
        root_gids = _df[_df.u.isin(soma_ids)].index
        for root_gid in root_gids:
            pts = PointLevel(_df.loc[root_gid, "points"], _df.loc[root_gid, "diameters"])
            sec = morph.append_root_section(pts, section_type=section_type)
            _df.loc[root_gid, "id"] = sec.id

        for gid in _df[_df.id == -1].index:
            if _df.loc[gid, "u"] in _df.loc[_df.id > -1, "v"].to_list():
                pts = PointLevel(_df.loc[gid, "points"], _df.loc[gid, "diameters"])
                _id = int(_df.loc[_df["v"] == _df.loc[gid, "u"], "id"])
                sec = morph.sections[_id].append_section(pts)
                _df.loc[gid, "id"] = sec.id
    return morph


def load_amira(filename):
    """Load amira morphology file into morphio.mut.Morphology object."""
    with open(filename) as f:
        lines = [line.rstrip() for line in f.readlines() if "&" not in line]

    sections = get_section_data(lines)
    labels = get_labels(lines)[0]
    dfs = create_dfs(sections, labels)
    morph = make_morph(dfs)
    return morph
