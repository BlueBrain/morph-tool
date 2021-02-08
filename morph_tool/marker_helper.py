"""Module to read/write/plot morphology markers."""
import yaml
from pathlib import Path
import numpy as np
import logging

logger = logging.getLogger(__name__)


class MarkerSet:
    """Container of several markers for a single morphology."""

    def __init__(self, markers):
        """Create list of markers.

        Either a path to a yaml file, or a list of Marker objects.
        """
        if isinstance(markers, str):
            with open(markers, "r") as f:
                self._from_dicts(yaml.safe_load(f))
        else:
            morph_names = {marker.morph_name for marker in markers}
            if len(morph_names) > 1:
                raise Exception("Markers of different morphologies were provided.")
            self.markers = markers
            self._morph_name = morph_names.pop()
            self.morph_path = markers[0].morph_path

    def _from_dicts(self, markers):
        """Load marker from dict."""
        self.markers = [Marker(**marker) for marker in markers["markers"]]
        self._morph_name = markers["morph_name"]
        self.morph_path = markers["morph_path"]

    @property
    def morph_name(self):
        """"""
        return self._morph_name

    def plot(self, filename="markers.html", with_plotly=True):
        """Plot morphology with markers."""
        if with_plotly:
            from neurom import load_neuron
            from neurom.geom import bounding_box
            from plotly_helper.neuron_viewer import NeuronBuilder
            from plotly_helper.object_creator import scatter, vector

            neuron = load_neuron(self.morph_path)
            builder = NeuronBuilder(neuron, "3d", line_width=4, title=f"{self.morph_name}")

            for marker in self.markers:
                data = np.array(marker._data)
                if marker.type == "points":
                    if len(np.shape(marker._data)) == 1:
                        data = data[np.newaxis]
                    builder.helper.add_data(
                        {
                            f"{marker._label}": scatter(
                                data, name=f"{marker._label}", showlegend=True, **marker.plot_style
                            )
                        }
                    )
                elif marker.type == "line":
                    point_x = np.array(marker._data[0])
                    point_y = np.array(marker._data[1])

                    bbox = bounding_box(neuron)
                    bbox[0] -= np.array([10, 10, 10])
                    bbox[1] += np.array([10, 10, 10])

                    fac = 1.0
                    point1 = point_x.copy()
                    while (point1 > bbox[0]).all() and (point1 < bbox[1]).all():
                        fac += 10.0
                        point1 = point_x - fac * (point_y - point_x)

                    fac = 1.0
                    point2 = point_x.copy()
                    while (point2 > bbox[0]).all() and (point2 < bbox[1]).all():
                        fac += 10.0
                        point2 = point_x + fac * (point_y - point_x)

                    builder.helper.add_data(
                        {
                            f"{marker._label}": vector(
                                point1,
                                point2,
                                name=f"{marker._label}",
                                showlegend=True,
                                **marker.plot_style,
                            )
                        }
                    )

                else:
                    logger.info(f"marker type {marker.type} not understood")
            builder.plot(filename=str(filename))
        else:
            raise Exception("Only plotly available for now")

    def to_dicts(self):
        """Create a list of dicts from marker class."""
        out_dict = {"morph_name": self._morph_name, "morph_path": self.morph_path, "markers": []}
        for marker in self.markers:
            _d = marker.to_dict()
            del _d["morph_name"]
            del _d["morph_path"]
            out_dict["markers"].append(_d)
        return out_dict

    def save(self, filepath=".", filename=None):
        """Save marker data.

        Args:
            filepath  (str): path to folder to save marker
            filename (str): if provided,  path to marker file

        TODO: append to the list if file exists.
        """
        if filename is None:
            filename = (Path(filepath) / ("markers_" + self.morph_name)).with_suffix(".yaml")
        with open(filename, "w") as f:
            yaml.dump(self.to_dicts(), f)


class Marker:
    """Container of a single marker."""

    def __init__(self, label, tpe, data, morph_name=None, morph_path=None, plot_style=None):
        """Create a marker.

        Args:
            label (str): label of marker
            tpe (str): marker types (points, line, plane, etc...)
            data (str): markere data
            morph_name (str): name of corresponding morphology
            morph_path (str): path to corresponding morphology
            plot_style (dict): custom dict for ploting arguments
        """
        self._label = label
        self._tpe = tpe
        self._data = data
        self._morph_name = morph_name
        self.morph_path = morph_path
        self._set_plot_style(plot_style)

        self._check_valid()

    def _set_plot_style(self, plot_style):
        """Set plotting style."""
        self.plot_style = {}
        if self._tpe == "points":
            self.plot_style = {"color": "black", "width": 10}
        if plot_style is not None:
            self.plot_style.update(plot_style)

    @property
    def label(self):
        """"""
        return self._label

    @property
    def type(self):
        """"""
        return self._tpe

    @property
    def data(self):
        """"""
        return self._data

    @property
    def morph_name(self):
        """"""
        return self._morph_name

    def _check_valid(self):
        """Check if the given marker data is valid."""
        if True:
            return True
        raise Exception("Given marker data is not valid.")

    @property
    def list_data(self):
        """Returns data without numpy types.

        WIP: so that we can yaml/json easily
        """
        return np.array(self._data, dtype=float).tolist()

    def to_dict(self):
        """Create a dict from marker class."""
        return {
            "morph_name": self._morph_name,
            "morph_path": self.morph_path,
            "label": self._label,
            "tpe": self._tpe,
            "data": self.list_data,
            "plot_style": self.plot_style,
        }
