"""Module to read/write/plot morphology markers."""
import yaml
from pathlib import Path
import numpy as np


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
            self._markers = markers
            self._morph_name = morph_names.pop()
            self.morph_path = markers[0].morph_path

    def _from_dicts(self, markers):
        """Load marker from dict."""
        self._markers = markers["markers"]
        self._morph_name = markers["morph_name"]
        self.morph_path = markers["morph_path"]

    @property
    def morph_name(self):
        """"""
        return self._morph_name

    def plot(self, with_plotly=True, filename="markers.html", plot_kwargs=None):
        """Plot morphology with markers."""
        if with_plotly:
            from neurom import load_neuron
            from plotly_helper.neuron_viewer import NeuronBuilder
            from plotly_helper.object_creator import scatter

            builder = NeuronBuilder(
                load_neuron(self.morph_path), "3d", line_width=4, title=f"{self.morph_name}"
            )
            for marker in self._markers:
                data = np.array(marker._data)
                if len(np.shape(marker._data)) == 1:
                    data = marker._data[np.newaxis]
                builder.helper.add_data(
                    {
                        f"{marker._label}": scatter(
                            data, name=f"{marker._label}", showlegend=True, **marker.plot_style
                        )
                    }
                )
            builder.plot(filename=filename)

    def to_dicts(self):
        """Create a list of dicts from marker class."""
        out_dict = {"morph_name": self._morph_name, "morph_path": self.morph_path, "markers": []}
        for marker in self._markers:
            _d = marker.to_dict()
            del _d["morph_name"]
            del _d["morph_path"]
            out_dict["markers"].append(_d)
        return out_dict

    def save(self, filepath="."):
        """Save marker data.

        TODO: append to the list if file exists.
        """
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
        self._type = tpe
        self._data = data
        self._morph_name = morph_name
        self.morph_path = morph_path
        self._set_plot_style(plot_style)

        self._check_valid()

    def _set_plot_style(self, plot_style):
        if plot_style is None:
            if self._type == "points":
                self.plot_style = {"color": "black", "width": 10}

    @property
    def label(self):
        """"""
        return self._label

    @property
    def type(self):
        """"""
        return self._type

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
        return self._data.tolist()

    def to_dict(self):
        """Create a dict from marker class."""
        return {
            "morph_name": self._morph_name,
            "morph_path": self.morph_path,
            "label": self._label,
            "type": self._type,
            "data": self.list_data,
        }
