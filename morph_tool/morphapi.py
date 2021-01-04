"""Module for morphology API.

some parts are copied/adapted from morphology_repair_workflow/morphdb.py
"""
from collections import defaultdict, namedtuple
import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd

NEURONDB_XML = "neuronDB.xml"
MTYPE_MSUBTYPE_SEPARATOR = ":"

# Note: mtype is the 'combined mtype'
XMLMorphInfo = namedtuple(
    "XMLMorphInfo", "name mtype layer use_axon use_dendrites raw_xml"
)
_FIX_NAMES = {"mtype": "fullmtype", "name": "morph_name"}


def morphdb_rename(db_name):
    """Remove any int at the front and set lowercase."""
    if db_name.split("_")[0].isdigit():
        db_name = "_".join(db_name.split("_")[1:])
    return db_name.lower()


def _morph_path_rename(df, db_name):
    """Clean the folder names.

    Remove any int at the front and set lowercase.
    """
    return df.rename(columns={"morph_path": "morph_path_" + morphdb_rename(db_name)})


def _get_combined_mtype(mtype, msubtype):
    """combine the mtype and msubtype if it exists

    Most morphologies don't have a subtype, but for those who do
    we need to handle the combination of type and subtype as unique, thus,
    they are joined them using the ':' separator which shouldn't exist
    in any morphology (sub)type names
    """
    if msubtype:
        mtype += MTYPE_MSUBTYPE_SEPARATOR + msubtype
    return mtype


def _xmlmorphinfo_from_xml(xml_morph):
    """extracts properties from a neurondb.xml <morphology> stanza"""
    name = xml_morph.findtext("name")
    mtype = xml_morph.findtext("mtype")
    msubtype = xml_morph.findtext("msubtype")

    mtype = _get_combined_mtype(mtype, msubtype)

    layer = xml_morph.findtext("layer")
    # According to Eilif, an empty use_axon (corresponding to a null in the database)
    # means that the axon is supposed to be used
    use_axon = xml_morph.findtext(".//use_axon") in (
        "",
        "true",
        "True",
    )

    # use the same behaviour for the dendrites; if it's blank, or doesn't exist, the
    # dendrite will be used
    use_dendrites = xml_morph.findtext(".//use_dendrites") in (
        None,
        "",
        "true",
        "True",
    )

    return XMLMorphInfo(
        name=name,
        mtype=mtype,
        layer=layer,
        use_axon=use_axon,
        use_dendrites=use_dendrites,
        raw_xml=xml_morph,
    )


def extract_morphinfo_from_xml(root, filter_=None):
    """returns a generator that contains all the morphologies from `root`"""
    for morph in root.findall(".//morphology"):
        morph = _xmlmorphinfo_from_xml(morph)
        if filter_ and not filter_(morph):
            continue
        # the XML created by the MorphologyTracker is ill-formed
        # and contains simply '<morphology>\n    </morphology>\n    '
        # skip those
        if not morph.name:
            continue
        yield morph


class MorphologyDB(object):
    """hold intermediate data that phases may need

    Currently this allows for maintenance of a neuronDB.xml, but if the intermediate format is
    changed, this can be updated
    """

    def __init__(self):
        self._root = ET.Element("neurondb")
        self._listing = ET.SubElement(self._root, "listing")
        self.lineage = {}

        # keep track of known morphologies, so that we don't add duplicates
        self._known_morphologies = set()

    @staticmethod
    def _morph_key(morph):
        """make a key for the _known_morphologies dict"""
        return morph.name, morph.mtype, morph.layer

    @classmethod
    def load_neurondb_xml(cls, path):
        """load an existing neurondb.xml"""
        morphdb = cls()

        with open(path) as fd:
            # pylint: disable=W0212
            morphdb._root = ET.fromstring(fd.read())
            morphdb._listing = morphdb._root.find("listing")
            morphdb._known_morphologies = set(
                MorphologyDB._morph_key(morph) for morph in morphdb
            )
        return morphdb

    def __iter__(self):
        return extract_morphinfo_from_xml(self._root)

    def remove_morphs(self, removed_morphs):
        """remove an iterable of morph names from the xml"""
        for name in removed_morphs:
            for morph_xml in self._listing.findall(".//morphology"):
                morph = _xmlmorphinfo_from_xml(morph_xml)
                if morph.name == name:
                    self._listing.remove(morph_xml)
                    self._known_morphologies.discard(MorphologyDB._morph_key(morph))
            self.lineage.pop(name, None)

    def add_morph_from_xml(self, morph):
        """add a morphology from a xml element"""
        self._listing.append(morph)

    def add_morph(self, name, mtype, layer, parent=None, dendrite=None, axon=None):
        """add a morphology based on name, mtype and layer

        Add lineage information: if parent, then dendrite and axon are both set,
        or they are set independently

        If the name, mtype, layer combination already exists, then nothing is added

        Return:
            bool: True if morphology added, False otherwise
        """
        morph_info = XMLMorphInfo(
            name=name,
            mtype=mtype,
            layer=layer,
            use_axon=None,
            use_dendrites=None,
            raw_xml=None,
        )
        morph_key = MorphologyDB._morph_key(morph_info)
        if morph_key in self._known_morphologies:
            return False
        self._known_morphologies.add(morph_key)

        add_morphology_to_xml(self._listing, morph_info)

        parents = {}
        if parent:
            parents["dendrite"] = parents["axon"] = parent
        elif dendrite and axon:
            parents["dendrite"] = dendrite
            parents["axon"] = axon
        elif dendrite:
            parents["dendrite"] = dendrite
        elif axon:
            parents["axon"] = axon

        if parents:
            self.lineage[name] = parents

        return True

    def write_xml(self, output_path, name=NEURONDB_XML):
        """output the neuronDB.xml to output_path/name"""
        full_path = joinp(output_path, name)
        L.debug("Writing XML: %s", full_path)

        overview = self._root.findall(".//overview")
        assert len(overview) <= 1, "Can only have zero or one overview section"

        if not overview:
            overview = ET.SubElement(self._root, "overview")
            count = ET.SubElement(overview, "count")
        else:
            overview = overview[0]
            count = overview.find(".//count")
        count.text = str(len(list(self)))

        with open(full_path, "w") as fd:
            fd.write(ET.tostring(self._root))

        return full_path

    def _write_neurondb(self, output_path, sep=" ", write_header=False, filt=None):
        if filt is None:
            filt = {}
        with open(output_path, "w") as fd:
            if write_header:
                fd.write(sep.join(["name", "layer", "mtype"]) + "\n")
            for morph in self:
                if ("use_axon" in filt) and (morph.use_axon != filt["use_axon"]):
                    continue
                if ("use_dendrites" in filt) and (
                    morph.use_dendrites != filt["use_dendrites"]
                ):
                    continue
                fd.write(sep.join([morph.name, morph.layer, morph.mtype]) + "\n")

    def write_neurondb_dat(self, output_path, name="neurondb.dat", filt=None):
        """write out the neurondb.dat file
        1 line per morphology in the xml:
            name        layer( starting at 1) mtype
            C300797C-P2 4                     L4_PC
        """
        result = joinp(output_path, name)
        self._write_neurondb(result, sep=" ", write_header=False, filt=filt)
        return result

    def write_neurondb_csv(self, output_path, name="neurondb.csv", filt=None):
        """write out the neurondb as a CSV file """
        result = joinp(output_path, name)
        self._write_neurondb(result, sep=",", write_header=True, filt=filt)
        return result

    def write_lineage(self, output_path, name="lineage.json"):
        """write the lineage to the output_path/name"""
        full_path = joinp(output_path, name)
        with open(full_path, "w") as fd:
            json.dump(self.lineage, fd, indent=2)
        return full_path


class MorphAPI:
    """Helper class to access and create database of morphologies."""

    def __init__(self, root_path):
        """Initialise a morphAPI from a root folder."""
        self._root_path = Path()
        self._morphdbs = {}
        self._morphsdfs = {}
        self._load(root_path)

    @property
    def dbs(self):
        """Returns the dict of morphology databases."""
        return self._morphsdbs

    def _load_from_single_folder(self, path):
        """Load morphologies from a single folder with neurondb.xml."""
        self._morphdbs[morphdb_rename(path.stem)] = MorphologyDB.load_neurondb_xml(
            path / NEURONDB_XML
        )

    def _load_from_folders(self, root_path):
        """Load several folders as a release with various neurondb.xml."""
        for folder in self._root_path.iterdir():
            if folder.is_dir():
                self._morphdbs[
                    morphdb_rename(folder.stem)
                ] = MorphologyDB.load_neurondb_xml(folder / NEURONDB_XML)

    def _load(self, root_path):
        """Load a morphology release."""
        self._root_path = Path(root_path).absolute()
        if (self._root_path / NEURONDB_XML).exists():
            self._load_from_single_folder(self._root_path)
        elif any(
            (folder / NEURONDB_XML).exists() for folder in self._root_path.iterdir()
        ):
            self._load_from_folders(self._root_path)
        else:
            raise Exception(f"We cannot load morphologies from path {root_path}")

    def as_dataframes(self):
        "Return a dataframe per morphology database." ""
        if self._morphsdfs:
            return self._morphsdfs

        _skip_columns = ["raw_xml"]
        columns = [col for col in XMLMorphInfo._fields if col not in _skip_columns]
        dfs = {}
        for db_name, db in self._morphdbs.items():
            _data = defaultdict(list)
            for morph in db:
                for column in columns:
                    _data[column].append(getattr(morph, column))

            _df = pd.DataFrame()
            for column in columns:
                _df[_FIX_NAMES.get(column, column)] = _data[column]
            _df["morph_path"] = str(self._root_path) + "/" + _df["morph_name"]

            dfs[db_name] = _df
        self._morphsdfs = dfs
        return dfs

    def as_collapsed_dataframes(self, _first=0):
        """Collapses the dataframes into a single one by duplicating morph_path columns."""
        dfs = self.as_dataframes()

        _other_dfs, _df_len = [], []
        first_db = list(dfs.keys())[_first]
        _df = _morph_path_rename(dfs[first_db], first_db)
        for i, (db_name, df) in enumerate(dfs.items()):
            if len(df.index) == len(dfs[first_db].index):
                _df = pd.merge(_df, _morph_path_rename(df, db_name))
            elif _first == 0 and len(df.index) not in _df_len:
                _other_dfs.append(self.as_collapsed_dataframes(_first=i))
                _df_len.append(len(df.index))
        if len(_other_dfs) > 0:
            return [_df] + _other_dfs
        return _df

    def as_single_dataframe(self, db_name="repairunravel-asc"):
        """Load a single morph db into a dataframe."""
        if db_name not in self._morphdbs:
            raise Exception(f"{db_name} is not a valid morphdb name")

        return self.as_dataframes()[db_name]
