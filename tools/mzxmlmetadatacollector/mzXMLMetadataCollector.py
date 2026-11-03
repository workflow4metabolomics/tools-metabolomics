# -*- coding: utf-8 -*-
"""
Created on Dec 04 2023
@author: quruin

Last update on Nov 14 2024

@author: quruin

Every metadata entry is extracted with the ElementTree XML parser. Headers are
streamed and parsing stops at the end of the first spectrum, so the base64
encoded spectral payloads of the remaining scans are never held in memory.
Element names are stripped of their namespace prefix, which lets the lookups
below ignore the difference between the mzML and mzXML schemas. Entries that
describe a single spectrum are read from the first spectrum itself, not from
the parsed tree, so that they do not depend on how far the parser read ahead.

The lookup helpers use the non raising ``find`` and ``get`` methods of the
parser and return the "Not found" placeholder when an entry is missing from a
file: an absent metadata entry is a normal case, not an error. Genuine errors
(missing file, malformed XML, unsupported extension, inconsistent arguments)
are not caught and make the tool fail.
"""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

NOT_FOUND = "Not found"
NOT_AVAILABLE = "Not available"
MZXML_NOT_AVAILABLE = "Not available in mzXML"

# Accessions describing the MS level, ordered by decreasing priority.
MS_LEVEL_ACCESSIONS = ("MS:1000580", "MS:1000579")
# Accession carrying the MS level as a plain number, used as a fallback for
# files which declare no MS1/MSn label at all.
MS_LEVEL_NUMERIC_ACCESSION = "MS:1000511"
# Accessions of the two acquisition modes of a spectrum in mzML.
SPECTRUM_TYPE_ACCESSIONS = ("MS:1000127", "MS:1000128")
# Accessions holding the instrument name, with the attribute to read from.
INSTRUMENT_MODEL_ACCESSIONS = (
    ("MS:1000703", "name"),
    ("MS:1001547", "name"),
    ("MS:1003094", "name"),
    ("MS:1003123", "name"),
    ("MS:1000495", "value"),
    ("MS:1000483", "name"),
)
# mzXML files carry the same instrument accessions, always in the name.
MZXML_INSTRUMENT_MODEL_ACCESSIONS = ("MS:1001547", "MS:1003123", "MS:1000495")

# Root elements each format is allowed to declare, checked against the file
# name so that a mislabelled file fails instead of producing an empty row.
MZML_ROOT_TAGS = ("mzML", "indexedmzML")
MZXML_ROOT_TAGS = ("mzXML",)
# Elements holding encoded data, and the attributes by which they declare the
# number of points: mzML uses defaultArrayLength, mzXML peaksCount. mzXML names
# its spectrum element "scan".
SPECTRUM_ELEMENTS = ("spectrum", "scan")
CHROMATOGRAM_ELEMENTS = ("chromatogram",)
POINTS_ATTRIBUTES = ("defaultArrayLength", "peaksCount")

HEADER = ("Name\tformat\tversion\tsize(ko)\tsize(Mo)\tsize(Go)\tMSlevel\t"
          "Spectrum type\tSource file\tAcquisition date\tSoftware(s) used\t"
          "Processing method(s)\tMachine\tNumber of scans\tNumber of points\t"
          "Number of chromatogram points\tEncoding\n")


def get_info(path, name):
    """Collect the metadata of one file as the fields of a table row."""
    extension = name.split(".")[-1]

    if extension.casefold() == "mzml":
        root, spectrum = parse_header(path, "spectrum")
        check_root(root, MZML_ROOT_TAGS, name)
        params = spectrum_params(root, spectrum)
        points, chromatogram_points = count_points(path)
        return (
            extension,
            get_version(root),
            get_kb_size(path),
            get_mb_size(path),
            get_gb_size(path),
            get_ms_level(params),
            get_spectrum_type(params),
            get_source_file(root),
            get_acquisition_date(root),
            get_softwares_list(root),
            get_processing_list(root),
            get_instrument_model(root),
            count_scans(root),
            points,
            chromatogram_points,
            get_encoding(params),
        )

    if extension.casefold() == "mzxml":
        root, spectrum = parse_header(path, "scan")
        check_root(root, MZXML_ROOT_TAGS, name)
        points, chromatogram_points = count_points(path)
        return (
            extension,
            get_mzxml_version(root),
            get_kb_size(path),
            get_mb_size(path),
            get_gb_size(path),
            MZXML_NOT_AVAILABLE,  # the MS level is not collected for mzXML
            get_mzxml_spectrum_type(spectrum),
            get_source_mzxml(root),
            MZXML_NOT_AVAILABLE,  # no acquisition date is stored in mzXML
            get_mzxml_software_list(root),
            MZXML_NOT_AVAILABLE,  # no processing method is stored in mzXML
            get_mzxml_instrument_model(root),
            count_scans_mzxml(root),
            points,
            chromatogram_points,
            get_encoding_mzxml(spectrum),
        )

    raise SystemExit(
        f"{name}: unsupported extension '{extension}', expected mzML or mzXML")


def parse_header(path, spectrum_tag):
    """Parse the metadata header of ``path``.

    Returns the ``(root, spectrum)`` pair, where ``spectrum`` is the first
    completed spectrum element -- ``spectrum_tag`` is "spectrum" for mzML and
    "scan" for mzXML -- or None when the file declares none.

    Streaming stops at the end of that first spectrum, which keeps the encoded
    payloads of the remaining scans out of memory. The parser reads in blocks,
    so a few elements after the first spectrum may still be in the tree: the
    per spectrum lookups below therefore search the returned element itself
    rather than the tree, so that their result cannot depend on how much of
    the document the parser happened to read ahead.
    """
    root = None
    spectrum = None
    with open(path, "rb") as handle:
        for event, element in ET.iterparse(handle,
                                           events=("start", "end")):
            if root is None:
                root = element
            if event == "end" and local_name(element.tag) == spectrum_tag:
                spectrum = element
                break
    return strip_namespaces(root), spectrum


def check_root(root, tags, name):
    """Fail when the document is not of the format the file name announces."""
    if root.tag not in tags:
        raise SystemExit(
            f"{name}: declares the '{name.split('.')[-1]}' extension but its"
            f" root element is '{root.tag}', not one of"
            f" {' or '.join(tags)}")


def spectrum_params(root, spectrum):
    """``cvParam`` elements describing the file content and first spectrum.

    This is the region the MS level, the spectrum type and the encoding are
    read from. Bounding it to the file level declaration and to the first
    spectrum keeps the reported values independent of the read ahead that
    ``parse_header`` describes.
    """
    params = []
    for element in (root.find(".//fileContent"), spectrum):
        if element is not None:
            params.extend(element.iter("cvParam"))
    return params


def local_name(tag):
    """Element name without its ``{namespace}`` prefix."""
    return tag.rpartition("}")[2]


def strip_namespaces(root):
    """Drop the namespace prefix of every element name of the tree.

    The mzML and mzXML schemas declare a default namespace, which ElementTree
    prefixes to every element name. Dropping the prefixes lets the lookups
    below use plain element names for both formats.
    """
    for element in root.iter():
        element.tag = local_name(element.tag)
    return root


def count_points(path):
    """Sum the declared points of the spectra and of the chromatograms.

    Returns a ``(spectra, chromatograms)`` pair. The whole file has to be
    walked, as the headers do not report these totals. Every element is
    discarded as soon as it is complete -- ``clear`` empties it and it is
    detached from its parent, which ``clear`` alone does not do -- so the
    memory used stays independent of the number of scans.
    """
    spectra = 0
    chromatograms = 0
    open_elements = []
    with open(path, "rb") as handle:
        for event, element in ET.iterparse(handle,
                                           events=("start", "end")):
            if event == "start":
                open_elements.append(element)
                continue
            name = local_name(element.tag)
            if name in SPECTRUM_ELEMENTS:
                spectra += declared_points(element)
            elif name in CHROMATOGRAM_ELEMENTS:
                chromatograms += declared_points(element)
            element.clear()
            open_elements.pop()
            if open_elements:
                open_elements[-1].remove(element)
    return spectra, chromatograms


def declared_points(element):
    """Number of encoded points declared by an element, 0 if it has none."""
    for attribute in POINTS_ATTRIBUTES:
        value = element.attrib.get(attribute)
        if value is not None:
            return int(value)
    return 0


def cv_param_name(element):
    """Name of the first ``cvParam`` of ``element``, None if it has none.

    Nested parameters are included, as some writers place the ``cvParam`` of a
    processing method or of a software inside a ``paramList`` child.
    """
    cv_param = element.find(".//cvParam")
    if cv_param is None:
        return None
    return cv_param.get("name")


def get_kb_size(path):
    return np.round(os.stat(path).st_size / 1024, 2)


def get_mb_size(path):
    return np.round(os.stat(path).st_size / 1024**2, 2)


def get_gb_size(path):
    return np.round(os.stat(path).st_size / 1024**3, 2)


def get_version(root):
    """Version declared by the ``mzML`` element of the file."""
    mzml = root if root.tag == "mzML" else root.find(".//mzML")
    if mzml is None:
        return NOT_FOUND
    return mzml.get("version", NOT_FOUND)


def get_ms_level(params):
    """Name of the MS level description of the file and first spectrum.

    Files which declare the level only as a number -- the ``ms level`` term
    carries it in its value -- are reported as "MS<number> spectrum".
    """
    for accession in MS_LEVEL_ACCESSIONS:
        for cv_param in params:
            if cv_param.get("accession") == accession:
                name = cv_param.get("name")
                if name:
                    return name
    for cv_param in params:
        if cv_param.get("accession") == MS_LEVEL_NUMERIC_ACCESSION:
            level = cv_param.get("value", "")
            if level:
                return f"MS{level} spectrum"
    return NOT_FOUND


def get_spectrum_type(params):
    """Name of the acquisition mode, centroid or profile, of the spectrum.

    The first declaration found wins, so a file level statement takes
    precedence over the one of the first spectrum.
    """
    for cv_param in params:
        if cv_param.get("accession") in SPECTRUM_TYPE_ACCESSIONS:
            return cv_param.get("name", NOT_FOUND)
    return NOT_FOUND


def get_source_file(root):
    """Location, or name, of the file the spectra were acquired in."""
    source_file = root.find(".//sourceFile")
    if source_file is None:
        return NOT_FOUND
    return source_file.get("location") or source_file.get("name") or NOT_FOUND


def get_acquisition_date(root):
    """Start time stamp of the run, date and time separated by a space."""
    run = root.find(".//run")
    if run is None:
        return NOT_FOUND
    timestamp = run.get("startTimeStamp")
    if not timestamp:
        return NOT_FOUND
    return timestamp.replace("T", " ", 1)


def get_softwares_list(root):
    """Every declared software as "name version", joined by " + "."""
    softwares = []
    for software in root.iter("software"):
        name = (software.get("id") or software.get("name")
                or cv_param_name(software))
        if not name:
            # Without a name, a bare version number carries no information.
            softwares.append(NOT_FOUND)
            continue
        version = software.get("version", "")
        softwares.append(f"{name} {version}".strip())
    if not softwares:
        return NOT_FOUND
    return " + ".join(softwares)


def get_processing_list(root):
    """Every processing method as "software method", joined by " + "."""
    methods = []
    for method in root.iter("processingMethod"):
        software = method.get("softwareRef", "")
        # Some writers place the name on the element rather than in a child.
        name = cv_param_name(method) or method.get("name") or NOT_FOUND
        methods.append(f"{software} {name}".strip())
    if not methods:
        return NOT_FOUND
    return " + ".join(methods)


def get_instrument_model(root):
    """Name of the instrument the spectra were acquired on."""
    for accession, attribute in INSTRUMENT_MODEL_ACCESSIONS:
        for cv_param in root.iter("cvParam"):
            if cv_param.get("accession") == accession:
                return cv_param.get(attribute) or NOT_AVAILABLE
    return NOT_AVAILABLE


def count_scans(root):
    """Number of spectra declared by the ``spectrumList`` element."""
    spectrum_list = root.find(".//spectrumList")
    if spectrum_list is None:
        return NOT_FOUND
    return spectrum_list.get("count", NOT_FOUND)


def get_encoding(params):
    """Name of the float encoding of the spectral data, e.g. 64-bit float."""
    for cv_param in params:
        name = cv_param.get("name", "")
        if name.endswith("-bit float"):
            return name
    return NOT_FOUND


def get_mzxml_version(root):
    """mzXML schema version, read from the attributes of the root element."""
    for value in root.attrib.values():
        if "mzXML_" in value:
            return value.split("mzXML_")[1].split(" ")[0]
    return root.get("version", NOT_FOUND)


def get_mzxml_spectrum_type(scan):
    """Acquisition mode of the first scan, "Not found" when it declares none.

    mzXML describes the centroiding mode as an attribute of every single scan,
    so only the first one is read -- the element ``parse_header`` returns.
    """
    if scan is None:
        return NOT_FOUND
    centroided = scan.get("centroided")
    if centroided == "0":
        return "profile"
    if centroided == "1":
        return "centroid"
    return NOT_FOUND


def get_source_mzxml(root):
    """Name of the native file the spectra were acquired in."""
    parent_file = root.find(".//parentFile")
    if parent_file is None:
        return NOT_FOUND
    return parent_file.get("fileName", NOT_FOUND)


def get_mzxml_software_list(root):
    """Every declared software as "type: name version", joined by " + "."""
    softwares = []
    for software in root.iter("software"):
        entry = ""
        type_ = software.get("type")
        if type_:
            entry = f"{type_}: "
        name = f"{software.get('name', '')} {software.get('version', '')}"
        softwares.append(entry + name.strip())
    if not softwares:
        return NOT_FOUND
    return " + ".join(softwares)


def get_mzxml_instrument_model(root):
    """msModel value, or the name of an instrument ``cvParam``."""
    ms_model = root.find(".//msModel")
    if ms_model is not None and ms_model.get("value"):
        return ms_model.get("value")
    for accession in MZXML_INSTRUMENT_MODEL_ACCESSIONS:
        for cv_param in root.iter("cvParam"):
            if cv_param.get("accession") == accession:
                return cv_param.get("name", NOT_FOUND)
    return NOT_FOUND


def count_scans_mzxml(root):
    """Number of scans declared by the ``msRun`` element."""
    ms_run = root.find(".//msRun")
    if ms_run is None:
        return NOT_FOUND
    return ms_run.get("scanCount", NOT_FOUND)


def get_encoding_mzxml(scan):
    """Precision of the encoded spectral data, e.g. 64-bit.

    The schema declares the precision on the ``peaks`` element, a child of the
    first scan, which is therefore searched instead of the whole document.
    """
    if scan is None:
        return NOT_AVAILABLE
    for element in scan.iter():
        precision = element.get("precision")
        if precision is not None:
            return f"{precision}-bit"
    return NOT_AVAILABLE


def update_entries(paths, names):
    """Build the content of the output: a header and one row per file."""
    rows = [HEADER]
    for path, name in zip(paths, names):
        fields = [name, *get_info(path, name)]
        rows.append("\t".join(str(field) for field in fields) + "\t\n")
    return "".join(rows)


def main():
    print("Parameters passed to the tool:\n")
    print(list(sys.argv))

    if len(sys.argv) < 5:
        raise SystemExit(
            "usage: mzXMLMetadataCollector.py OUTPUT_FILE SINGLE|MULTIPLE"
            " PATHS NAMES")

    output_file = sys.argv[1]
    selection = sys.argv[2]
    if selection == "SINGLE":
        paths = [sys.argv[3]]
    elif selection == "MULTIPLE":
        paths = sys.argv[3].split(",")
    else:
        raise SystemExit(f"unknown number of files selection '{selection}'")

    names = sys.argv[4:]
    if len(paths) != len(names):
        raise SystemExit(
            f"{len(paths)} file path(s) but {len(names)} file name(s) passed")

    with open(output_file, "w", encoding="utf-8") as result:
        result.write(update_entries(paths, names))


if __name__ == "__main__":
    main()
