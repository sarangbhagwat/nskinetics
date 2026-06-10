# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 21:24:59 2026

@author: sarangbhagwat
"""

from pathlib import Path


def antimony_txt_to_sbml_xml(
    antimony_txt_path: str,
    sbml_xml_path: str = None,
) -> Path:
    """
    Convert an Antimony .txt file to an SBML .xml file.

    Parameters
    ----------
    antimony_txt_path:
        Path to the input Antimony text file.
    sbml_xml_path:
        Optional output path. If omitted, writes next to the input file
        using the same stem and a .xml extension.

    Returns
    -------
    Path
        Path to the written SBML XML file.

    Raises
    ------
    FileNotFoundError
        If the Antimony input file does not exist.
    RuntimeError
        If Tellurium fails to convert the model.
    """
    try:
        import tellurium as te
    except ImportError as exc:
        raise ImportError(
            "Tellurium is required. Install it with: pip install tellurium"
        ) from exc

    antimony_txt_path = Path(antimony_txt_path)

    if not antimony_txt_path.exists():
        raise FileNotFoundError(f"Input file not found: {antimony_txt_path}")

    if sbml_xml_path is None:
        sbml_xml_path = antimony_txt_path.with_suffix(".xml")
    else:
        sbml_xml_path = Path(sbml_xml_path)

    antimony_text = antimony_txt_path.read_text(encoding="utf-8")

    try:
        sbml_text = te.antimonyToSBML(antimony_text)
    except Exception as exc:
        raise RuntimeError(
            f"Failed to convert Antimony file to SBML: {antimony_txt_path}"
        ) from exc

    sbml_xml_path.write_text(sbml_text, encoding="utf-8")

    return sbml_xml_path