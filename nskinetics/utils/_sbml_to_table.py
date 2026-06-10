# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 21:24:59 2026

@author: sarangbhagwat
"""

from __future__ import annotations

import html
import re
from pathlib import Path
from typing import Any, Dict, Optional, Union

import libsbml
import pandas as pd


PathLike = Union[str, Path]


# =============================================================================
# Public entry-point function
# =============================================================================

def export_sbml_to_tables(
    sbml_file_path: PathLike,
    output_excel_path: Optional[PathLike] = None,
    *,
    raise_on_sbml_warnings: bool = False,
    print_summary: bool = True,
) -> Dict[str, pd.DataFrame]:
    """
    Export an SBML kinetic model to tabular Excel sheets.

    Parameters
    ----------
    sbml_file_path:
        Path to the input SBML file, usually .xml or .sbml.

    output_excel_path:
        Path to the output Excel workbook. If None, the output will be created
        next to the SBML file as '<input_stem>_sbml_tables.xlsx'.

    raise_on_sbml_warnings:
        If True, SBML warnings are treated as errors. By default, only serious
        SBML errors stop execution.

    print_summary:
        If True, prints the output path and number of rows in each exported table.

    Returns
    -------
    Dict[str, pandas.DataFrame]
        Dictionary of exported tables. Keys are sheet/table names.
    """
    sbml_path = Path(sbml_file_path)

    if not sbml_path.exists():
        raise FileNotFoundError(f"SBML file not found: {sbml_path}")

    if output_excel_path is None:
        excel_path = sbml_path.with_name(f"{sbml_path.stem}_sbml_tables.xlsx")
    else:
        excel_path = Path(output_excel_path)

    document = read_sbml_document(
        sbml_path,
        raise_on_sbml_warnings=raise_on_sbml_warnings,
    )

    tables = extract_all_tables(document)
    export_tables_to_excel(tables, excel_path)

    if print_summary:
        print(f"Wrote Excel workbook: {excel_path}")
        for table_name, df in tables.items():
            print(f"  {table_name}: {len(df)} rows")

    return tables


# =============================================================================
# Helper functions
# =============================================================================

def has_method(obj: Any, method: str) -> bool:
    return hasattr(obj, method) and callable(getattr(obj, method))


def get_if_set(obj: Any, getter: str, is_setter: Optional[str] = None) -> Any:
    """
    Safely call a libSBML getter if present.

    If is_setter is provided, the getter is only called when the corresponding
    isSetX method returns True.
    """
    if not has_method(obj, getter):
        return ""

    if is_setter and has_method(obj, is_setter):
        try:
            if not getattr(obj, is_setter)():
                return ""
        except Exception:
            return ""

    try:
        value = getattr(obj, getter)()
        return "" if value is None else value
    except Exception:
        return ""


def math_to_string(math_ast: Any) -> str:
    """
    Convert an SBML Math ASTNode to a readable formula string.
    """
    if math_ast is None:
        return ""

    try:
        return libsbml.formulaToL3String(math_ast)
    except Exception:
        try:
            return libsbml.formulaToString(math_ast)
        except Exception:
            return ""


def clean_notes(notes_string: str) -> str:
    """
    Convert SBML notes XHTML into rough readable text.
    """
    if not notes_string:
        return ""

    text = re.sub(r"<[^>]+>", " ", notes_string)
    text = html.unescape(text)
    return re.sub(r"\s+", " ", text).strip()


def common_fields(obj: Any) -> Dict[str, Any]:
    """
    Fields common to most libSBML SBase-derived objects.
    """
    return {
        "id": get_if_set(obj, "getId", "isSetId"),
        "name": get_if_set(obj, "getName", "isSetName"),
        "meta_id": get_if_set(obj, "getMetaId", "isSetMetaId"),
        "sbo_term": get_if_set(obj, "getSBOTermID", "isSetSBOTerm"),
        "notes": clean_notes(get_if_set(obj, "getNotesString") or ""),
        "annotation_xml": get_if_set(obj, "getAnnotationString") or "",
    }


def read_sbml_document(
    sbml_path: Path,
    *,
    raise_on_sbml_warnings: bool = False,
) -> libsbml.SBMLDocument:
    """
    Read and validate an SBML document.
    """
    reader = libsbml.SBMLReader()
    document = reader.readSBML(str(sbml_path))

    serious_errors = []
    warnings = []

    for i in range(document.getNumErrors()):
        error = document.getError(i)
        message = error.getMessage()
        severity = error.getSeverity()

        if severity >= libsbml.LIBSBML_SEV_ERROR:
            serious_errors.append(message)
        else:
            warnings.append(message)

    if serious_errors:
        joined = "\n".join(f"- {msg}" for msg in serious_errors[:10])
        raise ValueError(
            f"SBML parsing failed with {len(serious_errors)} serious error(s):\n"
            f"{joined}"
        )

    if warnings and raise_on_sbml_warnings:
        joined = "\n".join(f"- {msg}" for msg in warnings[:10])
        raise ValueError(
            f"SBML parsing produced {len(warnings)} warning(s):\n"
            f"{joined}"
        )

    model = document.getModel()
    if model is None:
        raise ValueError("No SBML model found in file.")

    return document


def unit_definition_to_string(unit_def: libsbml.UnitDefinition) -> str:
    """
    Make a compact textual representation of an SBML UnitDefinition.

    Example:
        mole litre^-1 second^-1
    """
    parts = []

    for unit in unit_def.getListOfUnits():
        kind = libsbml.UnitKind_toString(unit.getKind())
        exponent = unit.getExponent()
        scale = unit.getScale()
        multiplier = unit.getMultiplier()

        term = kind

        if multiplier != 1:
            term = f"{multiplier:g}*{term}"

        if scale != 0:
            term = f"10^{scale}*{term}"

        if exponent != 1:
            term = f"{term}^{exponent:g}"

        parts.append(term)

    return " ".join(parts)


# =============================================================================
# Table extraction functions
# =============================================================================

def extract_all_tables(document: libsbml.SBMLDocument) -> Dict[str, pd.DataFrame]:
    """
    Extract all relevant SBML entities into pandas DataFrames.
    """
    model = document.getModel()

    tables = {
        "model": extract_model_metadata(document),
        "unit_definitions": extract_unit_definitions(model),
        "compartments": extract_compartments(model),
        "species": extract_species(model),
        "parameters": extract_parameters(model),
        "initial_assignments": extract_initial_assignments(model),
        "rules": extract_rules(model),
        "function_definitions": extract_function_definitions(model),
        "constraints": extract_constraints(model),
    }

    tables.update(extract_reactions(model))
    tables.update(extract_events(model))

    return tables


def extract_model_metadata(document: libsbml.SBMLDocument) -> pd.DataFrame:
    model = document.getModel()

    row = {
        "sbml_level": document.getLevel(),
        "sbml_version": document.getVersion(),
        "model_id": get_if_set(model, "getId", "isSetId"),
        "model_name": get_if_set(model, "getName", "isSetName"),
        "model_meta_id": get_if_set(model, "getMetaId", "isSetMetaId"),
        "substance_units": get_if_set(model, "getSubstanceUnits", "isSetSubstanceUnits"),
        "time_units": get_if_set(model, "getTimeUnits", "isSetTimeUnits"),
        "volume_units": get_if_set(model, "getVolumeUnits", "isSetVolumeUnits"),
        "area_units": get_if_set(model, "getAreaUnits", "isSetAreaUnits"),
        "length_units": get_if_set(model, "getLengthUnits", "isSetLengthUnits"),
        "extent_units": get_if_set(model, "getExtentUnits", "isSetExtentUnits"),
        "conversion_factor": get_if_set(
            model,
            "getConversionFactor",
            "isSetConversionFactor",
        ),
        "sbo_term": get_if_set(model, "getSBOTermID", "isSetSBOTerm"),
        "notes": clean_notes(get_if_set(model, "getNotesString") or ""),
        "annotation_xml": get_if_set(model, "getAnnotationString") or "",
    }

    return pd.DataFrame([row])


def extract_unit_definitions(model: libsbml.Model) -> pd.DataFrame:
    rows = []

    for unit_def in model.getListOfUnitDefinitions():
        base = common_fields(unit_def)
        base["unit_expression"] = unit_definition_to_string(unit_def)

        if unit_def.getNumUnits() == 0:
            rows.append(base)
            continue

        for i, unit in enumerate(unit_def.getListOfUnits()):
            row = dict(base)
            row.update(
                {
                    "unit_index": i,
                    "kind": libsbml.UnitKind_toString(unit.getKind()),
                    "exponent": unit.getExponent(),
                    "scale": unit.getScale(),
                    "multiplier": unit.getMultiplier(),
                }
            )
            rows.append(row)

    return pd.DataFrame(rows)


def extract_compartments(model: libsbml.Model) -> pd.DataFrame:
    rows = []

    for compartment in model.getListOfCompartments():
        row = common_fields(compartment)
        row.update(
            {
                "spatial_dimensions": get_if_set(
                    compartment,
                    "getSpatialDimensions",
                    "isSetSpatialDimensions",
                ),
                "size": get_if_set(compartment, "getSize", "isSetSize"),
                "units": get_if_set(compartment, "getUnits", "isSetUnits"),
                "constant": get_if_set(compartment, "getConstant", "isSetConstant"),
                "outside": get_if_set(compartment, "getOutside", "isSetOutside"),
            }
        )
        rows.append(row)

    return pd.DataFrame(rows)


def extract_species(model: libsbml.Model) -> pd.DataFrame:
    rows = []

    for species in model.getListOfSpecies():
        row = common_fields(species)
        row.update(
            {
                "compartment": get_if_set(species, "getCompartment", "isSetCompartment"),
                "initial_amount": get_if_set(
                    species,
                    "getInitialAmount",
                    "isSetInitialAmount",
                ),
                "initial_concentration": get_if_set(
                    species,
                    "getInitialConcentration",
                    "isSetInitialConcentration",
                ),
                "substance_units": get_if_set(
                    species,
                    "getSubstanceUnits",
                    "isSetSubstanceUnits",
                ),
                "has_only_substance_units": get_if_set(
                    species,
                    "getHasOnlySubstanceUnits",
                    "isSetHasOnlySubstanceUnits",
                ),
                "boundary_condition": get_if_set(
                    species,
                    "getBoundaryCondition",
                    "isSetBoundaryCondition",
                ),
                "constant": get_if_set(species, "getConstant", "isSetConstant"),
                "conversion_factor": get_if_set(
                    species,
                    "getConversionFactor",
                    "isSetConversionFactor",
                ),
            }
        )
        rows.append(row)

    return pd.DataFrame(rows)


def extract_parameters(model: libsbml.Model) -> pd.DataFrame:
    rows = []

    for parameter in model.getListOfParameters():
        row = common_fields(parameter)
        row.update(
            {
                "value": get_if_set(parameter, "getValue", "isSetValue"),
                "units": get_if_set(parameter, "getUnits", "isSetUnits"),
                "constant": get_if_set(parameter, "getConstant", "isSetConstant"),
            }
        )
        rows.append(row)

    return pd.DataFrame(rows)


def extract_initial_assignments(model: libsbml.Model) -> pd.DataFrame:
    rows = []

    for assignment in model.getListOfInitialAssignments():
        row = common_fields(assignment)
        row.update(
            {
                "symbol": get_if_set(assignment, "getSymbol", "isSetSymbol"),
                "math": math_to_string(assignment.getMath()),
            }
        )
        rows.append(row)

    return pd.DataFrame(rows)


def extract_rules(model: libsbml.Model) -> pd.DataFrame:
    rows = []

    for rule in model.getListOfRules():
        row = common_fields(rule)

        if rule.isAssignment():
            rule_type = "assignment"
        elif rule.isRate():
            rule_type = "rate"
        elif rule.isAlgebraic():
            rule_type = "algebraic"
        else:
            rule_type = "unknown"

        row.update(
            {
                "rule_type": rule_type,
                "variable": get_if_set(rule, "getVariable", "isSetVariable"),
                "math": math_to_string(rule.getMath()),
            }
        )
        rows.append(row)

    return pd.DataFrame(rows)


def extract_function_definitions(model: libsbml.Model) -> pd.DataFrame:
    rows = []

    for function_definition in model.getListOfFunctionDefinitions():
        row = common_fields(function_definition)
        row.update({"math": math_to_string(function_definition.getMath())})
        rows.append(row)

    return pd.DataFrame(rows)


def extract_reactions(model: libsbml.Model) -> Dict[str, pd.DataFrame]:
    reaction_rows = []
    participant_rows = []
    kinetic_law_rows = []
    local_parameter_rows = []

    for reaction in model.getListOfReactions():
        reaction_id = get_if_set(reaction, "getId", "isSetId")

        reaction_row = common_fields(reaction)
        reaction_row.update(
            {
                "reversible": get_if_set(reaction, "getReversible", "isSetReversible"),
                "fast": get_if_set(reaction, "getFast", "isSetFast"),
                "compartment": get_if_set(
                    reaction,
                    "getCompartment",
                    "isSetCompartment",
                ),
            }
        )
        reaction_rows.append(reaction_row)

        participant_rows.extend(extract_reaction_participants(reaction, reaction_id))

        if reaction.isSetKineticLaw():
            kinetic_law = reaction.getKineticLaw()

            kinetic_law_row = common_fields(kinetic_law)
            kinetic_law_row.update(
                {
                    "reaction_id": reaction_id,
                    "kinetic_law_math": math_to_string(kinetic_law.getMath()),
                    "substance_units": get_if_set(
                        kinetic_law,
                        "getSubstanceUnits",
                        "isSetSubstanceUnits",
                    ),
                    "time_units": get_if_set(
                        kinetic_law,
                        "getTimeUnits",
                        "isSetTimeUnits",
                    ),
                }
            )
            kinetic_law_rows.append(kinetic_law_row)

            local_parameter_rows.extend(
                extract_local_parameters(kinetic_law, reaction_id)
            )

    return {
        "reactions": pd.DataFrame(reaction_rows),
        "reaction_participants": pd.DataFrame(participant_rows),
        "kinetic_laws": pd.DataFrame(kinetic_law_rows),
        "local_parameters": pd.DataFrame(local_parameter_rows),
    }


def extract_reaction_participants(
    reaction: libsbml.Reaction,
    reaction_id: str,
) -> list[Dict[str, Any]]:
    rows = []

    participant_groups = [
        ("reactant", reaction.getListOfReactants()),
        ("product", reaction.getListOfProducts()),
        ("modifier", reaction.getListOfModifiers()),
    ]

    for role, species_refs in participant_groups:
        for species_ref in species_refs:
            row = common_fields(species_ref)
            row.update(
                {
                    "reaction_id": reaction_id,
                    "role": role,
                    "species": get_if_set(species_ref, "getSpecies", "isSetSpecies"),
                }
            )

            if role in {"reactant", "product"}:
                row.update(
                    {
                        "stoichiometry": get_if_set(
                            species_ref,
                            "getStoichiometry",
                            "isSetStoichiometry",
                        ),
                        "constant": get_if_set(
                            species_ref,
                            "getConstant",
                            "isSetConstant",
                        ),
                        "stoichiometry_math": (
                            math_to_string(species_ref.getStoichiometryMath().getMath())
                            if species_ref.isSetStoichiometryMath()
                            else ""
                        ),
                    }
                )
            else:
                row.update(
                    {
                        "stoichiometry": "",
                        "constant": "",
                        "stoichiometry_math": "",
                    }
                )

            rows.append(row)

    return rows


def extract_local_parameters(
    kinetic_law: libsbml.KineticLaw,
    reaction_id: str,
) -> list[Dict[str, Any]]:
    rows = []

    for parameter in kinetic_law.getListOfLocalParameters():
        row = common_fields(parameter)
        row.update(
            {
                "reaction_id": reaction_id,
                "value": get_if_set(parameter, "getValue", "isSetValue"),
                "units": get_if_set(parameter, "getUnits", "isSetUnits"),
                "parameter_scope": "local_parameter",
            }
        )
        rows.append(row)

    # Older SBML files may use local parameters as Parameters.
    if has_method(kinetic_law, "getListOfParameters"):
        for parameter in kinetic_law.getListOfParameters():
            row = common_fields(parameter)
            row.update(
                {
                    "reaction_id": reaction_id,
                    "value": get_if_set(parameter, "getValue", "isSetValue"),
                    "units": get_if_set(parameter, "getUnits", "isSetUnits"),
                    "constant": get_if_set(parameter, "getConstant", "isSetConstant"),
                    "parameter_scope": "kinetic_law_parameter",
                }
            )
            rows.append(row)

    return rows


def extract_events(model: libsbml.Model) -> Dict[str, pd.DataFrame]:
    event_rows = []
    event_assignment_rows = []

    for event in model.getListOfEvents():
        event_id = get_if_set(event, "getId", "isSetId")

        event_row = common_fields(event)
        event_row.update(
            {
                "use_values_from_trigger_time": get_if_set(
                    event,
                    "getUseValuesFromTriggerTime",
                    "isSetUseValuesFromTriggerTime",
                ),
                "trigger_math": (
                    math_to_string(event.getTrigger().getMath())
                    if event.isSetTrigger()
                    else ""
                ),
                "trigger_initial_value": (
                    get_if_set(event.getTrigger(), "getInitialValue", "isSetInitialValue")
                    if event.isSetTrigger()
                    else ""
                ),
                "trigger_persistent": (
                    get_if_set(event.getTrigger(), "getPersistent", "isSetPersistent")
                    if event.isSetTrigger()
                    else ""
                ),
                "delay_math": (
                    math_to_string(event.getDelay().getMath())
                    if event.isSetDelay()
                    else ""
                ),
                "priority_math": (
                    math_to_string(event.getPriority().getMath())
                    if event.isSetPriority()
                    else ""
                ),
            }
        )
        event_rows.append(event_row)

        for assignment in event.getListOfEventAssignments():
            assignment_row = common_fields(assignment)
            assignment_row.update(
                {
                    "event_id": event_id,
                    "variable": get_if_set(
                        assignment,
                        "getVariable",
                        "isSetVariable",
                    ),
                    "math": math_to_string(assignment.getMath()),
                }
            )
            event_assignment_rows.append(assignment_row)

    return {
        "events": pd.DataFrame(event_rows),
        "event_assignments": pd.DataFrame(event_assignment_rows),
    }


def extract_constraints(model: libsbml.Model) -> pd.DataFrame:
    rows = []

    for constraint in model.getListOfConstraints():
        row = common_fields(constraint)
        row.update(
            {
                "math": math_to_string(constraint.getMath()),
                "message": clean_notes(constraint.getMessageString() or ""),
            }
        )
        rows.append(row)

    return pd.DataFrame(rows)


# =============================================================================
# Excel writer
# =============================================================================

def export_tables_to_excel(
    tables: Dict[str, pd.DataFrame],
    output_excel_path: Path,
) -> None:
    """
    Write extracted SBML tables to an Excel workbook.
    """
    output_excel_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_excel_path, engine="openpyxl") as writer:
        sheet_name_map = {}

        for table_name, df in tables.items():
            sheet_name = make_excel_sheet_name(table_name, existing=set(sheet_name_map))
            sheet_name_map[table_name] = sheet_name
            df.to_excel(writer, sheet_name=sheet_name, index=False)

        autosize_excel_columns(writer, tables, sheet_name_map)


def make_excel_sheet_name(table_name: str, existing: set[str]) -> str:
    """
    Excel sheet names must be <= 31 characters and unique.
    """
    base = table_name[:31]
    candidate = base
    counter = 1

    while candidate in existing:
        suffix = f"_{counter}"
        candidate = f"{base[:31 - len(suffix)]}{suffix}"
        counter += 1

    return candidate


def autosize_excel_columns(
    writer: pd.ExcelWriter,
    tables: Dict[str, pd.DataFrame],
    sheet_name_map: Dict[str, str],
) -> None:
    """
    Make the Excel output easier to inspect.
    """
    for table_name, df in tables.items():
        sheet_name = sheet_name_map[table_name]
        worksheet = writer.sheets[sheet_name]

        if df.empty:
            worksheet.cell(row=1, column=1, value="No entries found")
            continue

        for col_idx, column in enumerate(df.columns, start=1):
            values = df[column].head(200).fillna("").astype(str).tolist()
            max_len = max([len(str(column))] + [len(value) for value in values])
            width = min(max(max_len + 2, 12), 80)

            column_letter = worksheet.cell(row=1, column=col_idx).column_letter
            worksheet.column_dimensions[column_letter].width = width

        worksheet.freeze_panes = "A2"
