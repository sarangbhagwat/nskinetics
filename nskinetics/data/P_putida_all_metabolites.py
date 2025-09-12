# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 15:31:10 2025

@author: sarangbhagwat

Extract and clean compound–SMILES pairs from a tab-separated file.

- Reads a tab-separated text file with two columns: compound identifier/name and SMILES.
- Cleans whitespace, drops empty rows, and de-duplicates.
- Validates and canonicalizes SMILES with RDKit if available.
- Computes InChIKey if RDKit + InChI support are present.
- Saves a clean rd_inchi next to the source file.

"""
from __future__ import annotations
import os
import csv
from pathlib import Path
from typing import Optional

import pandas as pd

from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi  # InChI support is optional in some builds
    
from ugropy import Groups

import thermosteam as tmo

def load_putida_smiles(path: Path) -> pd.DataFrame:
    """
    Load the putida SMILES table.

    Parameters
    ----------
    path : Path
        Path to the input TXT file (tab-separated).

    Returns
    -------
    DataFrame
        Columns: ['compound', 'smiles'] plus optional ['smiles_canonical', 'inchi_key', 'valid_smiles']
    """
    # Read as tab-separated; keep raw strings
    df = pd.read_csv(
        path,
        sep="\t",
        header=0,
        dtype=str,
        quoting=csv.QUOTE_NONE,
        engine="python",  # tolerant to odd characters
    )

    # Normalize column names safely
    # Expect something like: "(Compounds Polynucleotides)" and "SMILES"
    df.columns = [str(c).strip().lower().replace(" ", "_") for c in df.columns]

    # Map to standard names when possible
    col_map = {}
    for c in df.columns:
        if "smile" in c:
            col_map[c] = "smiles"
        elif "compound" in c or "name" in c or "id" in c:
            col_map[c] = "compound"
    if col_map:
        df = df.rename(columns=col_map)

    # If columns aren’t exactly two after rename, try best-effort selection
    if "compound" not in df.columns or "smiles" not in df.columns:
        # Fall back: assume first col = compound, last col = smiles
        df = df.iloc[:, [0, -1]]
        df.columns = ["compound", "smiles"]

    # Trim whitespace
    for col in ["compound", "smiles"]:
        df[col] = df[col].astype(str).str.strip()

    # Drop rows missing SMILES or compound
    df = df.replace({"": pd.NA}).dropna(subset=["compound", "smiles"])

    # Drop exact duplicates
    df = df.drop_duplicates(subset=["compound", "smiles"]).reset_index(drop=True)

    # Optional: Validate and canonicalize SMILES with RDKit
    valid_flags = []
    canon_smiles = []
    inchi_keys = []
    for s in df["smiles"]:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            valid_flags.append(False)
            canon_smiles.append(pd.NA)
            inchi_keys.append(pd.NA)
        else:
            valid_flags.append(True)
            # Canonical (Kekule/implicit H standard form)
            canon_smiles.append(Chem.MolToSmiles(mol, canonical=True))
            # InChIKey if supported
            if rd_inchi is not None:
                try:
                    inchi_keys.append(rd_inchi.MolToInchiKey(mol))
                except Exception:
                    inchi_keys.append(pd.NA)
            else:
                inchi_keys.append(pd.NA)
    df["valid_smiles"] = valid_flags
    df["smiles_canonical"] = canon_smiles
    df["inchi_key"] = inchi_keys

    return df


def load_putida_chems(path: Path):
    """
    Load thermosteam-compatible Putida metabolite chemicals.
    """
    df = load_putida_smiles(path)
    
    DortmundGroupCounts = tmo.equilibrium.unifac.DortmundGroupCounts
    chems = tmo.Chemicals([])
    n_smiles_imported = 0
    smiles_not_loaded = []
    chems_with_all_props = tmo.Chemicals([])
    ugropy_found_groups = {}
    for i, smiles in zip(df['valid_smiles'], df['smiles_canonical']): 
        if i: 
            n_smiles_imported+=1
            try: 
                c = tmo.Chemical(f'SMILES={smiles}')
                chems.append(c)
                if c.Dortmund and (c.Tb is not None) and c.Hvap and c.Psat:
                    try:
                        c = tmo.Chemical(ID=c.iupac_name[0], #!!! if single-char ID error, check iupac_name
                                         search_ID=f'SMILES={smiles}')
                    except:
                        breakpoint()
                    chems_with_all_props.append(c)
                elif (not c.Dortmund) and (c.Tb is not None) and c.Hvap and c.Psat: # only Dortmund missing
                    # Access Dortmund groups using ugropy
                    print(f'\n\nTrying ugropy for {smiles} ...\n')
                    molecule_groups = Groups(smiles, 'smiles')
                    dortmund_groups = molecule_groups.dortmund.subgroups
                    if dortmund_groups: 
                        print(f'Found groups: {dortmund_groups}')
                        ugropy_found_groups[smiles] = dortmund_groups
                        # breakpoint()
                        reformatted_dortmund_groups = {v: k for k,v in dortmund_groups.items()}
                        c._Dortmund = DortmundGroupCounts.from_dict(reformatted_dortmund_groups)
                        try:
                            c = tmo.Chemical(ID=c.iupac_name[0], search_ID=f'SMILES={smiles}')
                        except:
                            breakpoint()
                        chems_with_all_props.append(c)
            except:
                smiles_not_loaded.append(smiles)
    
    print(ugropy_found_groups)
    
    return chems_with_all_props