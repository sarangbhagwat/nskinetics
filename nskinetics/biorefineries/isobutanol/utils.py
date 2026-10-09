# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Helpers for authoring the parameter-distribution workbooks.

:func:`generate_save_kinetic_parameter_distributions` writes the kinetic rows
of a parameter-distribution table (the format of the workbooks under
``data/parameter_distributions``) for every ``k_*``/``K_*`` parameter of a
kinetic model.
"""

import csv

__all__ = ('generate_save_kinetic_parameter_distributions',)


def generate_save_kinetic_parameter_distributions(nsk_kinetic_model, filename,
                                                  fermentation_reactor_ID='V406',
                                                  shape='triangular', factors=(0.8, 1.0, 1.2)):
    """Write a CSV of triangular distributions for every kinetic parameter.

    One row per global parameter of ``nsk_kinetic_model._te`` whose ID starts
    with ``k_``/``K_`` (case-insensitive), with lower/midpoint/upper values
    of ``baseline * factors`` and a workbook load statement
    ``'<fermentation_reactor_ID>.nsk_kinetic_model._te.<ID> = x'``.

    Parameters
    ----------
    nsk_kinetic_model : nskinetics.KineticModel
        Kinetic model whose RoadRunner instance (``_te``) is read.
    filename : str
        Output CSV path (relative paths resolve against the working
        directory).
    fermentation_reactor_ID : str, optional
        Flowsheet ID of the fermentor in the load statements. Default 'V406'.
    shape : str, optional
        Distribution shape written to the 'Shape' column. Default
        'triangular'.
    factors : tuple[float, float, float], optional
        Multipliers of the baseline for the lower, midpoint and upper values.
        Default (0.8, 1.0, 1.2), the +-20 % kinetic convention.

    Notes
    -----
    The 'Units' column is a heuristic: any ID containing a lowercase ``k`` is
    labelled a rate (``g_per_l_per_h``), else ``g_per_l``. Check the units
    before using the rows in a workbook.
    """
    r_te = nsk_kinetic_model._te
    all_params = r_te.getGlobalParameterIds()
    kinetic_params = [i for i in all_params if i[:2].lower()=='k_']

    data = []
    data.append(['Parameter name', 'Element', 'Kind', 'Units',
                 'Baseline', 'Shape', 'Lower', 'Midpoint', 'Upper',
                 'References', 'Load statement'])

    for k in kinetic_params:
        if 'k' in k:
            units = 'g_per_l_per_h'
        elif 'K' in k:
            units = 'g_per_l'
        else:
            units = None
        baseline = r_te.__getattribute__(k)
        data.append([k, 'Fermentation Kinetics', 'coupled', units,
                     baseline, shape, baseline*factors[0], baseline*factors[1], baseline*factors[2],
                     None, f'{fermentation_reactor_ID}.nsk_kinetic_model._te.{k} = x'])

    with open(filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
