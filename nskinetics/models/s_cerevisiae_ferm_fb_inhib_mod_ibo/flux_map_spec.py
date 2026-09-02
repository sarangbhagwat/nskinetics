# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Flux-map layout for the shipped *S. cerevisiae* ethanol/isobutanol model.

Node placement follows ``conceptual_diagram.py``, re-fitted to one 88 mm-wide
panel (the Ehrlich column is raised so its entry edge is not horizontal, and
the network is shifted down to sit under the panel header). Edges are the
kinetic reactions (dilution/outflow reactions are omitted -- D = 0 in fed-batch use).
The inhibition mapping is the canonical set of product-inhibition coefficients:
exponential terms (k_*ie/ia/ii) and the saturable denominator self-inhibition
terms (K_6e on r6, K_16i on r16). r10 (active-biomass decay) is
product-ENHANCED, listed in enhancement_reactions.

r10 is deliberately absent from ``edges``: biomass decay has no network edge to
hang a strip on, so no strip is drawn for it. It stays in ``inhibition_map`` so
that a caller who passes ``reactions=None`` (or a list including ``'r10'``) to
:func:`~nskinetics.compute_flux_summary` still gets its enhancement numbers.
"""

from ...visualization import FluxMapSpec

__all__ = ('FLUX_MAP_SPEC',)

# (x_mm, y_mm, label) within an 88 x 92 mm panel. The top ~10 mm carry the
# panel header, so the network sits below y = 82. The Ehrlich column (right)
# starts level with glucose so that r13 leaves pyruvate as a steep diagonal:
# ``draw_flux_map`` puts a flux label 1.5 mm to the left of an edge midpoint,
# which lands *on* a horizontal edge.
_NODES = {
    'glu': (44, 78, 'Glucose'),
    'pyr': (44, 56, 'Pyruvate'),
    'ald': (44, 32, 'Acetald.'),
    'eth': (44, 8, 'Ethanol'),
    'ace': (22, 32, 'Acetate'),
    'tca': (12, 56, 'TCA'),
    'bm':  (12, 8, 'Biomass'),
    'al':  (72, 78, 'AL'),
    'dhi': (72, 54, 'DHIV'),
    'kiv': (72, 30, 'KIV'),
    'ibo': (72, 8, 'Isobutanol'),
}
_EDGES = {
    'r1': ('glu', 'pyr'),   # glycolysis
    'r3': ('pyr', 'ald'),   # PDC
    'r6': ('ald', 'eth'),   # ADH
    'r2': ('pyr', 'tca'),   # pyruvate -> TCA
    'r4': ('ald', 'ace'),   # acetaldehyde -> acetate
    'r5': ('ace', 'tca'),   # acetate -> TCA
    'r7': ('glu', 'bm'),    # growth on glucose
    'r8': ('ace', 'bm'),    # growth on acetate
    'r13': ('pyr', 'al'),   # ALS
    'r14': ('al', 'dhi'),   # KARI
    'r15': ('dhi', 'kiv'),  # DHAD
    'r16': ('kiv', 'ibo'),  # KDC+ADH
}
# {inhibitor: color} (Okabe-Ito, matching conceptual_diagram)
_INHIBITORS = {
    'ethanol': '#D55E00',
    'acetate': '#CC79A7',
    'isobutanol': '#0072B2',
}
# {parameter_id: (reaction_id, inhibitor_name)}
_INHIBITION_MAP = {
    'k_1ie': ('r1', 'ethanol'),   'k_1ia': ('r1', 'acetate'),   'k_1ii': ('r1', 'isobutanol'),
    'k_4ie': ('r4', 'ethanol'),   'k_4ia': ('r4', 'acetate'),   'k_4ii': ('r4', 'isobutanol'),
    'K_6e':  ('r6', 'ethanol'),   'k_6ia': ('r6', 'acetate'),   'k_6ii': ('r6', 'isobutanol'),
    'k_7ie': ('r7', 'ethanol'),   'k_7ia': ('r7', 'acetate'),   'k_7ii': ('r7', 'isobutanol'),
    'k_10ie': ('r10', 'ethanol'), 'k_10ia': ('r10', 'acetate'), 'k_10ii': ('r10', 'isobutanol'),
    'k_16ie': ('r16', 'ethanol'), 'k_16ia': ('r16', 'acetate'), 'K_16i': ('r16', 'isobutanol'),
}

#: Shipped :class:`~nskinetics.visualization.FluxMapSpec` for this model; its
#: ``inhibition_map`` is the canonical coefficient mapping to hand to
#: :func:`~nskinetics.compute_flux_summary`.
FLUX_MAP_SPEC = FluxMapSpec(
    nodes=_NODES,
    edges=_EDGES,
    inhibitors=_INHIBITORS,
    inhibition_map=_INHIBITION_MAP,
    enhancement_reactions={'r10'},
    products=['s_EtOH', 's_IBO'],
    # Strips hang down from their anchor. Only inhibited reactions get one, so
    # r2/r5/r8 need no entry; r1/r4/r7 are moved off the diagonals, node boxes
    # and flux labels they would otherwise cross.
    strip_offsets={'r1': (-19.0, 0.0), 'r4': (-5.0, -5.0), 'r7': (5.0, 8.0)},
    size_mm=(88.0, 92.0),
)
