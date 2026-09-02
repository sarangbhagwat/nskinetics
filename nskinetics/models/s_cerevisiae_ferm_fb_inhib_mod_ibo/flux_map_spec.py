# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.
"""Flux-map layout for the shipped *S. cerevisiae* ethanol/isobutanol model.

Node placement follows ``conceptual_diagram.py``, re-fitted to one 88 x 92 mm
panel and re-arranged so no edge crosses another: the central glycolytic
column runs down x = 44, its two by-product branches leave to the left
(biomass at the top, acetate/TCA in the middle), and the engineered Ehrlich
column runs down x = 72. Edges are the kinetic reactions (dilution/outflow
reactions are omitted -- D = 0 in fed-batch use). The inhibition mapping is
the canonical set of product-inhibition coefficients: exponential terms
(k_*ie/ia/ii) and the saturable denominator self-inhibition terms (K_6e on r6,
K_16i on r16). r10 (active-biomass decay) is product-ENHANCED, listed in
enhancement_reactions.

r10 is deliberately absent from ``edges``: biomass decay has no network edge to
hang a strip on, so no strip is drawn for it. It stays in ``inhibition_map``,
and :attr:`~nskinetics.visualization.FluxMapSpec.reactions` appends it after
the drawn edges, so the shipped path (:func:`draw_scenario_flux_map`) still
computes its enhancement numbers.
"""

from ...visualization import FluxMapSpec

__all__ = ('FLUX_MAP_SPEC', 'draw_scenario_flux_map')

# (x_mm, y_mm, label) within an 88 x 92 mm panel. The top ~8 mm carry the
# panel header, so the network sits below y = 82. Value labels go above a
# horizontal edge / right of a vertical one, and strips below / left, so the
# left-hand branches (r7 at y = 78, r4 at y = 32) keep clear space on their
# lower-left.
_NODES = {
    'glu': (44, 78, 'Glucose'),
    'pyr': (44, 56, 'Pyruvate'),
    'ald': (44, 32, 'Acetald.'),
    'eth': (44, 8, 'Ethanol'),
    'bm':  (12, 78, 'Biomass'),
    'tca': (26, 56, 'TCA'),
    'ace': (12, 32, 'Acetate'),
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
#: :func:`~nskinetics.compute_flux_summary`, and its ``reactions`` property is
#: the reaction list that path should compute (drawn edges plus r10).
FLUX_MAP_SPEC = FluxMapSpec(
    nodes=_NODES,
    edges=_EDGES,
    inhibitors=_INHIBITORS,
    inhibition_map=_INHIBITION_MAP,
    enhancement_reactions={'r10'},
    products=['s_EtOH', 's_IBO'],
    product_labels={'s_EtOH': 'ethanol', 's_IBO': 'isobutanol'},
    # Strips hang down from their anchor (4 rows, ~5.3 mm tall). Only
    # inhibited reactions get one, so r2/r5/r8/r13-r15 need no entry; the
    # entries below move r1/r4/r6/r7/r16 off the node boxes, edges and value
    # labels they would otherwise cross.
    strip_offsets={
        'r7': (-8.5, -2.5),
        'r1': (-13.0, 2.5),
        'r4': (-4.5, -2.5),
        'r6': (-13.0, 2.5),
        'r16': (-13.0, 2.5),
    },
    size_mm=(88.0, 92.0),
)


def draw_scenario_flux_map(unit, save_dir=None, spec=FLUX_MAP_SPEC,
                           **draw_kwargs):
    """Simulate both scenarios on ``unit`` and draw the two-panel flux map.

    Applies scenario A to the shipped ``te_r``, re-simulates ``unit.system``,
    summarizes; repeats for scenario B; then restores scenario A and
    re-simulates in a ``finally`` so the caller's system is left as it was
    found. Only the kinetic Ehrlich rate constants change between panels: the
    fed-batch feeding strategy (spike count, thresholds, tau_max) is NOT
    touched, so the two panels differ solely by the engineered pathway.

    Parameters
    ----------
    unit : NSKBatchReactor
        A reactor in an already-built system (e.g. ``V406``); its
        ``system.simulate()`` is called three times.
    save_dir : str, optional
        Passed to :func:`~nskinetics.visualization.draw_flux_map`; writes
        ``flux_map.png`` / ``flux_map.pdf`` there.
    spec : FluxMapSpec, optional
        Layout and inhibition mapping; defaults to :data:`FLUX_MAP_SPEC`.
        ``spec.reactions`` sets the reactions summarized, so mapped-but-undrawn
        reactions (r10) are computed too.
    **draw_kwargs
        Forwarded to :func:`~nskinetics.visualization.draw_flux_map`
        (``formats``, ``dpi``, ``annotate_values``, ``show``).

    Returns
    -------
    (matplotlib.figure.Figure, list of matplotlib.axes.Axes, tuple)
        The figure, the panel axes, and ``(summary_A, summary_B)``.
    """
    from ...engine import compute_flux_summary
    from ...visualization import draw_flux_map
    from .s_cerevisiae_ferm_fb_inhib_mod_ibo import te_r
    from .scenarios import apply_scenario_A, apply_scenario_B

    imap = spec.inhibition_map
    reactions = spec.reactions
    try:
        apply_scenario_A(te_r)
        unit.system.simulate()
        summary_A = compute_flux_summary(unit, imap, reactions=reactions,
                                         label='Scenario A')
        apply_scenario_B(te_r)
        unit.system.simulate()
        summary_B = compute_flux_summary(unit, imap, reactions=reactions,
                                         label='Scenario B')
    finally:
        apply_scenario_A(te_r)
        unit.system.simulate()
    fig, axes = draw_flux_map([summary_A, summary_B], spec,
                              save_dir=save_dir, **draw_kwargs)
    return fig, axes, (summary_A, summary_B)
