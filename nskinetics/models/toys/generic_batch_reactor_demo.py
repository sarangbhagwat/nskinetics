# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

"""
Minimal, non-isobutanol demonstration that ``NSKBatchReactor`` is general:
a generic Tellurium model with a single substrate consumed to product, plus one
``FeedSpike``, driving a biosteam batch reactor. Run directly to build the
reactor, simulate it end-to-end (kinetics -> effluent -> biosteam sizing/cost),
and print a short summary.
"""
import tellurium as te

import biosteam as bst
import nskinetics as nsk

__all__ = ('build_demo_reactor',)

_MODEL = """
model generic_batch()
  compartment env; species S in env, P in env;
  S = 50; P = 0; env = 1; k = 0.5;
  threshold = 5; target = 50; feed_conc = 500;
  n_spk = 0; n_max = 3; last_vol = 0; tot_vol = 0;
  n_spk has dimensionless; n_max has dimensionless;
  J: S => P; k*S*env;
  curr_env := env;
end
"""


def build_demo_reactor():
    """Build (do not simulate) an ``NSKBatchReactor`` on a generic kinetic model.

    Returns
    -------
    (reactor, te_r) : tuple
    """
    r = te.loadAntimonyModel(_MODEL)
    te_r = nsk.KineticModel(r, units={'time': 'h', 'conc': 'g/L'})
    spike = nsk.FeedSpike(
        species='S', when='S <= threshold', target='target',
        feed_conc='feed_conc', volume_var='env', max_count='n_max',
        count_var='n_spk', last_vol_var='last_vol', tot_vol_var='tot_vol',
        priority=5, name='S_spike')
    for event in spike.expand():
        te_r.add_event(event)
    te_r.compile_events()

    bst.settings.set_thermo(['Water', 'Glucose', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=1000., Glucose=50., units='kg/hr')
    seed = bst.Stream('seed', Water=10., units='kg/hr')
    spike_feed = bst.Stream('spike_feed', Glucose=500., units='kg/hr')

    reactor = nsk.units.NSKBatchReactor(
        'R_demo', ins=(feed, seed, spike_feed), outs=('vent', 'effluent'),
        nsk_kinetic_model=te_r,
        map_species_to_chemicals={'S': 'Glucose', 'P': 'Ethanol'},
        track_vars=['[S]', '[P]'],  # record concentrations for plotting
        tau=48., tau_max=120., volume_var='curr_env',
        feed_volume_added_var='tot_vol',  # exercise the fed-batch volume correction
        spike_feed_index=2)
    return reactor, te_r


if __name__ == '__main__':
    import warnings

    reactor, te_r = build_demo_reactor()

    # Standalone `simulate()` runs biosteam's design/cost over several internal
    # passes and does not re-clear results between them, so a later pass sees the
    # previous pass's costs and trips biosteam's "_run added unit results" check.
    # The kinetic `_run` adds no costs (early passes are clean), so this warning
    # is a benign standalone-simulation artifact; silence just that one message.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore', message=r'.*method added unit results.*',
            category=RuntimeWarning)
        reactor.simulate()

    reactor.show(N=100)
    print(reactor.results())
