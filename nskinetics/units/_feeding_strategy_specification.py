# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the MIT open-source license. See
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

from flexsolve import IQ_interpolation

from ..exceptions import FeedingStrategyError

__all__ = ('FedBatchStrategySpecification', 'SpikeControlVariables',
           'ConcentrationActuator')

#%%
class ConcentrationActuator:
    """
    One adjustable knob on an upstream unit that a
    :class:`FedBatchStrategySpecification` may tune to reach a desired
    concentration in a terminal stream.

    Parameters
    ----------
    unit : biosteam.Unit
        The unit operation carrying the adjustable attribute.
    attr : str
        Name of the attribute to set (e.g. ``'V'`` on an evaporator,
        ``'water_to_sugar_mol_ratio'`` on a mixer).
    lb : float
        Lower bound of the solve bracket used when tuning this actuator.
    ub : float
        Upper bound of the solve bracket used when tuning this actuator.
    """
    def __init__(self, unit, attr, lb, ub):
        self.unit = unit
        self.attr = attr
        self.lb = lb
        self.ub = ub

    def set(self, value):
        """Set the actuated attribute to ``value``."""
        setattr(self.unit, self.attr, value)

    def get(self):
        """Return the current value of the actuated attribute."""
        return getattr(self.unit, self.attr)

    def __repr__(self):
        unit_ID = getattr(self.unit, 'ID', self.unit)
        return (f'{self.__class__.__name__}(unit={unit_ID}, attr={self.attr!r}, '
                f'lb={self.lb}, ub={self.ub})')


#%%
class SpikeControlVariables:
    """
    Names of the kinetic-model variables and result columns through which a
    :class:`FedBatchStrategySpecification` controls fed-batch spike feeding.

    Parameters
    ----------
    spike_conc_var : str
        Model variable holding the spike-feed concentration.
    target_conc_var : str
        Model variable holding the concentration a spike restores.
    threshold_conc_var : str
        Model variable holding the concentration that triggers a spike.
    volume_col : str, optional
        Result column for the working volume. Defaults to the fermentation
        reactor's ``volume_var``.
    feed_volume_added_col : str, optional
        Result column for the cumulative spike-feed volume added. Defaults to
        the fermentation reactor's ``feed_volume_added_var``.
    max_n_spikes_var : str, optional
        Model variable capping the number of feed spikes. Defaults to the
        fermentation reactor's ``spike_retry.max_count_var``.
    default_max_n_spikes_attr : str, optional
        Attribute on the kinetic model holding the cap its reset function
        restores at the start of every reactor run. Defaults to
        ``f'default_{max_n_spikes_var}'``.
    """
    def __init__(self, *, spike_conc_var, target_conc_var, threshold_conc_var,
                 volume_col=None, feed_volume_added_col=None,
                 max_n_spikes_var=None, default_max_n_spikes_attr=None):
        self.spike_conc_var = spike_conc_var
        self.target_conc_var = target_conc_var
        self.threshold_conc_var = threshold_conc_var
        self.volume_col = volume_col
        self.feed_volume_added_col = feed_volume_added_col
        self.max_n_spikes_var = max_n_spikes_var
        self.default_max_n_spikes_attr = default_max_n_spikes_attr

    def _resolve(self, explicit, reactor, reactor_attr, what):
        if explicit is not None:
            return explicit
        col = getattr(reactor, reactor_attr, None)
        if col is None:
            raise ValueError(
                f'No {what} result column: none was given and the '
                f'fermentation reactor declares no {reactor_attr!r}.')
        return col

    def resolve_volume_col(self, reactor):
        """Working-volume result column (explicit, else ``reactor.volume_var``)."""
        return self._resolve(self.volume_col, reactor, 'volume_var',
                             'working-volume')

    def resolve_feed_volume_added_col(self, reactor):
        """Cumulative-feed-volume column (explicit, else
        ``reactor.feed_volume_added_var``)."""
        return self._resolve(self.feed_volume_added_col, reactor,
                             'feed_volume_added_var', 'feed-volume-added')

    def resolve_max_n_spikes_var(self, reactor):
        """Spike-cap model variable.

        Parameters
        ----------
        reactor : NSKBatchReactor
            Reactor supplying the fallback ``spike_retry.max_count_var``.

        Returns
        -------
        str
            The explicit ``max_n_spikes_var`` when given, else the reactor's
            ``spike_retry.max_count_var``.
        """
        if self.max_n_spikes_var is not None:
            return self.max_n_spikes_var
        spike_retry = getattr(reactor, 'spike_retry', None)
        var = getattr(spike_retry, 'max_count_var', None)
        if var is None:
            raise ValueError(
                'No spike-cap model variable: none was given and the '
                'fermentation reactor declares no spike_retry.max_count_var.')
        return var

    def resolve_default_max_n_spikes_attr(self, reactor):
        """Kinetic-model attribute holding the spike cap restored on reset.

        Parameters
        ----------
        reactor : NSKBatchReactor
            Reactor used to resolve the underlying model variable name.

        Returns
        -------
        str
            The explicit ``default_max_n_spikes_attr`` when given, else
            ``f'default_{resolve_max_n_spikes_var(reactor)}'``.
        """
        if self.default_max_n_spikes_attr is not None:
            return self.default_max_n_spikes_attr
        return f'default_{self.resolve_max_n_spikes_var(reactor)}'


#%%
class FedBatchStrategySpecification:
    """
    Specification and orchestration of a fed-batch feeding strategy based on
    concentration control of selected species for fermentation processes.

    Defines target and threshold concentrations and translates them into
    operating conditions for upstream unit operations that supply both a base
    feed and a spike feed to a fermentation reactor. Chemical- and
    kinetic-model-agnostic: the controlled species, the solvent used as the
    concentration volume basis, the kinetic-model variable names, and the
    adjustable upstream knobs are all injected via parameters.

    The specification logic includes:

    - Setting desired concentrations for the initial feed and spike feed by
      tuning one concentrator and one diluter actuator per feed train.
    - Defining a threshold concentration that triggers spike feeding.
    - Updating the maximum reactor residence time (``tau_max``) and
      propagating its effects upstream through splitter and feed unit
      simulations.
    - Coordinated simulation of all upstream units to ensure consistency
      between specified targets and achievable process conditions.

    Parameters
    ----------
    target_conc : float
        Desired concentration of the controlled species in the initial feed
        to the fermentation reactor.
    threshold_conc : float
        Concentration in the reactor environment that triggers a spike
        feeding event.
    spike_conc : float
        Concentration of the spike feed stream used during
        threshold-triggered feeding.
    tau_max : float
        Maximum residence time of the fermentation reactor.
    max_n_spikes : float, optional
        Cap on the number of fed-batch spikes. ``None`` (default) leaves the
        kinetic model's own cap untouched; ``0`` means no spikes. Imposed as an
        upper bound on the first simulated run — a reactor's
        :class:`~nskinetics.units.SpikeReduceRetry` may reduce it further
        within that run.
    fermentation_reactor : NSKBatchReactor
        Fermentation reactor unit; must expose a kinetic model via its
        ``nsk_kinetic_model`` attribute and support simulation with
        spike feeding logic.
    splitter : biosteam.Unit
        Splitter dividing flow between the base-feed and spike-feed trains.
    control_variables : SpikeControlVariables
        Kinetic-model variable and result-column names through which the
        strategy is imposed on the model.
    feed_concentrator, feed_diluter : ConcentrationActuator
        Adjustable knobs on the initial-feed train (concentrating and
        diluting, respectively).
    spike_concentrator, spike_diluter : ConcentrationActuator
        Adjustable knobs on the spike-feed train.
    species_IDs : list of str
        biosteam chemical IDs of the controlled species; their total mass
        defines the controlled concentration.
    solvent_ID : str, optional
        Chemical ID whose stream volume is the concentration basis
        (default ``'Water'``).
    feed_units_sequential : list of biosteam.Unit, optional
        Units simulated in order to propagate initial-feed changes.
    spike_units_sequential : list of biosteam.Unit, optional
        Units simulated in order to propagate spike-feed changes.
    baseline_specifications : dict, optional
        Baseline values of the five specifications, keyed by
        ``'target_conc'``, ``'threshold_conc'``, ``'spike_conc'``,
        ``'tau_max'``, ``'max_n_spikes'``.

    Notes
    -----
    Set fed-batch strategy parameters through this specification rather than by
    writing kinetic-model variables directly. Three entry points exist, in
    descending order of precedence:

    1. a model specification — keyword arguments passed to
       :meth:`load_specifications`;
    2. the system specification a system factory attaches, which forwards its
       keyword arguments to :meth:`load_specifications`;
    3. the values stored on this object (e.g. ``spec.max_n_spikes = 16``).

    An explicit keyword argument at 1 or 2 always wins, because
    :meth:`load_specifications` falls back to a stored value only when the
    corresponding argument is ``None``.
    """

    def __init__(
        self,
        target_conc: float,
        threshold_conc: float,
        spike_conc: float,
        tau_max: float,
        fermentation_reactor,
        splitter,
        control_variables,
        feed_concentrator,
        feed_diluter,
        spike_concentrator,
        spike_diluter,
        species_IDs,
        solvent_ID='Water',
        feed_units_sequential=None,
        spike_units_sequential=None,
        baseline_specifications=None,
        max_n_spikes=None,
        ):
        self.target_conc = target_conc
        self.threshold_conc = threshold_conc
        self.spike_conc = spike_conc
        self.tau_max = tau_max
        self.max_n_spikes = max_n_spikes

        self._spec_names = ['target_conc',
                            'threshold_conc',
                            'spike_conc',
                            'tau_max',
                            'max_n_spikes',
                            ]

        self._baseline_specifications = baseline_specifications

        self.fermentation_reactor = fermentation_reactor
        self.splitter = splitter
        self.control_variables = control_variables
        self.feed_concentrator = feed_concentrator
        self.feed_diluter = feed_diluter
        self.spike_concentrator = spike_concentrator
        self.spike_diluter = spike_diluter
        self.feed_units_sequential = feed_units_sequential
        self.spike_units_sequential = spike_units_sequential

        self.species_IDs = list(species_IDs)
        self.solvent_ID = solvent_ID

        self._validate_parameters()

    def _validate_parameters(self) -> None:
        """
        Basic validation for initialization parameters.
        """
        if self.target_conc <= 0:
            raise ValueError("target_conc must be positive.")

        if self.threshold_conc < 0:
            raise ValueError("threshold_conc must be positive.")

        if self.spike_conc <= 0:
            raise ValueError("spike_conc must be positive.")

        if self.max_n_spikes is not None and self.max_n_spikes < 0:
            raise ValueError("max_n_spikes must be non-negative.")

        if self.threshold_conc > self.target_conc:
            raise ValueError("threshold_conc should not exceed target_conc.")

        if self.target_conc > self.spike_conc:
            raise ValueError("target_conc should not exceed spike_conc.")

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"target_conc={self.target_conc}, "
            f"threshold_conc={self.threshold_conc}, "
            f"spike_conc={self.spike_conc})"
        )

    def load_specifications(self,
                            target_conc=None,
                            spike_conc=None,
                            threshold_conc=None,
                            tau_max=None,
                            max_n_spikes=None,
                            ):
        if target_conc is None:
            target_conc = self.target_conc

        if threshold_conc is None:
            threshold_conc = self.threshold_conc

        if spike_conc is None:
            spike_conc = self.spike_conc

        if tau_max is None:
            tau_max = self.tau_max

        if max_n_spikes is None:
            max_n_spikes = self.max_n_spikes

        if not (threshold_conc < target_conc and target_conc < spike_conc):
            raise ValueError(f'Specifications do not meet required condition: threshold_conc ({threshold_conc}) < target_conc ({target_conc}) < spike_conc ({spike_conc}).\n')

        if max_n_spikes is not None and max_n_spikes < 0:
            raise ValueError(f'Specification max_n_spikes ({max_n_spikes}) must be non-negative.\n')

        self.threshold_conc = threshold_conc
        self.target_conc = target_conc
        self.spike_conc = spike_conc
        self.max_n_spikes = max_n_spikes

        # Imposed BEFORE the concentration actuators and
        # load_threshold_conc_and_tau_max: the latter simulates the reactor and
        # derives the splitter split from the spike volume that run added, so
        # that run must already honor the cap.
        self.load_max_n_spikes(max_n_spikes)

        self.load_desired_concs(target_conc=target_conc,
                                spike_conc=spike_conc)

        self.load_threshold_conc_and_tau_max(threshold_conc=threshold_conc,
                                             tau_max=tau_max)

    def load_max_n_spikes(self, max_n_spikes):
        """Impose the fed-batch spike cap on the kinetic model.

        Writes both places the cap must live: the live model variable the spike
        events read, and the kinetic model attribute its reset function restores
        at the start of every reactor run. Writing only the former is silently
        undone by that reset.

        Parameters
        ----------
        max_n_spikes : float or None
            The cap. ``None`` is a no-op — the model keeps the cap it has.
        """
        if max_n_spikes is None:
            return
        reactor = self.fermentation_reactor
        km = reactor.nsk_kinetic_model
        cv = self.control_variables
        km.set_value(cv.resolve_max_n_spikes_var(reactor), max_n_spikes)
        setattr(km, cv.resolve_default_max_n_spikes_attr(reactor), max_n_spikes)

    def load_desired_concs(self, target_conc, spike_conc):
        km = self.fermentation_reactor.nsk_kinetic_model
        cv = self.control_variables
        km.set_value(cv.spike_conc_var, spike_conc)
        km.set_value(cv.target_conc_var, target_conc)

        fermentation_reactor = self.fermentation_reactor
        spike_index = getattr(fermentation_reactor, 'spike_feed_index', None)
        if spike_index is None:
            spike_index = 2

        self._estimate_set_actuators(desired_conc=target_conc,
                               concentrator=self.feed_concentrator,
                               diluter=self.feed_diluter,
                               f_simulate_units_sequential=self._simulate_feed_units,
                               terminal_stream=fermentation_reactor.ins[0])

        self._estimate_set_actuators(desired_conc=spike_conc,
                               concentrator=self.spike_concentrator,
                               diluter=self.spike_diluter,
                               f_simulate_units_sequential=self._simulate_spike_units,
                               terminal_stream=fermentation_reactor.ins[spike_index])

    def load_threshold_conc_and_tau_max(self, threshold_conc, tau_max):
        fermentation_reactor = self.fermentation_reactor
        km = fermentation_reactor.nsk_kinetic_model
        cv = self.control_variables

        self.tau_max = tau_max
        fermentation_reactor.tau_max = tau_max

        km.set_value(cv.threshold_conc_var, threshold_conc)

        self._simulate_upstream_units()
        fermentation_reactor.simulate()

        d = fermentation_reactor.nsk_results_specific_tau_dict
        final_env_vol = d[cv.resolve_volume_col(fermentation_reactor)]
        vol_spike_added = d[cv.resolve_feed_volume_added_col(fermentation_reactor)]

        initial_env_vol = final_env_vol - vol_spike_added
        species_in_initial_feed = initial_env_vol * self.target_conc
        species_in_spikes = vol_spike_added * self.spike_conc
        self.splitter.split = species_in_initial_feed/(species_in_initial_feed+species_in_spikes) # split to initial feed

    def _solve_actuator(self, actuator, obj_f, desired_conc):
        try:
            IQ_interpolation(obj_f, actuator.lb, actuator.ub, ytol=1e-3)
        except Exception as e:
            raise FeedingStrategyError(
                f'Could not reach desired concentration {desired_conc} by '
                f'adjusting {actuator!r} within its bounds.') from e

    def _estimate_set_actuators(self,
                               desired_conc,
                               concentrator, diluter,
                               f_simulate_units_sequential,
                               terminal_stream):

        get_conc = self.get_conc
        f_simulate_units_sequential()

        if diluter.unit.outs[0].F_vol:
            def _concentrator_obj_f(value):
                concentrator.set(value)
                f_simulate_units_sequential()
                return get_conc(terminal_stream) - desired_conc

            def _diluter_obj_f(value):
                diluter.set(value)
                f_simulate_units_sequential()
                return get_conc(terminal_stream) - desired_conc

            _concentrator_obj_f(concentrator.lb)

            if _diluter_obj_f(diluter.lb) < 0: # if there is too low a conc even with no dilution
                self._solve_actuator(concentrator, _concentrator_obj_f, desired_conc)
            elif _diluter_obj_f(diluter.ub) > 0:
                self._solve_actuator(diluter, _diluter_obj_f, desired_conc)
            else:
                _concentrator_obj_f(concentrator.lb)
                self._solve_actuator(diluter, _diluter_obj_f, desired_conc)

    def _simulate_feed_units(self):
        self.splitter.simulate()
        for i in self.feed_units_sequential:
            i.simulate()

    def _simulate_spike_units(self):
        self.splitter.simulate()
        for i in self.spike_units_sequential:
            i.simulate()

    def _simulate_upstream_units(self):
        self._simulate_feed_units()
        self._simulate_spike_units()

    def get_conc(self, stream):
        return stream.imass[self.species_IDs].sum()/stream.ivol[self.solvent_ID]

    def get_feed_conc(self):
        return self.get_conc(stream=self.feed_diluter.unit.outs[0])

    def get_spike_conc(self):
        return self.get_conc(stream=self.spike_diluter.unit.outs[0])

    @property
    def current_specifications(self):
        return {k: self.__getattribute__(k) for k in self._spec_names}
    @property
    def baseline_specifications(self):
        return self._baseline_specifications
