# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

__all__ = ('ReactionSystemGUI',)

# GUI

# scrollable frame helper

# scrollable frame helper
class ScrollableFrame(ttk.Frame):
    def __init__(self, container, width=400, height=700, *args, **kwargs):
        super().__init__(container, *args, **kwargs)
        style = ttk.Style()
        bg_color = style.lookup("TFrame", "background") or "#f0f0f0"
        # bg_color = "#3a3a3a"
        
        self.canvas = tk.Canvas(self, borderwidth=0, width=width, height=height,
                                background=bg_color, highlightthickness=0)
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        
        self.scrollable_frame = ttk.Frame(self.canvas)
        self.scrollable_frame.configure(style="TFrame")
        
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(
                scrollregion=self.canvas.bbox("all")
            )
        )
        
        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=scrollbar.set)
            
        self.canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
            
        # enable mousewheel scroll
        self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)
        
    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(int(-1*(event.delta/120)), "units")


class ReactionSystemGUI:
    def __init__(self, root, system,
                 initial_conc_scrollparams=(0, 5.0, 0.05),
                 rxn_kinetic_param_scrollparams=(0, 500, 5),
                 tspan_ub_scrollparams=(10, 7*24*3600, 10),
                 timeout_solve=0.2):
        
        self.root = root
        self.system = system
        system._timeout_solve_ivp = timeout_solve
        
        style = ttk.Style()
        style.theme_use("clam")
        style.configure(".", font=("Arial", 10))
        
        root.title(f"NSKinetics: Reaction System {system.ID}")
        root.geometry("1300x750")
        
        # top bar
        topbar = ttk.Frame(root, padding=5)
        topbar.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(topbar, text="Reaction System Control Panel", font=("Arial", 12, "bold")).pack(side=tk.LEFT)
        ttk.Button(topbar, text="Quit", command=root.destroy).pack(side=tk.RIGHT)
        
        # left side controls with a single vertical scroll
        controls_outer = ttk.Frame(root, padding=5)
        controls_outer.pack(side=tk.LEFT, fill=tk.BOTH, padx=5, pady=5)
        
        controls_scroll = ScrollableFrame(controls_outer, width=500, height=700)
        controls_scroll.pack(fill=tk.BOTH, expand=True)
        controls_frame = controls_scroll.scrollable_frame
        
        # ========== Initial Concentrations & Time ==========
        init_frame = ttk.LabelFrame(controls_frame, text="Initial Conditions", padding=5)
        init_frame.pack(fill="x", anchor="n", expand=True, padx=5, pady=5)
        
        self.init_conc_vars = []
        sp_sys = system.species_system
        self.all_sp_IDs = sp_sys.all_sp_IDs
        
        self.t_span = (0., 300.)
        self.t_span_ub_var = tk.DoubleVar(value=self.t_span[1])
        
        ttk.Label(init_frame, text="Simulation Duration (s)").pack(anchor=tk.W, pady=2)
        entry = ttk.Entry(init_frame, textvariable=self.t_span_ub_var, width=10)
        entry.pack(anchor=tk.W, pady=2, fill="x")
        ttk.Scale(init_frame, variable=self.t_span_ub_var,
                  from_=tspan_ub_scrollparams[0],
                  to=tspan_ub_scrollparams[1],
                  orient=tk.HORIZONTAL,
                  length=300).pack(fill="x", pady=2)
        # entry.bind("<Return>", lambda e: self.on_t_span_ub_change())
        # entry.bind("<FocusOut>", lambda e: self.on_t_span_ub_change())
        self.t_span_ub_var.trace_add("write", self.on_t_span_ub_change)
        
        ttk.Separator(init_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=5)
        ttk.Label(init_frame, text="Initial Concentrations").pack(anchor=tk.W)
        
        self.initial_concentrations = sp_sys.concentrations.copy()
        
        for species, val in zip(self.all_sp_IDs, sp_sys.concentrations):
            var = tk.DoubleVar(value=val)
            self.init_conc_vars.append(var)
            
            subframe = ttk.Frame(init_frame)
            subframe.pack(fill="x", pady=2)
            
            ttk.Label(subframe, text=species).pack(side=tk.LEFT)
            entry = ttk.Entry(subframe, textvariable=var, width=16)
            entry.pack(side=tk.RIGHT, padx=2)
            ttk.Scale(init_frame, variable=var,
                      from_=initial_conc_scrollparams[0],
                      to=initial_conc_scrollparams[1],
                      orient=tk.HORIZONTAL,
                      length=300).pack(fill="x", pady=2)
            # entry.bind("<Return>", lambda e: self.on_init_conc_change())
            # entry.bind("<FocusOut>", lambda e: self.on_init_conc_change())
            var.trace_add("write", self.on_init_conc_change)
        
        # ========== Kinetic Parameters ==========
        param_frame = ttk.LabelFrame(controls_frame, text="Kinetic Parameters", padding=5)
        param_frame.pack(fill="x", anchor="n", expand=True, padx=5, pady=5)
        
        self.param_vars = []
        param_vals = system._get_reaction_kinetic_params()
        param_keys = system._reaction_kinetic_param_keys
        
        current_rxn = None
        
        for k, v in zip(param_keys, param_vals):
            rxn, param = k
            var = tk.DoubleVar(value=v)
            self.param_vars.append(var)
            
            if rxn != current_rxn:
                header = ttk.Label(param_frame, text=f"Reaction: {rxn}", font=("Arial", 10, "bold"))
                header.pack(anchor="w", pady=(5,2))
                current_rxn = rxn
                
            subframe = ttk.Frame(param_frame)
            subframe.pack(fill="x", pady=2)
            ttk.Label(subframe, text=param).pack(side=tk.LEFT)
            entry = ttk.Entry(subframe, textvariable=var, width=16)
            entry.pack(side=tk.RIGHT, padx=2)
            ttk.Scale(param_frame, variable=var,
                      from_=rxn_kinetic_param_scrollparams[0],
                      to=rxn_kinetic_param_scrollparams[1],
                      orient=tk.HORIZONTAL,
                      length=300).pack(fill="x", pady=2)
            # entry.bind("<Return>", lambda e: self.on_param_change())
            # entry.bind("<FocusOut>", lambda e: self.on_param_change())
            var.trace_add("write", self.on_param_change)
        
        # ========== Plot Species ==========
        species_frame = ttk.LabelFrame(controls_frame, text="Plot Species", padding=5)
        species_frame.pack(fill="x", anchor="n", expand=True, padx=5, pady=5)
        
        self.sps_to_include_vars = []
        self.sps_to_include = self.all_sp_IDs.copy()
        
        for species in self.all_sp_IDs:
            var = tk.BooleanVar(value=True)
            self.sps_to_include_vars.append(var)
            ttk.Checkbutton(species_frame, text=species, variable=var).pack(anchor="w", pady=2)
            var.trace_add("write", self.on_sps_included_change)
        
        # ========== Plot ==========
        fig, ax = plt.subplots(figsize=(6,5))
        
        self.system.solve(t_span=self.t_span, events=None, spikes=None, save_events_df=False)
        system.plot_solution(fig=fig, ax=ax, sps_to_include=self.sps_to_include,
                             auto_ticks=False, show=False)
        
        self.fig = fig
        self.ax = ax
        
        self.canvas = FigureCanvasTkAgg(fig, master=root)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
    # event handlers
    def on_t_span_ub_change(self, *_):
        self.t_span = (0., self.t_span_ub_var.get())
        self.simulate_and_update_plot()

    def on_init_conc_change(self, *_):
        self.initial_concentrations = np.array([v.get() for v in self.init_conc_vars])
        self.simulate_and_update_plot()
        
    def on_param_change(self, *_):
        self.system.set_reaction_kinetic_params(np.array([v.get() for v in self.param_vars]))
        self.simulate_and_update_plot()
    
    def on_sps_included_change(self, *_):
        self.sps_to_include = [sp for v, sp in zip(self.sps_to_include_vars, self.all_sp_IDs) if v.get()]
        self.update_only_plot()
    
    def load_initial_concentrations(self):
        self.system.species_system.concentrations = self.initial_concentrations
        
    def simulate(self):
        self.load_initial_concentrations()
        self.system.solve(t_span=self.t_span, events=None, spikes=None, save_events_df=False)
        
    def update_only_plot(self):
        self.ax.clear()
        self.system.plot_solution(fig=self.fig, ax=self.ax,
                                  sps_to_include=self.sps_to_include,
                                  auto_ticks=False, show=False)
        self.canvas.draw()
    
    def simulate_and_update_plot(self):
        self.root.focus()   # force any Entry to lose focus and validate
        self.simulate()
        self.update_only_plot()



