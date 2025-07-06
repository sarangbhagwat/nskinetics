import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import time

__all__ = ('ReactionSystemGUI',)

# GUI
class ReactionSystemGUI:
    def __init__(self, root, system,
                 initial_conc_scrollparams=(0, 5.0, 0.01),
                 rxn_kinetic_param_scrollparams=(0, 500, 0.01),
                 tspan_ub_scrollparams=(10, 7*24*3600, 10),
                 log_transform_concs=True,
                 timeout_solve=0.2,
                 ):
        self.root = root
        self.system = system
        self.log_transform_concs = log_transform_concs
        system._timeout_solve_ivp = timeout_solve
        # self.initial_concentration_scrollparams = initial_concentration_scrollparams
        # self.rxn_kinetic_param_scrollparams = rxn_kinetic_param_scrollparams
        
        
        # top bar frame
        topbar = ttk.Frame(root)
        topbar.pack(side=tk.TOP, fill=tk.X)
        
        # place quit button on the right
        quit_btn = ttk.Button(topbar, text="Quit", command=root.destroy)
        quit_btn.pack(side=tk.RIGHT, padx=5, pady=5)

        ### Change Initial Concentrations, Time Span ###
        # frame
        init_frame = ttk.LabelFrame(root, 
                                    # text="Initial concentrations",
                                    )
        init_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        
        self.init_conc_vars = init_conc_vars = []
        sp_sys = system.species_system
        self.all_sp_IDs = all_sp_IDs = sp_sys.all_sp_IDs
        
        #!!!
        self.events = None
        self.spikes = None
        
        # Time span upper bound
        
        self.t_span = t_span = (0., 300.)
        self.t_span_ub_var = t_span_ub_var = tk.DoubleVar(value=t_span[1])
        # time span ub label
        row = ttk.Frame(init_frame)
        row.pack(fill=tk.X, pady=2)
        label = ttk.Label(row, text='Duration')
        label.pack(side=tk.LEFT)
        # time span ub entry box
        entry = ttk.Entry(row, textvariable=t_span_ub_var, width=10)
        entry.pack(side=tk.RIGHT)
        # time span ub scroll bar
        from_, to_, res_ = tspan_ub_scrollparams
        scale = tk.Scale(
            init_frame, variable=t_span_ub_var,
            from_=from_, to=to_, resolution=res_,
            orient=tk.HORIZONTAL,
            length=150)
        scale.pack()
        # update time span ub
        t_span_ub_var.trace_add("write", self.on_t_span_ub_change)
        #
        
        # Initial concentrations
        # sp_sys.concentrations[1] = 2.
        self.initial_concentrations = sp_sys.concentrations.copy()
        
        for species, val in zip(all_sp_IDs, sp_sys.concentrations):
            
            var = tk.DoubleVar(value=val)
            
            # map tkinter vars
            init_conc_vars.append(var)

            # species label
            row = ttk.Frame(init_frame)
            row.pack(fill=tk.X, pady=2)
            label = ttk.Label(row, text=species)
            label.pack(side=tk.LEFT)

            # species initial concentration entry box
            entry = ttk.Entry(row, textvariable=var, width=10)
            entry.pack(side=tk.RIGHT)

            # species initial concentration entry scroll bar
            
            from_, to_, res_ = initial_conc_scrollparams
            scale = tk.Scale(
                init_frame, variable=var,
                from_=from_, to=to_, resolution=res_,
                orient=tk.HORIZONTAL,
                length=150
            )
            scale.pack()
            
            # update the initial concentration
            var.trace_add("write", self.on_init_conc_change)
        ###
        
        ### Change Kinetic Parameters ###
        # frame
        param_frame = ttk.LabelFrame(root, 
                                     text="Kinetic parameters",
                                     )
        param_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)

        param_vals = system._get_reaction_kinetic_params()
        param_keys = system._reaction_kinetic_param_keys
        
        # self.param_vars = param_vars =\
        #     {k:v for k,v in zip (param_keys, param_vals)}
        self.param_vars = param_vars = []
        
        prev_rxn = None
        
        for k, v in zip(param_keys, param_vals):
            r, param = k
            
            var = tk.DoubleVar(value=v)
            
            # map tkinter vars
            param_vars.append(var)
            
            # reaction label
            if not r==prev_rxn:
                row = ttk.Frame(param_frame)
                row.pack(fill=tk.X, pady=2)
                label = ttk.Label(row, text=r)
                label.pack(side=tk.LEFT)
                prev_rxn = r
                
            # parameter label
            row = ttk.Frame(param_frame)
            row.pack(fill=tk.X, pady=2)
            label = ttk.Label(row, text=param)
            label.pack(side=tk.LEFT)
            
            # parameter entry box
            entry = ttk.Entry(row, textvariable=var, width=10)
            entry.pack(side=tk.RIGHT)
            
            # parameter scroll bar
            from_, to_, res_ = rxn_kinetic_param_scrollparams
            scale = tk.Scale(
                param_frame, variable=var,
                from_=from_, to=to_, resolution=res_,
                orient=tk.HORIZONTAL,
                length=150
            )
            scale.pack()
            
            # update the param value
            var.trace_add("write", self.on_param_change)
        ###
        
        ### Change Species Shown in Plot ###
        # frame
        init_frame = ttk.LabelFrame(root, text="Plot Species")
        init_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)

        self.sps_to_include_vars = sps_to_include_vars = []
        self.sps_to_include = sps_to_include = all_sp_IDs.copy()
        # sp_sys = system.species_system
        # all_sp_IDs = sp_sys.all_sp_IDs
        
        for species in all_sp_IDs:
            
            var = tk.BooleanVar(value=True)
            
            # map tkinter vars
            sps_to_include_vars.append(var)

            # species initial concentration entry scroll bar
            check = ttk.Checkbutton(param_frame, 
                                    text=species, 
                                    variable=var)
            check.pack(pady=5)
            
            # update the initial concentration
            var.trace_add("write", self.on_sps_included_change)
        ###
        
        
        ### Change Events ###
        
        ###
        
        
        ### Change Spikes ###
        
        ###
        
        
        ### Build Matplotlib Figure
        fig, ax = plt.subplots(figsize=(5,4))
        
        self.system.solve(
                          t_span=self.t_span,
                          events=self.events,
                          spikes=self.spikes,
                          log_transform_concs=self.log_transform_concs,
                          save_events_df=False,
                          )
        
        system.plot_solution(fig=fig, ax=ax, sps_to_include=sps_to_include,
                             auto_ticks=False, show=False)
        self.fig = fig
        self.ax = ax
        
        self.canvas = FigureCanvasTkAgg(fig, master=root)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        # self.simulate_and_update_plot()
        ###
    
    def on_t_span_ub_change(self, *_):
        # try:
        self.t_span = (0., self.t_span_ub_var.get())
        # except:
        #     pass
        self.simulate_and_update_plot()

    def on_init_conc_change(self, *_):
        # try:
        self.initial_concentrations =\
            np.array([v.get() for v in self.init_conc_vars])
        # except:
        #     pass
        self.simulate_and_update_plot()
            
    def on_param_change(self, *_):
        # try:
        self.system.set_reaction_kinetic_params(np.array([v.get() 
                                                 for v in self.param_vars]))
        # except:
        #     pass
        self.simulate_and_update_plot()
    
    def on_sps_included_change(self, *_):
        # try:
        all_sp_IDs = self.all_sp_IDs
        self.sps_to_include = sps_to_include = []
        for v, sp in zip(self.sps_to_include_vars, self.all_sp_IDs):
            if v.get():
                sps_to_include.append(sp)
        # except:
        #     pass
        self.update_only_plot()
    
    def load_initial_concentrations(self):
        self.system.species_system.concentrations = self.initial_concentrations
        
    def simulate(self):
        self.load_initial_concentrations()
        self.system.solve(
                          t_span=self.t_span,
                          events=self.events,
                          spikes=self.spikes,
                          log_transform_concs=self.log_transform_concs,
                          save_events_df=False
                          )
    
    def update_only_plot(self):
        self.ax.clear()
        self.system.plot_solution(fig=self.fig, ax=self.ax, 
                                  sps_to_include=self.sps_to_include,
                                  auto_ticks=False, show=False)
        self.canvas.draw()
        
    def simulate_and_update_plot(self):
        # t, y = self.system.simulate()
        # self.ax.clear()
        # self.ax.plot(t, y)
        # self.ax.set_xlabel("Time")
        # self.ax.set_ylabel("Concentration")
        self.simulate()
        self.update_only_plot()
        # print(self.system.__str__())
        # print(self.system.species_system.concentrations)
        
# if __name__ == "__main__":
#     root = tk.Tk()
#     root.title("NSKinetics - Reaction System GUI")
#     app = ReactionSystemGUI(root, system)
#     root.mainloop()
