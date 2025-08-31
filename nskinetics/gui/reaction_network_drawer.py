# -*- coding: utf-8 -*-
"""
Created on Sat Aug 30 22:35:15 2025

@author: sarangbhagwat
"""

import tkinter as tk
from tkinter import simpledialog, messagebox
from ..reactions import Reaction
from ..species import SpeciesSystem

__all__ = ('ReactionNetworkDrawerGUI',)

# nskinetics/gui/reaction_network_drawer.py
# -*- coding: utf-8 -*-
"""
ReactionNetworkDrawerGUI

Tkinter canvas to draw reaction networks for NSKinetics.

Modes:
- Select/Move (default): left-click empty = add species; left-drag node = move.
- Arrow (A→B): click source node, then destination node to create an arrow. (Also supports right/middle drag.)
- Build Reaction (multi): select reactants → Next → select products → Finish → stoich & kinetics.

UI:
- Toolbar with radio buttons to switch modes, plus Build Reaction flow buttons.
- Kinetic labels on all arrows (simple & hyper).
- Stoichiometry labels near nodes for multi-substrate/product reactions (hidden when coeff=1).

Saving:
- Converts everything to ReactionSystem.add_reaction(...) strings.
"""

import tkinter as tk
from tkinter import simpledialog, messagebox

class ReactionNetworkDrawerGUI:
    NODE_R = 22  # circle radius
    
    def __init__(self, root, system, *, node_radius: int = None):
        self.root = root
        self.system = system
        self.r = int(node_radius) if node_radius else self.NODE_R
        
        # Canvas
        self.canvas = tk.Canvas(root, bg="white", width=1000, height=700, highlightthickness=0)
        self.canvas.pack(fill="both", expand=True)
        
        # ===== Data =====
        # nodes: id -> {"species", "x","y", "circle","label"}
        self.nodes = {}
        # edges (simple): id -> {"src","dst","line","kin":{"kf","kb"},"kin_label"}
        self.edges = {}
        # hyper rxns: id -> {
        #   "reactants":[nids], "products":[nids],
        #   "r_stoich":{nid:coeff}, "p_stoich":{nid:coeff},
        #   "kin":{"kf","kb"}, "main","connectors":[ids],
        #   "kin_label", "r_txt":{nid:tid}, "p_txt":{nid:tid}
        # }
        self.hyper = {}
        
        # interaction state
        self.mode = tk.StringVar(value="select")  # "select" | "arrow" | "build"
        self.dragging_node = None
        
        # arrow (sequential) mode
        self.arrow_src_node = None
        
        # right/middle drag convenience
        self.drag_temp_line = None
        self.drag_start_node = None
        
        # build reaction (multi) state
        self.build_phase = None  # None | "reactants" | "products"
        self.build_reactants = set()
        self.build_products  = set()
        
        # ===== Bindings =====
        self.canvas.bind("<Button-1>", self._on_left_down)
        self.canvas.bind("<B1-Motion>", self._on_left_drag)
        self.canvas.bind("<ButtonRelease-1>", self._on_left_up)
        
        # Convenience: right/middle drag arrow (works in any mode except while picking build)
        for b in (2, 3):
            self.canvas.bind(f"<Button-{b}>", self._on_alt_down)
            self.canvas.bind(f"<B{b}-Motion>", self._on_alt_drag)
            self.canvas.bind(f"<ButtonRelease-{b}>", self._on_alt_up)
            
        # Double-click arrows to edit kinetics; eat single-clicks on arrows
        self.canvas.tag_bind("arrow", "<Double-Button-1>", self._on_arrow_double)
        self.canvas.tag_bind("arrow", "<Button-1>", lambda e: "break")
        
        # ===== Toolbar =====
        toolbar = tk.Frame(root)
        toolbar.pack(side="bottom", fill="x")
        
        # Left side: mode
        left = tk.Frame(toolbar)
        left.pack(side="left", padx=6, pady=6)
        
        tk.Label(left, text="Mode: ").pack(side="left")
        tk.Radiobutton(left, text="Select/Move", variable=self.mode, value="select",
                       command=self._sync_buttons).pack(side="left", padx=2)
        tk.Radiobutton(left, text="Arrow (A→B)", variable=self.mode, value="arrow",
                       command=self._sync_buttons).pack(side="left", padx=2)
        tk.Button(left, text="Build Reaction", command=self._start_build).pack(side="left", padx=8)
        
        # Build flow
        self.next_btn   = tk.Button(left, text="Next: Products", command=self._build_next)
        self.finish_btn = tk.Button(left, text="Finish Reaction", command=self._build_finish)
        self.cancel_btn = tk.Button(left, text="Cancel", command=self._build_cancel)
        self.next_btn.pack(side="left", padx=2)
        self.finish_btn.pack(side="left", padx=2)
        self.cancel_btn.pack(side="left", padx=2)
        
        # Right side: help/save
        right = tk.Frame(toolbar)
        right.pack(side="right", padx=6, pady=6)
        tk.Button(right, text="Help", command=self._show_help).pack(side="right", padx=6)
        tk.Button(right, text="Save Network", command=self._save).pack(side="right", padx=6)
        
        # Status
        self.status = tk.Label(toolbar, text="Ready", anchor="w")
        self.status.pack(side="left", padx=10)
        
        self._sync_buttons()

    # ===== Utilities =====

    def _hit_item(self, x, y):
        items = self.canvas.find_overlapping(x, y, x, y)
        return items[-1] if items else None

    def _node_of_item(self, item):
        for nid, n in self.nodes.items():
            if item in (n["circle"], n["label"]):
                return nid
        return None

    @staticmethod
    def _midpoint(p, q):
        return ((p[0]+q[0])/2.0, (p[1]+q[1])/2.0)

    @staticmethod
    def _towards(p_from, p_to, offset):
        import math
        dx, dy = p_to[0]-p_from[0], p_to[1]-p_from[1]
        d = math.hypot(dx, dy)
        if d == 0:
            return (p_from[0]+offset, p_from[1])
        ux, uy = dx/d, dy/d
        return (p_from[0]+ux*offset, p_from[1]+uy*offset)

    def _centroid(self, node_ids):
        if not node_ids:
            return None
        xs = [self.nodes[i]["x"] for i in node_ids]
        ys = [self.nodes[i]["y"] for i in node_ids]
        return (sum(xs)/len(xs), sum(ys)/len(ys))

    @staticmethod
    def _fmt_kin(kf, kb):
        return f"kf={kf:g}" if kb is None else f"kf={kf:g}, kb={kb:g}"

    @staticmethod
    def _fmt_coeff(c):
        return None if abs(c-1.0) < 1e-12 else f"{c:g}×"

    # ===== Button state / mode =====

    def _sync_buttons(self):
        # Arrow mode cancels any sequential arrow source and build flow
        if self.mode.get() != "arrow":
            self.arrow_src_node = None
        if self.mode.get() != "select":
            self.dragging_node = None
        if self.mode.get() != "build":
            self._build_reset_visuals(clear_phase=True)
            
        # Build flow buttons
        if self.mode.get() == "build":
            if self.build_phase == "reactants":
                self.next_btn.config(state=tk.NORMAL if self.build_reactants else tk.DISABLED)
                self.finish_btn.config(state=tk.DISABLED)
                self.cancel_btn.config(state=tk.NORMAL)
            elif self.build_phase == "products":
                self.next_btn.config(state=tk.DISABLED)
                self.finish_btn.config(state=tk.NORMAL if self.build_products else tk.DISABLED)
                self.cancel_btn.config(state=tk.NORMAL)
            else:
                self.next_btn.config(state=tk.DISABLED)
                self.finish_btn.config(state=tk.DISABLED)
                self.cancel_btn.config(state=tk.DISABLED)
        else:
            self.next_btn.config(state=tk.DISABLED)
            self.finish_btn.config(state=tk.DISABLED)
            self.cancel_btn.config(state=tk.DISABLED)
            
        self.status.config(text=f"Mode: {self.mode.get().capitalize()}")

    # ===== Node creation / dragging =====

    def _on_left_down(self, e):
        item = self._hit_item(e.x, e.y)
        
        # Click on arrow? consume so we don't create nodes accidentally
        if item and "arrow" in self.canvas.gettags(item):
            return "break"
        
        if self.mode.get() == "select":
            # click on node → start drag; else → create species
            nid = self._node_of_item(item) if item else None
            if nid is not None:
                self.dragging_node = nid
                return "break"
            # empty → create new species
            sp = simpledialog.askstring("New Species", "Enter species ID:")
            if not sp:
                return "break"
            self._create_node(sp.strip(), e.x, e.y)
            return "break"
        
        elif self.mode.get() == "arrow":
            nid = self._node_of_item(item) if item else None
            if nid is None:
                # click empty space: ignore
                return "break"
            if self.arrow_src_node is None:
                self.arrow_src_node = nid
                self.status.config(text=f"Arrow: source = {self.nodes[nid]['species']}. Select destination.")
            else:
                if nid != self.arrow_src_node:
                    self._create_simple_edge(self.arrow_src_node, nid, ask_kinetics=True)
                    self.status.config(text="Arrow created.")
                self.arrow_src_node = None
            return "break"

        elif self.mode.get() == "build":
            # toggle membership in current phase
            nid = self._node_of_item(item) if item else None
            if nid is None:
                return "break"
            if self.build_phase == "reactants":
                self._toggle_reactant(nid)
            elif self.build_phase == "products":
                self._toggle_product(nid)
            self._sync_buttons()
            return "break"
        
        return "break"

    def _on_left_drag(self, e):
        if self.mode.get() != "select":
            return
        if self.dragging_node is None:
            return
        self._move_node(self.dragging_node, e.x, e.y)

    def _on_left_up(self, e):
        self.dragging_node = None

    # ===== Alt/right drag arrows (optional convenience) =====

    def _on_alt_down(self, e):
        if self.mode.get() == "build":
            return
        item = self._hit_item(e.x, e.y)
        nid = self._node_of_item(item) if item else None
        if nid is None:
            return
        self.drag_start_node = nid
        p = (self.nodes[nid]["x"], self.nodes[nid]["y"])
        self.drag_temp_line = self.canvas.create_line(p[0], p[1], e.x, e.y,
                                                      arrow="last", dash=(4, 2),
                                                      width=1.6, fill="#333", tags=("arrow", "temp"))

    def _on_alt_drag(self, e):
        if self.drag_temp_line is None or self.drag_start_node is None:
            return
        p = (self.nodes[self.drag_start_node]["x"], self.nodes[self.drag_start_node]["y"])
        self.canvas.coords(self.drag_temp_line, p[0], p[1], e.x, e.y)

    def _on_alt_up(self, e):
        if self.drag_temp_line is None or self.drag_start_node is None:
            return
        item = self._hit_item(e.x, e.y)
        dst = self._node_of_item(item) if item else None
        self.canvas.delete(self.drag_temp_line)
        self.drag_temp_line = None
        if dst is not None and dst != self.drag_start_node:
            self._create_simple_edge(self.drag_start_node, dst, ask_kinetics=True)
        self.drag_start_node = None

    # ===== Arrow double-click (edit kinetics) =====

    def _on_arrow_double(self, e):
        arrow = self.canvas.find_withtag("current")
        if not arrow:
            return "break"
        aid = arrow[0]

        # simple edge?
        for ed in self.edges.values():
            if ed["line"] == aid:
                if self._edit_kinetics(ed["kin"]):
                    self._update_simple_edge(ed)
                return "break"
            
        # hyper main?
        for rxn in self.hyper.values():
            if rxn["main"] == aid:
                if self._edit_kinetics(rxn["kin"]):
                    self._update_hyper_geometry(rxn)
                return "break"
            
        return "break"

    # ===== Creation helpers =====

    def _create_node(self, species, x, y):
        r = self.r
        circle = self.canvas.create_oval(x-r, y-r, x+r, y+r,
                                         fill="#cae8ff", outline="#195a8a", width=1.6,
                                         tags=("node",))
        label  = self.canvas.create_text(x, y, text=species, tags=("label",))
        nid = len(self.nodes)
        self.nodes[nid] = {"species": species, "x": x, "y": y, "circle": circle, "label": label}
        return nid

    def _move_node(self, nid, x, y):
        n = self.nodes[nid]
        r = self.r
        self.canvas.coords(n["circle"], x-r, y-r, x+r, y+r)
        self.canvas.coords(n["label"],  x, y)
        n["x"], n["y"] = x, y
        
        # update simple edges
        for ed in self.edges.values():
            self._update_simple_edge(ed)
            
        # update hyper
        for rxn in self.hyper.values():
            self._update_hyper_geometry(rxn)

    def _create_simple_edge(self, src, dst, *, ask_kinetics=False):
        p = (self.nodes[src]["x"], self.nodes[src]["y"])
        q = (self.nodes[dst]["x"], self.nodes[dst]["y"])
        line = self.canvas.create_line(p[0], p[1], q[0], q[1],
                                       arrow="last", width=1.8, fill="#333",
                                       tags=("arrow", "simple"))
        kin = {"kf": 1.0, "kb": None}
        if ask_kinetics:
            if not self._edit_kinetics(kin):
                # user canceled → delete line and abort
                self.canvas.delete(line)
                return None
        mid = self._midpoint(p, q)
        kin_label = self.canvas.create_text(mid[0], mid[1]-14,
                                            text=self._fmt_kin(kin["kf"], kin["kb"]),
                                            fill="#444", font=("TkDefaultFont", 9),
                                            tags=("kin_label",))
        eid = len(self.edges)
        self.edges[eid] = {"src": src, "dst": dst, "line": line, "kin": kin, "kin_label": kin_label}
        return eid

    # ===== Simple edge updates =====

    def _update_simple_edge(self, ed):
        p = (self.nodes[ed["src"]]["x"], self.nodes[ed["src"]]["y"])
        q = (self.nodes[ed["dst"]]["x"], self.nodes[ed["dst"]]["y"])
        self.canvas.coords(ed["line"], p[0], p[1], q[0], q[1])
        mid = self._midpoint(p, q)
        self.canvas.coords(ed["kin_label"], mid[0], mid[1]-14)
        self.canvas.itemconfig(ed["kin_label"], text=self._fmt_kin(ed["kin"]["kf"], ed["kin"]["kb"]))

    # ===== Build reaction flow =====

    def _start_build(self):
        self.mode.set("build")
        self.build_phase = "reactants"
        self.build_reactants.clear()
        self.build_products.clear()
        self.status.config(text="Build: select REACTANTS (click nodes), then Next.")
        self._sync_buttons()

    def _toggle_reactant(self, nid):
        if nid in self.build_reactants:
            self.build_reactants.remove(nid)
            self._highlight(nid, False)
        else:
            self.build_reactants.add(nid)
            self._highlight(nid, True, color="#0b7d45")
        self.status.config(text=f"Reactants: {len(self.build_reactants)} selected.")

    def _toggle_product(self, nid):
        if nid in self.build_products:
            self.build_products.remove(nid)
            self._highlight(nid, False)
        else:
            self.build_products.add(nid)
            self._highlight(nid, True, color="#6e36bf")
        self.status.config(text=f"Products: {len(self.build_products)} selected.")

    def _build_next(self):
        if self.mode.get() != "build" or self.build_phase != "reactants":
            return
        if not self.build_reactants:
            messagebox.showwarning("Select reactants", "Pick at least one reactant.")
            return
        self.build_phase = "products"
        # remove green outlines; let products get purple
        for nid in list(self.build_reactants):
            self._highlight(nid, False)
        self.status.config(text="Build: select PRODUCTS (click nodes), then Finish.")
        self._sync_buttons()

    def _build_finish(self):
        if self.mode.get() != "build" or self.build_phase != "products":
            return
        if not self.build_products:
            messagebox.showwarning("Select products", "Pick at least one product.")
            return
        if self.build_reactants & self.build_products:
            messagebox.showerror("Invalid selection", "A species cannot be both reactant and product.")
            return
        
        reactants = sorted(self.build_reactants)
        products  = sorted(self.build_products)
        
        r_sto = self._prompt_stoich(reactants, role="reactant")
        if r_sto is None:
            return
        p_sto = self._prompt_stoich(products, role="product")
        if p_sto is None:
            return
        
        kin = {"kf": 1.0, "kb": None}
        if not self._edit_kinetics(kin):
            return
        
        self._create_hyper_reaction(reactants, products, r_sto, p_sto, kin)
        # clear visuals / exit
        for nid in reactants + products:
            self._highlight(nid, False)
        self._build_reset_visuals(clear_phase=True)
        self.mode.set("select")
        self.status.config(text="Reaction added.")
        self._sync_buttons()

    def _build_cancel(self):
        for nid in list(self.build_reactants) + list(self.build_products):
            self._highlight(nid, False)
        self._build_reset_visuals(clear_phase=True)
        self.mode.set("select")
        self.status.config(text="Build canceled.")
        self._sync_buttons()

    def _build_reset_visuals(self, *, clear_phase=False):
        self.build_reactants.clear()
        self.build_products.clear()
        if clear_phase:
            self.build_phase = None

    def _prompt_stoich(self, nids, role="reactant"):
        out = {}
        for nid in nids:
            sp = self.nodes[nid]["species"]
            v = simpledialog.askstring("Stoichiometry", f"{role.capitalize()} '{sp}' coefficient (default 1):", initialvalue="1")
            if v is None:
                return None
            v = v.strip()
            if v == "":
                c = 1.0
            else:
                try:
                    c = float(v)
                except ValueError:
                    messagebox.showerror("Invalid input", f"Coefficient for {sp} must be a number.")
                    return None
            if c <= 0:
                messagebox.showerror("Invalid input", f"Coefficient for {sp} must be > 0.")
                return None
            out[nid] = c
        return out

    # ===== Hyper reaction (multi) =====

    def _create_hyper_reaction(self, reactants, products, r_sto, p_sto, kin):
        rc = self._centroid(reactants)
        pc = self._centroid(products)
        if rc is None or pc is None:
            return
        
        main = self.canvas.create_line(rc[0], rc[1], pc[0], pc[1],
                                       arrow="last", width=2.2, fill="#222",
                                       tags=("arrow", "hyper"))
        connectors = []
        r_txt = {}
        p_txt = {}
        
        # connectors + stoich labels
        for nid in reactants:
            n = self.nodes[nid]
            cid = self.canvas.create_line(n["x"], n["y"], rc[0], rc[1], dash=(3, 2), fill="#888", tags=("connector",))
            connectors.append(cid)
            text = self._fmt_coeff(r_sto.get(nid, 1.0))
            if text:
                pos = self._towards((n["x"], n["y"]), rc, self.r + 10)
                tid = self.canvas.create_text(pos[0], pos[1], text=text, fill="#555",
                                              font=("TkDefaultFont", 9), tags=("stoich_label",))
                r_txt[nid] = tid
                
        for nid in products:
            n = self.nodes[nid]
            cid = self.canvas.create_line(pc[0], pc[1], n["x"], n["y"], dash=(3, 2), fill="#888", tags=("connector",))
            connectors.append(cid)
            text = self._fmt_coeff(p_sto.get(nid, 1.0))
            if text:
                pos = self._towards((n["x"], n["y"]), pc, self.r + 10)
                tid = self.canvas.create_text(pos[0], pos[1], text=text, fill="#555",
                                              font=("TkDefaultFont", 9), tags=("stoich_label",))
                p_txt[nid] = tid
                
        mid = self._midpoint(rc, pc)
        kin_label = self.canvas.create_text(mid[0], mid[1]-14, text=self._fmt_kin(kin["kf"], kin["kb"]),
                                            fill="#444", font=("TkDefaultFont", 9), tags=("kin_label",))
        
        hid = len(self.hyper)
        self.hyper[hid] = {
            "reactants": reactants, "products": products,
            "r_stoich": r_sto, "p_stoich": p_sto,
            "kin": kin, "main": main, "connectors": connectors,
            "kin_label": kin_label, "r_txt": r_txt, "p_txt": p_txt,
        }

    def _update_hyper_geometry(self, rxn):
        rc = self._centroid(rxn["reactants"])
        pc = self._centroid(rxn["products"])
        if rc is None or pc is None:
            return
        
        self.canvas.coords(rxn["main"], rc[0], rc[1], pc[0], pc[1])
        
        # remove & redraw connectors
        for cid in rxn["connectors"]:
            try:
                self.canvas.delete(cid)
            except Exception:
                pass
        rxn["connectors"].clear()

        # reactants
        for nid in rxn["reactants"]:
            n = self.nodes[nid]
            cid = self.canvas.create_line(n["x"], n["y"], rc[0], rc[1], dash=(3, 2), fill="#888", tags=("connector",))
            rxn["connectors"].append(cid)
            # stoich label
            coeff = rxn["r_stoich"].get(nid, 1.0)
            text = self._fmt_coeff(coeff)
            if text:
                pos = self._towards((n["x"], n["y"]), rc, self.r + 10)
                tid = rxn["r_txt"].get(nid)
                if tid is None:
                    tid = self.canvas.create_text(pos[0], pos[1], text=text, fill="#555",
                                                  font=("TkDefaultFont", 9), tags=("stoich_label",))
                    rxn["r_txt"][nid] = tid
                else:
                    self.canvas.coords(tid, pos[0], pos[1])
                    self.canvas.itemconfig(tid, text=text)
            else:
                tid = rxn["r_txt"].pop(nid, None)
                if tid:
                    try:
                        self.canvas.delete(tid)
                    except Exception:
                        pass

        # products
        for nid in rxn["products"]:
            n = self.nodes[nid]
            cid = self.canvas.create_line(pc[0], pc[1], n["x"], n["y"], dash=(3, 2), fill="#888", tags=("connector",))
            rxn["connectors"].append(cid)
            coeff = rxn["p_stoich"].get(nid, 1.0)
            text = self._fmt_coeff(coeff)
            if text:
                pos = self._towards((n["x"], n["y"]), pc, self.r + 10)
                tid = rxn["p_txt"].get(nid)
                if tid is None:
                    tid = self.canvas.create_text(pos[0], pos[1], text=text, fill="#555",
                                                  font=("TkDefaultFont", 9), tags=("stoich_label",))
                    rxn["p_txt"][nid] = tid
                else:
                    self.canvas.coords(tid, pos[0], pos[1])
                    self.canvas.itemconfig(tid, text=text)
            else:
                tid = rxn["p_txt"].pop(nid, None)
                if tid:
                    try:
                        self.canvas.delete(tid)
                    except Exception:
                        pass

        # kinetic label
        mid = self._midpoint(rc, pc)
        self.canvas.coords(rxn["kin_label"], mid[0], mid[1]-14)
        self.canvas.itemconfig(rxn["kin_label"], text=self._fmt_kin(rxn["kin"]["kf"], rxn["kin"]["kb"]))

    # ===== Edit kinetics dialog (returns True if confirmed) =====

    def _edit_kinetics(self, kin_dict):
        kf = simpledialog.askfloat("Kinetics", "Enter kf:", initialvalue=kin_dict["kf"])
        if kf is None:
            return False
        kb_str = simpledialog.askstring("Kinetics", "Enter kb (leave blank for irreversible):",
                                        initialvalue="" if kin_dict["kb"] is None else str(kin_dict["kb"]))
        if kb_str is None:
            return False
        kb = None if kb_str.strip() == "" else self._safe_float(kb_str)
        if kb_str.strip() != "" and kb is None:
            messagebox.showerror("Error", "kb must be a number or blank.")
            return False
        kin_dict["kf"], kin_dict["kb"] = kf, kb
        return True

    @staticmethod
    def _safe_float(s):
        try:
            return float(s)
        except Exception:
            return None

    # ===== Visual highlight =====

    def _highlight(self, nid, on, color=None):
        circle = self.nodes[nid]["circle"]
        if on:
            self.canvas.itemconfig(circle, outline=color if color else "#d45500", width=3)
        else:
            self.canvas.itemconfig(circle, outline="#195a8a", width=1.6)

    # ===== Save =====

    @staticmethod
    def _fmt_species(sp, coeff):
        if coeff == 1 or abs(coeff-1.0) < 1e-12:
            return sp
        return f"{coeff:g} {sp}"

    def _save(self):
        # species
        for k, v in self.nodes.items():
            self.system.species_system.add_species(v['species'])
        
        # simple edges
        for ed in self.edges.values():
            a = self.nodes[ed["src"]]["species"]
            b = self.nodes[ed["dst"]]["species"]
            kf = ed["kin"]["kf"]; kb = ed["kin"]["kb"]
            if kb is None:
                eq = f"{a} -> {b}; kf = {kf}"
            else:
                eq = f"{a} <-> {b}; kf = {kf}, kb = {kb}"
            self.system.add_reaction(eq)

        # hyper reactions
        for rxn in self.hyper.values():
            r_terms = [self._fmt_species(self.nodes[n]["species"], rxn["r_stoich"].get(n, 1.0))
                       for n in rxn["reactants"]]
            p_terms = [self._fmt_species(self.nodes[n]["species"], rxn["p_stoich"].get(n, 1.0))
                       for n in rxn["products"]]
            lhs = " + ".join(r_terms); rhs = " + ".join(p_terms)
            kf = rxn["kin"]["kf"]; kb = rxn["kin"]["kb"]
            if kb is None:
                eq = f"{lhs} -> {rhs}; kf = {kf}"
            else:
                eq = f"{lhs} <-> {rhs}; kf = {kf}, kb = {kb}"
            self.system.add_reaction(eq)

        messagebox.showinfo("Saved", "Reaction network saved to ReactionSystem.")
        self.root.destroy()

    # ===== Help =====

    def _show_help(self):
        msg = (
            "Usage\n\n"
            "Select/Move (default):\n"
            "  • Left-click empty space: add a species\n"
            "  • Left-drag a node: move it (arrows & labels follow)\n\n"
            "Arrow (A→B):\n"
            "  • Click a source node, then a destination node to create an arrow\n"
            "  • You will be prompted for kf/kb; a label appears on the arrow\n"
            "  • (Optional) Right/Middle-drag from a node to another to create quickly\n\n"
            "Build Reaction (multi-substrate/product):\n"
            "  • Click 'Build Reaction' to enter build mode\n"
            "  • Phase 1: click nodes to select REACTANTS, then 'Next: Products'\n"
            "  • Phase 2: click nodes to select PRODUCTS, then 'Finish Reaction'\n"
            "  • You’ll set stoichiometries for each selected species, then kinetics\n"
            "  • A grouped arrow (centroid→centroid) is drawn, with stoich labels near species\n\n"
            "Double-click any arrow text/path to edit kinetics. Click 'Save Network' to write reactions."
        )
        messagebox.showinfo("Help", msg)


