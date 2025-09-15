# -*- coding: utf-8 -*-
# NSKinetics: simulation of Non-Steady state enzyme Kinetics and inhibitory phenomena
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
# 
# This module is under the MIT open-source license. See 
# https://github.com/sarangbhagwat/nskinetics/blob/main/LICENSE
# for license details.

from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Sequence, Tuple, Dict, Any, List, Optional

import numpy as np
import pickle
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool

try:
    from scipy.optimize import differential_evolution, basinhopping, minimize
except Exception:
    differential_evolution = None
    basinhopping = None
    minimize = None

__all__ = ('hybrid_global_optimize',)

#%%
@dataclass
class GlobalOptResult:
    x: np.ndarray
    fun: float
    success: bool
    message: str
    nfev: int
    history: List[Dict[str, Any]]
    meta: Dict[str, Any]


def _ensure_bounds(bounds: Sequence[Tuple[float, float]]) -> np.ndarray:
    B = np.asarray(bounds, dtype=float)
    if B.ndim != 2 or B.shape[1] != 2:
        raise ValueError("bounds must be a sequence of (low, high) pairs.")
    if np.any(B[:, 0] >= B[:, 1]):
        raise ValueError("Each bound must satisfy low < high.")
    return B


def hybrid_global_optimize(
    objective: Callable[[np.ndarray], float],
    bounds: Sequence[Tuple[float, float]],
    *,
    # Parallelism control
    parallel: str = "threads",           # "none" | "threads" | "processes"
    workers: int | None | str = None,    # -1, int, None, or "auto"
    # Differential Evolution
    de_popsize: int = 15,
    de_maxiter: int = 2000,
    de_tol: float = 1e-6,
    de_strategy: str = "best1bin",
    de_mutation: tuple = (0.5, 1.0),
    de_recombination: float = 0.7,
    # Basin Hopping
    do_basinhop: bool = True,
    bh_niter: int = 10,
    bh_stepsize: float = 0.5,
    # Local polish
    local_method: str = "L-BFGS-B",
    local_maxiter: int = 1000,
    polish_top_k: int = 5,
    # General
    random_state: Optional[int] = None,
    callback: Optional[Callable[[np.ndarray, float, Dict[str, Any]], None]] = None,
    return_all: bool = False,
    **minimize_kwargs,
) -> GlobalOptResult:
    """
    Hybrid global optimization using DE → (optional) Basin-Hopping → local polish.

    `parallel`:
        "threads" (safe in notebooks/Windows; default), "processes" (requires __main__ guard),
        or "none" for serial.
    `workers`:
        None/1 (serial), int for #threads/processes, or -1/"auto" for all CPUs.
    """
    if differential_evolution is None or minimize is None:
        raise ImportError("SciPy (`scipy.optimize`) is required.")

    rng = np.random.default_rng(random_state)
    B = _ensure_bounds(bounds)
    dim = len(B)

    # ---- Normalize workers
    w = workers
    if w == "auto" or (isinstance(w, int) and w < 0):
        try:
            w = max(1, cpu_count())
        except Exception:
            w = None
    if w in (None, 1):
        parallel = "none"

    # ---- Detect picklability (relevant for processes)
    picklable = True
    try:
        pickle.dumps(objective)
    except Exception:
        picklable = False

    history: List[Dict[str, Any]] = []
    total_evals = 0

    # ---------- 1) Differential Evolution ----------
    pool = None
    de_kwargs = dict(
        func=objective,
        bounds=bounds,
        strategy=de_strategy,
        maxiter=de_maxiter,
        popsize=de_popsize,
        tol=de_tol,
        mutation=de_mutation,
        recombination=de_recombination,
        seed=random_state,
        updating="deferred" if parallel != "none" else "immediate",
    )

    if parallel == "threads" and w not in (None, 1):
        pool = ThreadPool(processes=w)
        de_kwargs["workers"] = pool.map
    elif parallel == "processes" and w not in (None, 1):
        if not picklable:
            # fall back safely
            parallel = "none"
        else:
            de_kwargs["workers"] = w

    try:
        de_res = differential_evolution(**de_kwargs)
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    total_evals += int(getattr(de_res, "nfev", 0))

    if return_all:
        history.append({
            "stage": "differential_evolution",
            "x": de_res.x.copy(),
            "fun": float(de_res.fun),
            "nfev": int(getattr(de_res, "nfev", 0)),
            "nit": int(getattr(de_res, "nit", 0)),
            "message": str(de_res.message),
        })
    if callback:
        callback(de_res.x, float(de_res.fun), {"stage": "differential_evolution", "result": de_res})

    best_x = de_res.x.copy()
    best_f = float(de_res.fun)

    # ---------- 2) Basin-Hopping (optional) ----------
    if do_basinhop and basinhopping is not None:
        minimizer_kwargs = dict(method=local_method, bounds=bounds, options={"maxiter": local_maxiter})
        minimizer_kwargs.update(minimize_kwargs)
        bh = basinhopping(objective, best_x, niter=bh_niter, stepsize=bh_stepsize,
                          minimizer_kwargs=minimizer_kwargs, seed=random_state, disp=False)
        bh_x = np.asarray(bh.x, dtype=float)
        bh_f = float(bh.fun)
        if bh_f < best_f:
            best_x, best_f = bh_x.copy(), bh_f
        if return_all:
            history.append({"stage": "basin_hopping", "x": bh_x.copy(), "fun": bh_f,
                            "niter": bh_niter, "message": "Basin-hopping finished"})
        if callback:
            callback(bh_x, bh_f, {"stage": "basin_hopping"})

    # ---------- 3) Local polish from top-k DE candidates ----------
    candidates: List[Tuple[float, np.ndarray]] = [(float(de_res.fun), de_res.x.copy())]
    try:
        pop = np.asarray(de_res.population)
        pop_f = np.asarray(de_res.population_energies)
        idx = np.argsort(pop_f)[:max(polish_top_k, 1)]
        for i in idx:
            candidates.append((float(pop_f[i]), pop[i].copy()))
    except Exception:
        for _ in range(max(polish_top_k - 1, 0)):
            r = rng.random(dim); x = B[:, 0] + r * (B[:, 1] - B[:, 0])
            f = float(objective(x)); total_evals += 1
            candidates.append((f, x))

    # Deduplicate and sort
    uniq: Dict[tuple, float] = {}
    for f, x in candidates:
        key = tuple(np.round(x, 12))
        if key not in uniq or f < uniq[key]:
            uniq[key] = f
    cand_sorted = sorted(((f, np.array(x)) for x, f in uniq.items()), key=lambda t: t[0])[:max(polish_top_k, 1)]

    local_best = (best_f, best_x.copy())
    local_records = []
    for rank, (f0, x0) in enumerate(cand_sorted, 1):
        res_local = minimize(objective, x0, method=local_method, bounds=bounds,
                             options={"maxiter": local_maxiter}, **minimize_kwargs)
        total_evals += int(getattr(res_local, "nfev", 0))
        local_records.append({
            "rank": rank,
            "x0": x0.copy(), "f0": float(f0),
            "x": res_local.x.copy(), "fun": float(res_local.fun),
            "success": bool(res_local.success),
            "message": str(res_local.message),
            "nfev": int(getattr(res_local, "nfev", 0)),
            "nit": int(getattr(res_local, "nit", 0)) if hasattr(res_local, "nit") else None,
        })
        if res_local.fun < local_best[0]:
            local_best = (float(res_local.fun), res_local.x.copy())
        if callback:
            callback(res_local.x, float(res_local.fun), {"stage": "local_polish", "rank": rank})

    if return_all:
        history.append({"stage": "local_polish", "records": local_records,
                        "message": f"Polished {len(local_records)} starts with {local_method}"})

    best_f2, best_x2 = local_best
    if best_f2 < best_f:
        best_f, best_x = best_f2, best_x2

    return GlobalOptResult(
        x=np.asarray(best_x, float),
        fun=float(best_f),
        success=True,
        message=f"Finished (DE → {'BH → ' if do_basinhop else ''}{local_method})",
        nfev=int(total_evals),
        history=history,
        meta={"parallel": parallel, "workers": w, "dim": dim,
              "options": {"de_popsize": de_popsize, "de_maxiter": de_maxiter,
                          "bh_niter": bh_niter, "local_maxiter": local_maxiter}}
    )


#%%
if __name__ == "__main__":
    
    # Noisy, multi-modal Rastrigin-like objective
    def obj(x):
        A = 10
        return A * len(x) + np.sum(x**2 - A * np.cos(2 * np.pi * x))
    
    bounds = [(-5.12, 5.12)] * 5
    
    res = hybrid_global_optimize(
        obj, bounds,
        parallel='threads',
        de_popsize=12, de_maxiter=300, workers=-1,   # use all CPUs
        do_basinhop=True, bh_niter=10, bh_stepsize=0.25,
        local_method="L-BFGS-B", polish_top_k=7,
        random_state=42, return_all=True
    )
    
    print("Best f:", res.fun)
    print("Best x:", res.x)
    print("Function evals:", res.nfev)
    print("Message:", res.message)
