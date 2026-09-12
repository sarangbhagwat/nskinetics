# Justification for `anaerobic_growth_mult = 0.75`

**Model:** `s_cerevisiae_ferm_fb_inhib_mod_ibo` (Antimony source
`nskinetics/models/s_cerevisiae_ferm_fb_inhib_mod_ibo/s_cerevisiae_ferm_fb_inhib_mod_ibo_antimony.txt`)
**Date:** 2026-09-12 · **Probed at** `eb34820` (dev); the shipped default was changed from `1.0` to `0.75` in the same change that adds this report

## 1. What the parameter does

Growth on glucose (r7) is the only biomass-forming reaction that is allowed to run
without oxygen. Its rate law is

```
r7 = (anaerobic_growth_mult + (1 - anaerobic_growth_mult) * f_O2) * v7_0
```

where `v7_0` is the unscaled aerobic growth rate and `f_O2` drops to 0 when
`is_aerobic` is cleared at the end of stage 1. `anaerobic_growth_mult` is therefore
the **ratio of the anaerobic to the aerobic specific growth rate on excess glucose**,
μ_anaerobic / μ_aerobic. The stoichiometry of r7 (0.732 g biomass per g glucose) is
fixed, so the parameter scales a rate, not a yield. The shipped default of 1.0 asserts
that yeast grows exactly as fast without oxygen as with it.

## 2. Literature

**Specific growth rate.** In their review of anaerobic yeast cultivation, Hakkaart et
al. (2021) state that the maximum specific growth rate of *S. cerevisiae* in sterol- and
unsaturated-fatty-acid (UFA)-supplemented, glucose-grown anaerobic batch cultures is
typically about 25 % lower than in corresponding aerobic cultures [1]. That is a
ratio of ≈ 0.75, and it is precisely the quantity this parameter encodes.

**Absolute rates.** The reference laboratory strain CEN.PK113-7D grows at roughly
0.37–0.40 h⁻¹ aerobically on glucose [2]. Fully anaerobic batch cultures of the same
strain supplemented with ergosterol and Tween 80 grow at about 0.28–0.33 h⁻¹, consistent
with the 25 % reduction above. With ergosterol but **no** UFA supplement, Dekker et al.
(2019) measured 0.14–0.20 h⁻¹ over consecutive anaerobic bioreactor batches [3], i.e. a
ratio of only 0.4–0.5 relative to aerobic growth. An unsupplemented industrial medium is
closer to that lower case.

**Why anaerobic growth is slower and finite.** Sterol and UFA biosynthesis are
O₂-dependent (12 mol O₂ per mol ergosterol; the Ole1 desaturase for UFAs), which is why
strictly anaerobic growth requires supplementation [3, 4]. Cells leaving an aerobic
stage carry sterol/UFA reserves and keep dividing for a limited number of generations,
then stall. Growth in the anaerobic stage should therefore be slower than aerobic growth
and should not be treated as unconstrained.

**Biomass yield (consistency check).** Verduyn et al. (1990) report a maximal anaerobic
biomass yield of 0.10 g g⁻¹ glucose for *S. cerevisiae* CBS 8066 [5], versus ≈ 0.5 g g⁻¹
for respiratory aerobic growth. The multiplier does not set the yield directly, but the
model's anaerobic yield should remain of this order.

## 3. Model behaviour across candidate values

The shipped model was probed at `eb34820` in two settings: an anaerobic 100 g/L glucose
batch (x₀ = 0.1 g/L, no feed spikes, `is_aerobic = 0`), and the shipped fed-batch process
(`create_sugar_prep_and_fermentation_system`, aeration off at 5 g/L biomass).

| `anaerobic_growth_mult` | Batch anaerobic μ_max (h⁻¹) | Batch Y_x/s (g g⁻¹) | Process τ (h) | Harvest x (g/L) | Harvest EtOH (g/L) | Glucose spikes |
| --- | --- | --- | --- | --- | --- | --- |
| 1.0 (old default) | 0.69 | 0.103 | 42.9 | 15.7 | 118.0 | 4 |
| **0.75** | **≈ 0.49** | **≈ 0.08** | **61.8** | **13.0** | **114.8** | **3** |
| 0.5 | 0.29 | 0.053 | 66.0 | 10.2 | 106.9 | 1 |

For reference, the aerobic glucose-excess batch gives μ_max = 0.69 h⁻¹ at the shipped
`k_7`; the values at 0.75 are interpolated from probes at 0.8 (0.53 h⁻¹, 0.084 g g⁻¹) and
0.7 (0.45 h⁻¹, 0.074 g g⁻¹).

Two conclusions follow:

1. **The parameter must be applied as a ratio, not tuned to an absolute anaerobic rate.**
   The model's aerobic μ_max (0.69 h⁻¹) already exceeds the literature by ≈ 1.6×. Forcing
   the anaerobic rate down to the literature's ≈ 0.3 h⁻¹ would need a multiplier near
   0.5, which halves the anaerobic biomass yield to 0.05 g g⁻¹ because glycolysis (r1)
   keeps consuming glucose at its own rate while r7 is slowed. The absolute rate is a
   `k_7` calibration question and is out of scope here.
2. **The yield check supports the upper part of the plausible range.** At 0.75 the
   anaerobic yield (≈ 0.08 g g⁻¹) remains within the scatter of reported anaerobic yields
   [5]; at 0.5 it is clearly too low.

## 4. Decision

- **Default:** `anaerobic_growth_mult = 0.75`, the literature ratio for sterol/UFA-
  sufficient anaerobic growth [1].
- **Uncertainty range for downstream analyses:** 0.5–0.8. The lower bound represents
  sterol-/UFA-poor media [3]; the upper bound leaves headroom above the central estimate.
- **Consequence:** the change moves the shipped fed-batch baseline substantially
  (τ 42.9 → 61.8 h, one fewer glucose spike, ≈ 3 g/L less ethanol at harvest). This is an
  intended model change, not a regression, and the downstream isobutanol smoke-test pins
  must be re-pinned before the validation gate passes again.

## References

1. Hakkaart, X., Liu, Y., Hulst, M., El Masoudi, A., Peuscher, E., Pronk, J., van Gulik,
   W., Daran-Lapujade, P. (2021). Critical parameters and procedures for anaerobic
   cultivation of yeasts in bioreactors and anaerobic chambers. *FEMS Yeast Research*
   21(5), foab035. https://doi.org/10.1093/femsyr/foab035
2. van Dijken, J. P., et al. (2000). An interlaboratory comparison of physiological and
   genetic properties of four *Saccharomyces cerevisiae* strains. *Enzyme and Microbial
   Technology* 26(9–10), 706–714. https://doi.org/10.1016/S0141-0229(00)00162-9
3. Dekker, W. J. C., Wiersma, S. J., Bouwknegt, J., Mooiman, C., Pronk, J. T. (2019).
   Anaerobic growth of *Saccharomyces cerevisiae* CEN.PK113-7D does not depend on
   synthesis or supplementation of unsaturated fatty acids. *FEMS Yeast Research* 19(6),
   foz060. https://doi.org/10.1093/femsyr/foz060
4. Andreasen, A. A., Stier, T. J. B. (1953). Anaerobic nutrition of *Saccharomyces
   cerevisiae*. I. Ergosterol requirement for growth in a defined medium. *Journal of
   Cellular and Comparative Physiology* 41(1), 23–36.
   https://doi.org/10.1002/jcp.1030410103
5. Verduyn, C., Postma, E., Scheffers, W. A., van Dijken, J. P. (1990). Physiology of
   *Saccharomyces cerevisiae* in anaerobic glucose-limited chemostat cultures. *Journal
   of General Microbiology* 136(3), 395–403.
   https://doi.org/10.1099/00221287-136-3-395
6. Pham, H. T. B., Larsson, G., Enfors, S.-O. (1998). Growth and energy metabolism in
   aerobic fed-batch cultures of *Saccharomyces cerevisiae*: simulation and model
   verification. *Biotechnology and Bioengineering* 60(4), 474–482. (Source of the
   r1–r11 kinetic structure.)
