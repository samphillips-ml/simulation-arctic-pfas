# Simulation Design Decisions
# Arctic PFOA Transport Simulation -- Methods Reference
# Updated: June 2026

---

## 1. Southern Boundary Concentrations (60N)

**Decision:** Atlantic sector 0.032 ng/L, Pacific sector 0.025 ng/L, zero below 200m.

**Previous values:** Atlantic 0.093 ng/L, Pacific 0.083 ng/L (Benskin 2012, Li 2018).
Those were point measurements not representative of the full 60N boundary, and produced
boundary inflow of ~45,000 kg/yr PFOA -- 3x the riverine total and inconsistent with
the Joerss 2020 Fram Strait export budget.

**Justification:** Yeung et al. 2017 (doi:10.1021/acs.est.7b00788) report mean PFOA
concentration in the Arctic polar mixed layer of 32 +/- 15 pg/L (0.032 ng/L).
This represents the interior Arctic concentration set by Atlantic inflow, and is used
here as the Atlantic boundary value. Pacific sector set slightly lower (0.025 ng/L)
consistent with Bering Strait observations being somewhat lower than Atlantic corridor.
Boundary applied as static climatology (no temporal variation); temporal variability
delegated to PINN bias correction.

**Limitation:** Boundary concentration varies spatially and temporally along 60N.
Applying a uniform sector-mean underrepresents the Norwegian Atlantic Current inflow
pathway and overrepresents low-concentration regions. A time-varying spatially resolved
boundary from a global PFAS model (e.g. Zhang et al. 2017) would be preferable but
is not available operationally for PFOA.

---

## 2. Riverine Injection Box Size

**Decision:** BOX_HALF = 16 (33x33 cell injection box, ~4 deg x 4 deg).

**Previous value:** BOX_HALF = 2 (5x5 box, ~0.6 deg x 0.6 deg).

**Justification:** At 0.125 deg resolution the Ob and Yenisei deltas are effectively
single cells. The 5x5 box produced peak mouth concentrations of ~11 ng/L/month at
the Ob (June peak), ~200x observed values, causing numerical blowup by year 2.
A 33x33 box (~400km x 400km) represents the coastal mixing zone across the broad
Kara Sea shelf and brings peak injection to ~0.37 ng/L/month, consistent with
near-mouth observations. This is standard practice in coarse-resolution ocean tracer
models where river plume mixing is sub-grid scale.

**Limitation:** The 33x33 box spreads riverine signal further than the actual mixing
zone, potentially diluting near-mouth signal. No direct PFOA measurements exist at
Ob, Yenisei, Lena, or Mackenzie river mouths to constrain this parameter.

---

## 3. Riverine Flux Magnitude

**Decision:** 14,800 kg/yr total pan-Arctic PFOA, partitioned by ArcticGRO discharge
fraction, with ArcticGRO monthly seasonal scaling.

**Source:** Stemmler & Lammel 2010 (doi:10.5194/acp-10-9965-2010).

**Justification:** No direct mainstem measurements exist for PFOA flux from Mackenzie,
Ob, Yenisei, or Lena. Stemmler & Lammel 2010 provides the only pan-Arctic modeled
estimate. Yeung et al. 2017 suggests the lower bound of riverine estimates fits Arctic
observations better; 14,800 kg/yr is treated as an upper bound. Sensitivity to this
parameter is high and should be reported.

**Limitation:** Modeled estimate with large uncertainty. No validation data available
at river mouths. Seasonal scaling via ArcticGRO discharge fractions assumes PFOA
concentration is uniform year-round (only discharge varies seasonally).

---

## 4. Sink Terms -- First-Order Decay

**Decision:** First-order decay with k = ln(2) / (5 yr * 365.25 * 86400 s) applied
each 6-hour timestep as a numerical stabilization term.

**Justification:** PFOA does not measurably degrade in seawater (hydrolysis half-life
>92 years; Scheringer et al. 2012). Real removal mechanisms (sediment burial, sea ice
partitioning, advective export) are not parameterized due to lack of Arctic-specific
rate constants. The decay term is a numerical stabilization prior, not a physical
process representation. Its effect is subsumed into PINN bias correction. This follows
the Ice-BCNet approach (Yuan et al. 2024, doi:10.1016/j.ocemod.2024.102326) of
delegating unresolved sink processes to the learned correction.

**Limitation:** Not physically motivated. Half-life of 5 years is arbitrary. Sensitivity
analysis across half-life values (2, 5, 10, 50 yr) should be reported. Fram Strait
net export of 6.4 +/- 1.0 t/yr PFOA (Joerss et al. 2020) used as Phase 3 validation
target for simulated advective export rather than as a prescribed sink.

---

## 5. Polar Cap Mask

**Decision:** Grid cells above 88N treated as land (land_mask = True).

**Justification:** Regular lat/lon grids have a coordinate singularity at 90N where
zonal cell width dx -> 0. Despite velocity capping (CFL <= 1), tracer accumulates
in the polar convergence zone due to the velocity sink at the pole -- all currents
converge but none escape. Diagnostic output confirmed exponential blowup originating
at 89-90N starting month 13, with the polar cap acting as a tracer trap. Masking
above 88N removes the singularity. TOPAZ4 uses a bipolar grid natively to avoid
this issue; the regular lat/lon interpolated output re-introduces it.

**Limitation:** Removes ~0.5% of domain area. No training observations exist above
88N so PINN is unaffected. Tipton 2025 expedition includes North Pole depth profiles
whose latitude is unknown; if any are above 88N they cannot be used for validation.
Future work should implement advection on TOPAZ4 native bipolar grid.

---

## 6. Vertical Structure -- Two-Layer Scheme

**Decision:** Layer 1: 0-200m (depth-weighted mean velocity). Layer 2: 200-4000m
(depth-weighted mean velocity). Vertical exchange via diapycnal diffusion kappa_v
= 1e-5 m2/s.

**Justification:** Yeung et al. 2017 report PFOA detectable only above 150m depth
in the Arctic polar mixed layer. A two-layer split at 200m captures the primary
transport layer while allowing vertical exchange with the deep reservoir. Full 3D
output avoided (~270GB storage) in favor of monthly layer-mean 2D fields sufficient
for PINN training against sparse observations.

---

## 7. Atmospheric Deposition

**Decision:** MacInnis et al. 2017 (doi:10.1039/c6em00593d) Table S10 Devon Ice Cap
values 1993-2007, held constant at 2004-2007 mean (~72 ng/m2/yr) post-2007.

**Justification:** Only available multi-decadal Arctic atmospheric PFOA deposition
record. Post-2007 constant follows Yeung et al. 2017 precedent. Atmospheric input
~1 t/yr total, minor relative to riverine (~14.8 t/yr) and boundary inflow.

---

## 8. Fram Strait Export -- Validation Target Only

**Decision:** Joerss et al. 2020 net PFOA export through Fram Strait (6.4 +/- 1.0
t/yr) used as Phase 3 validation target, not prescribed as a sink.

**Justification:** Joerss 2020 is a single-year point measurement (2018), not a
climatological constraint. Prescribing it as a boundary condition would be circular
-- it is an emergent result of the circulation, not an independent forcing. Instead,
simulated PFOA flux through Fram Strait cells (78-80N, 10W-10E) will be computed
from C * vy * dx * dz and compared to 6.4 t/yr as a diagnostic validation check
in Phase 3.