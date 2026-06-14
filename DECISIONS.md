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

## 4. Sink Terms -- None Applied (Conservative Tracer)

**Decision:** No sink terms applied. PFOA treated as a fully conservative tracer
(HALF_LIFE_YR = 99.0, effectively no decay).

**Previous approach (superseded):** First-order decay with a 5-year half-life applied
as a "numerical stabilization term."

**Justification:** PFOA does not measurably degrade in seawater (hydrolysis half-life
>92 years; Scheringer et al. 2012), consistent with treating it as conservative. Real
removal mechanisms (sediment burial, sea ice partitioning, advective export) lack
Arctic-specific rate constants and are not parameterized. Rather than impose an
arbitrary numerical decay term to mask this, all sink processes are omitted entirely
and delegated to PINN bias correction -- following the Ice-BCNet approach (Yuan et al.
2024, doi:10.1016/j.ocemod.2024.102326) of delegating unresolved processes to the
learned correction, applied here as a clean omission rather than a placeholder decay
rate.

**Limitation:** Simulation will tend to overestimate concentrations relative to
observations in regions/depths where real removal processes are significant, since
nothing in the model removes mass except advective export at boundaries. This is
expected and intentional; the bias correction is responsible for learning the
spatial/temporal pattern of this overestimation. Fram Strait net export of 6.4 +/- 1.0
t/yr PFOA (Joerss et al. 2020) remains the Phase 3 validation target for simulated
advective export (see Section 8), computed diagnostically rather than prescribed.

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

## 7. Atmospheric Deposition -- Not Represented

**Decision:** No atmospheric deposition term included.

**Previous approach (superseded):** MacInnis et al. 2017 (doi:10.1039/c6em00593d)
Table S10 Devon Ice Cap values 1993-2007, held constant at 2004-2007 mean (~72
ng/m2/yr) post-2007, characterized as "minor relative to riverine (~14.8 t/yr)."

**Justification:** A single ice core record (Devon) extrapolated pan-Arctic and held
constant for ~18 years is not a defensible spatially- or temporally-resolved flux
estimate, and the "minor" characterization does not hold up against independent
evidence: Yeung et al. 2017 (doi:10.1021/acs.est.7b00788) modeled atmospheric inputs
as accounting for 34-59% (~11-19 pg/L) of measured PFOA concentration in the polar
mixed layer (mean 32 +/- 15 pg/L). This is a substantial fraction of surface-layer
PFOA, not a second-order term -- so including a weak, single-source deposition
estimate would add a spatially uniform term whose own magnitude is highly uncertain,
without clearly improving the prior. The residual (obs - sim) already reflects the
true missing atmospheric contribution whether or not a deposition term is explicitly
modeled. Atmospheric deposition is therefore omitted entirely, consistent with the
Section 4 sink-term philosophy: processes lacking a defensible flux estimate are
delegated to bias correction rather than included as placeholders.

**Limitation:** Atmospheric deposition is a known, literature-supported, but currently
unquantifiable input, expected to contribute to underestimation of surface-layer
(0-200m) PFOA concentrations and is a plausible partial explanation for the observed
basin-mean underestimate (~2-2.5x). Should be discussed explicitly in the limitations
section. If PINN residuals show depth- or season-structured patterns consistent with
atmospheric input (e.g. concentrated in the surface layer, correlated with
sea-ice/precipitation seasonality), this would be indirect evidence the omission
matters -- a citable finding in its own right.

**Note:** Sections 4 and 7 pull in opposite directions on the basin-mean bias: no
sink terms pushes the simulation toward overestimation, while no atmospheric
deposition pushes toward underestimation (especially in the surface layer). These
two omissions partially offset in their net effect on the mean, though they operate
through different mechanisms and at different depths -- worth stating explicitly in
the methods/limitations text.

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

---

### 9a. River Mouth Coordinates and TOPAZ4 Mask Verification

All ten river injection coordinates were checked against the TOPAZ4 land
mask (0.125 deg resolution) using check_dvina_mask.py and
check_new_river_mouths.py. At this resolution, narrow deltas and bays are
frequently masked as land; where the natural delta/bay coordinate fell on
land, the nearest wet cell was used instead. All shifts are well within
the BOX_HALF=16 (~2 deg) injection box half-width (sec. 2), so the
injected signal still lands within the intended shelf region regardless
of the shift.

| River          | Coordinate used      | Notes                                          | Shift from natural mouth |
|-----------------|----------------------|-------------------------------------------------|---------------------------|
| Lena            | 74.000N, 126.500E    | Laptev Sea outlet (original config)             | --                          |
| Yenisei         | 73.750N,  82.000E    | Kara Sea outlet (original config)               | --                          |
| Ob              | 72.500N,  73.625E    | Kara Sea outlet (original config)               | --                          |
| Mackenzie       | 70.000N, -134.000E   | Beaufort Sea outlet (original config)           | --                          |
| Pechora         | 69.125N,  54.000E    | Pechora Sea, Barents                            | 0.13 deg                    |
| Yana            | 71.750N, 136.500E    | Yana Bay, Laptev Sea                            | 0.25 deg                    |
| Olenek          | 73.500N, 122.500E    | Olenek Bay, Laptev Sea (west of Lena delta)     | 0.71 deg                    |
| Kolyma          | 70.000N, 162.125E    | Kolyma Gulf, East Siberian Sea                  | 0.80 deg                    |
| Indigirka       | 71.750N, 151.750E    | East Siberian Sea, off delta                    | 1.21 deg                    |
| NorthernDvina   | 70.500N,  38.250E    | Southernmost Barents cell; White Sea masked     | n/a -- see sec. 9 note      |

**Indigirka note:** the 1.21 deg shift is the largest of the five new
rivers. The natural delta coordinate (~71.0N, 150.8E) and surrounding area
are masked as land in TOPAZ4 -- likely reflecting a broad shallow/deltaic
region not resolved at 0.125 deg. The nearest wet cell at 71.75N, 151.75E
is still on the East Siberian Sea shelf and within the injection box
range of the intended outlet, so the spatial intent (ESS shelf signal) is
preserved.

**Olenek note:** placed west of the Lena delta (73.5N, 122.5E) vs Lena's
own mouth at 74.0N, 126.5E -- a separation of roughly 4 deg longitude,
sufficient that the two rivers' 33x33 injection boxes do not fully
overlap, preserving them as distinct spatial sources.

No formal literature citation is required for mouth/bay coordinates
themselves (geographic fact, not a contested estimate). The methods
section should state that all injection coordinates were verified against
the TOPAZ4 land mask, with the table above available as supplementary
detail if needed.