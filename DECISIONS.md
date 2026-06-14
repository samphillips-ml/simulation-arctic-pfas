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
advective export (see Section 8), computed diagnostically rather than prescribed --
though see Section 8 for an important caveat on what this comparison can and cannot
establish given this decision.

The largest single omitted sink -- shelf sediment sorption/burial -- would be spatially
concentrated on the same Arctic shelves (Kara, Laptev, East Siberian, Beaufort) where
riverine injection occurs and where in-situ observational coverage is weakest (the
Russian-sector data gap). The spatial pattern of this omission is therefore plausibly
correlated with the spatial pattern of the training data gap, which may limit the
bias correction's ability to learn this term even where it matters most.

---

## 5. Polar Cap Mask

**Decision:** Grid cells above 89N treated as land (land_mask = True).

**Previous value:** 88N. Raised to 89N after confirming the Transpolar Drift signal
was largely absent from the tracer field at 88N -- the mask was removing too much of
the polar cap and blocking a major Arctic transport pathway. At 89N the TPD signal is
visible in the tracer field. Numerical stability confirmed at 89N for the full run.

**Justification:** Regular lat/lon grids have a coordinate singularity at 90N where
zonal cell width dx -> 0. Despite velocity capping (CFL <= 1), tracer accumulates
in the polar convergence zone due to the velocity sink at the pole -- all currents
converge but none escape. Diagnostic output confirmed exponential blowup originating
at the pole in early test configurations, with the polar cap acting as a tracer trap.
Masking above 89N removes the singularity while preserving the TPD pathway. TOPAZ4
uses a bipolar grid natively to avoid this issue; the regular lat/lon interpolated
output re-introduces it.

**Limitation:** Removes a small fraction of domain area (<0.5%) immediately
surrounding the pole. No training observations exist above 89N so PINN training is
unaffected. Tipton 2025 expedition includes North Pole depth profiles whose exact
latitude is unknown; if any are above 89N they cannot be used for validation.
Representing the TPD pathway fully (i.e. removing the mask entirely) would require
implementing advection on TOPAZ4's native bipolar grid -- flagged as future work, not
pursued here given the scope of regridding the advection scheme relative to remaining
project time. The current 89N mask is judged to capture the TPD pathway adequately
for this purpose.

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

**Caveat (given Section 4 -- no sink terms):** Bias correction operates on the
concentration field at observation locations; it does not enforce a global mass
balance. With no sink terms, the simulation's total PFOA inventory has no removal
mechanism other than advective export at the open boundaries, so the basin-wide
inventory is expected to be systematically high. The bias-corrected concentration
field may still carry this elevated-inventory signature into the Fram Strait export
calculation, particularly in regions upstream of Fram Strait with sparse training
observations to constrain the correction. Agreement (or disagreement) between the
corrected simulated export and the Joerss 6.4 +/- 1.0 t/yr target should therefore be
interpreted primarily as a check on transport pathway/timing fidelity, conditional on
the no-sink inventory bias -- not as an independent validation of the corrected
field's absolute accuracy. Additionally, Joerss 2020 measures total (dissolved +
particulate/ice-associated) export; the simulation represents dissolved-phase
advection only, so sea-ice-mediated export (a real pathway not captured here) is a
further reason the two quantities may not be directly comparable even with a perfect
correction.

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

---

## 10. Advection Scheme -- First-Order Upwind

**Decision:** First-order upwind finite-difference advection, applied per layer
(Section 6), with velocity capping for CFL <= 1.

**Justification:** First-order upwind is numerically diffusive but unconditionally
monotonic -- it cannot produce negative concentrations or spurious overshoots near
sharp gradients (e.g. at river-mouth injection cells, Section 2). This matters
specifically because PFOA concentration is physically non-negative and is expected
to be log-transformed for PINN training (concentrations span orders of magnitude);
a scheme that occasionally produces small negative values near steep gradients would
be incompatible with that pipeline without ad hoc clipping.

Higher-order TVD schemes with flux limiters (e.g. MUSCL, Superbee, van Leer) are
more accurate for tracer transport and are arguably more "standard" for this class
of problem -- they substantially reduce numerical diffusion while retaining
monotonicity via the limiter. These were considered but not adopted given remaining
project time: the river injection box size (Section 2, BOX_HALF=16) was tuned against
upwind's diffusive behavior, and switching schemes this late would require rerunning
the full 264-month simulation and re-verifying stability (a less diffusive scheme
could behave differently near the injection boxes and polar cap mask, in ways not
yet characterized).

**Limitation:** Numerical diffusion from the upwind scheme is expected to smear sharp
tracer features (e.g. river plume filaments) over a wider area than physically
realistic, lowering peak concentrations and broadening gradients. This is a smooth,
spatially-coherent bias (more pronounced where gradients are steep, e.g. near river
mouths and the ice edge) and is treated as part of the systematic error the bias
correction is designed to absorb, consistent with the overall framing of the
simulation as an intentionally imperfect physics prior. TVD schemes are flagged as
future work if a later iteration of the simulation is undertaken.