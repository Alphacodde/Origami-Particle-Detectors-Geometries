# Mark 1 → Mark 2: Why the Gun Changed

## Mark 1 (forward cone, single-sided beam)

`PrimaryGeneratorAction` sampled a start point on a transverse disk in front
of the detector (`gunZmm`, `diskRadiusMm`), then fired inward with an
incidence angle θ drawn uniformly in solid angle up to a cone half-angle
(`randomAngleMaxDeg`, default 60°), azimuth φ uniform over the full 360°.

This models a **one-sided test-beam** scenario: particles approaching from
a single general direction, spread over a cone. It is a reasonable model
for a fixed-target/test-beam setup, or for validating a single sensor
module in isolation (see `validate_flat_plate.mac`).

It is **not** a good model for how a real barrel tracker layer is used.
HEP trackers (ALICE ITS3, CMS, ATLAS ITk, etc.) are wrapped around a
collision point and are hit by particles emerging from that vertex across
the *full* solid angle — forward, backward, and everywhere between — not
from one external cone. Testing a barrel-relevant fold structure with a
Mark-1-style forward cone silently restricts the angular acceptance study
to a slice of the angles the structure will actually see in the real
use case implied by the proposal (Sec 6.4's "angular acceptance", Sec 4.1's
framing against ALICE ITS3 / CMS Phase-2 / LHCb).

Symptom this produced: sparse/null bins at large θ in the angle scan
output (see the `path_X0_by_angle_bin` nulls at 45°/55° in early Kresling/
Yoshimura scans) — solid-angle-uniform sampling within a 60° cone puts very
few events at the high end of that cone to begin with, and hit efficiency
was already dropping there, so some bins landed on literally zero hits by
chance. That is a **statistics artifact of the sampling geometry**, not
evidence that the structure has no acceptance past 35–45°.

## Mark 2 (isotropic 4π vertex source)

`PrimaryGeneratorAction` now places the vertex at (or near, if a beamspot
smear is enabled) the origin — the detector's nominal center — and samples
the **direction** isotropically over the full 4π solid angle: cosθ uniform
in [-1, 1], φ uniform in [0, 2π]. This is the standard "particles emerging
from a common vertex, detector wrapped around it" model appropriate for a
barrel-tracker acceptance study.

Two angle quantities are now recorded per event, because they answer
different questions and the proposal's own material-budget formula
(Sec 4.2, `B(θ) = t(θ)/X₀`) genuinely wants the *local* one:

- `labThetaDeg` — polar angle of the primary's momentum direction relative
  to the fixed lab z-axis. This is the physics variable relevant to
  detector coverage vs. pseudorapidity η, and is what most naturally maps
  onto "angular acceptance" in the collider sense.
- `localIncidenceDeg` — angle between the primary's momentum direction and
  the local sensor-facet normal at first entry into silicon. This is the
  variable that actually drives path length via the `1/cosθ_local` slant
  factor, and is the physically correct quantity for the B(θ) material-
  budget curve.

For a flat single panel at normal orientation the two coincide; for a
folded, multi-facet structure they generally do not, and conflating them
was a latent issue in the Mark 1 analysis (`path_X0_by_angle_bin` was
implicitly binning by lab-frame angle while the write-up's framing implies
per-facet incidence).

## What did NOT change between Mark 1 and Mark 2

- `DetectorConstruction` (STL import, thickening validation, Kapton
  substrate, structure/fold-state metadata) — unchanged, still correct.
- `SensorSD` / hit recording — unchanged.
- `RunAction` filename/ntuple bookkeeping — unchanged, with two new
  ntuple columns appended (never inserted in the middle) so old analysis
  scripts that index by position for columns 0–10 still work; anything
  reading by name should pick up `labThetaDeg` / `localIncidenceDeg`
  directly.
- Physics list, particle type (pi+), momentum (5 GeV/c) — unchanged.
- Mark 1 macros are kept under `macros/mark1/` for anyone who wants the
  old single-sided-cone behavior (e.g. reproducing a test-beam-style
  cross-check) — they still work against Mark 2 binaries; Mark 2 only
  adds new messenger commands, it does not remove old ones except where
  noted below.

## Breaking change to be aware of

`randomAngleMaxDeg`, `gunZmm`, and `diskRadiusMm` are **no longer used** in
Mark 2's default isotropic mode — see `PrimaryGeneratorAction.hh` for the
new `/origami/gun/mode` switch (`isotropic4pi` vs `mark1cone`) if you need
to run the old mode for a specific comparison. Mark 1 macros still set
those commands; they are accepted (so old macros don't fatal out) but are
silently ignored unless `mark1cone` mode is selected. A startup log line
states which mode is active on every run so this can't go unnoticed in a
sweep.
