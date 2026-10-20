# Frozen Analysis Specification (B0)

**Frozen:** 2026-10-18, before any new (post-Gate-1) analysis. Code: `experiments/common.py`. Plan: `refine-logs/EXPERIMENT_PLAN.md`.
Any change after this date is logged at the bottom with reason and date; no setting is chosen by maximum ICC.

## Data
- GazeBaseVR TEX, `data/collected/gazebasevr/full/S_<round><pid>_S<session>_<n>_TEX.csv`, 250 Hz, columns lx, ly, rx, ry (deg).
- RAN stays inside `gazebasevr.zip` until the B5 parameters below are frozen (they are, as of this file) and M4 starts.

## Event extraction (shared by all passes)
- Raw NaN mask (any eye, any of lx/ly/rx/ry) recorded **before** interpolation; linear interpolation of gaps ≤ 12 samples (inside only).
- Per-eye velocity: Savitzky–Golay, window 7, order 2, first derivative, dt = 4 ms, NaN → 0 before filtering; per-eye 2-D speed s_L, s_R.
- Detector (primary): I-VT on mean speed ½(s_L + s_R) > 40 °/s AND all four channels finite after interpolation; runs ≥ 3 samples.
- Event fields: onset s, offset e (exclusive end − 1 = e1), per-eye amplitudes from sample s−1 to e1, ampx = mean horizontal amplitude, per-eye PV = max speed over [s, e1], duration = (e1 − s)·4 ms.
- Missingness flag: any raw NaN in [s − 6, e1 + 20] (≈ −24 ms / +80 ms).

## Event pool, classes and abstention
- Exclusions: duration ≤ 150 ms and missingness flag false.
- Classes: RS = ampx < −15°; progressive = 0.5° < ampx < 10°; others discarded.
- l = log(PV_L / PV_R).
- Centring pool = RS ∪ progressive events of the recording after exclusions. Keep iff |l − median_pool(l)| < log 2.
- Count floor (primary): ≥ 5 RS and ≥ 20 progressive surviving events per recording; sensitivity floors RS ≥ 1 and RS ≥ 10 (progressive ≥ 20 throughout). Unfiltered rows apply the same floors to unfiltered counts.

## Marker and statistics
- A_r = ½[median_RS l − median_prog l] per recording (= participant-session).
- Round summary = mean of the two sessions; longitudinal comparisons require both sessions in both rounds.
- Primary: ICC(A,1) (two-way, absolute agreement, single). Secondary: ICC(C,1), Spearman, mean difference, SD of differences, Bland–Altman limits.
- CIs: participant bootstrap, 5000 resamples, shared indices across every compared method within a table; percentile intervals; seed 20261018.
- Every ΔICC is computed on intersection pairs (eligible under both compared methods).

## B2 random baseline
- 200 reps; per recording × class, keep a uniformly random subset of the same size as the centred rule keeps (seed 20261018 + rep). Random estimate = mean ICC(A,1) over reps; recomputed on each bootstrap draw's resampled pairs.

## B3a synthetic failure injection (frozen grid)
- Targets: R1 TEX recordings. Placement from the baseline pass. Injected eye uniform L/R. Seed per (condition, file) = crc32(condition + file).
- Per recording per condition: up to 3 RS + 3 progressive target events (F2, F3, N1) or up to 6 fixation intervals ≥ 300 ms (F1), injections ≥ 500 ms apart.
- F1 one-eye pseudo-saccade: at fixation start + 20 samples, horizontal step a ∈ {2, 4, 8}°, raised-cosine rise over 5 samples, then exponential return τ = 150 ms for 5τ (750 ms), then 0.
- F2a slow-down: injected eye's H and V on [s − 1, e1 + 1] time-stretched by f ∈ {1.5, 2, 3} (linear interpolation of the original segment, start-anchored; the stretched segment overwrites the following samples), so PV ≈ PV/f with amplitude preserved.
- F2b glitch: + b ∈ {3, 6}° on 3 consecutive horizontal samples of the injected eye centred on its speed peak within [s, e1].
- F3 dropout: injected eye set to NaN for k ~ U{3..12} samples starting at s + U{−10..10}.
- N1 coherent mild asymmetry: injected eye x' = x + (κ − 1)(x − x[s − 1]) on [s, e1] for H and V, then the end offset decays exponentially (τ = 150 ms) for 5τ, then 0; κ ∈ {0.85, 1.15}.
- Contamination: F2a f = 3 on a random fraction p ∈ {.05, .10, .20, .30} of a recording's eligible RS/prog events (no per-recording cap).
- Recovery: 5% of R1 participants (seeded), one random session each; F2a f = 3 on 30% of eligible events + F1 a = 4° on 10% of fixation intervals ≥ 300 ms.
- Matching: injected window ([s − 2, e1 + 2] for event faults; [c, c + 8] for F1) ↔ re-detected event with greatest sample overlap. Untouched events = events ≥ 25 samples from every injection window, matched to baseline by onset (±1 sample).
- Criteria (C2): abstention recall among admitted ≥ .80 with participant-cluster LB > .70 for F1 a ∈ {4, 8}, F2a f = 3, F2b b = 6°; N1 excess rejection ≤ 5 pp (UB < 10 pp); untouched decision flips ≤ 1%; recovery |ICC_centred(inj) − ICC_centred(base)| ≤ .05 and centred − unfiltered ΔICC CI > 0 on injected data. Floors (F1 2°, F2a f ∈ {1.5, 2}, F2b 3°) reported, not judged. Admitted rate < 30% ⇒ inconclusive for that fault.

## B4 per-eye matched detector
- Each eye: s_eye > 40 °/s & finite, runs ≥ 3 samples. Candidate pair = overlapping intervals AND |onset gap| ≤ 3 samples. 1:1 assignment by greatest overlap, then smallest onset gap; exact ties discarded. Event interval = union; amplitudes over union; per-eye PV over own interval. Same exclusions/classes/abstention/floors.
- Independent agreement: a mean-speed event is "matched" iff some 1:1 per-eye pair overlaps it. Criterion: OR(unmatched | rejected vs retained) with participant-cluster LB > 2; Cohen κ reported.

## B6 gain grid
- Per eye about the recording median of that eye's channel: x' = med + g (x − med) for H and V; g_L = q e^{δ/2}, g_R = q e^{−δ/2}; q ∈ {.9, 1, 1.1}, δ ∈ {−.30, −.15, 0, .15, .30}; plus session-asymmetric (q = 1; S1 δ = +.15 / S2 δ = −.15, and reverse). Full re-extraction. Criteria on |δ| ≤ .15 (incl. asymmetric): median |ΔA| ≤ .01, P95 ≤ .03, |ΔICC| ≤ .05, coverage loss ≤ 5 pp, all on records/pairs eligible in both runs.

## B5 matched-direction controls (frozen before RAN is unzipped)
- TEX: |amp_x| 1–10°, angle ≤ 15° from horizontal; within-recording greedy matching left↔right (Δ|amp| ≤ 1°, Δmid-x ≤ 2°, Δmid-y ≤ 1°, reversed-angle Δ ≤ 10°).
- RAN: first response 80–600 ms after target jump, |amp_x| 3–25°, angle ≤ 30°, direction agrees with target; matching Δ|amp| ≤ 2°, target-mid Δx ≤ 3°, Δy ≤ 2°, angle ≤ 10°; bands 3–10° and 15–25°.
- Same duration/missingness exclusions; no ratio filter in primary B5 (centred variant = sensitivity, centring within task pool).
- D_j = ½(l_left − l_right); session median; participant D = mean of the two R1 session medians; eligibility ≥ 5 pairs total. Transfer: Spearman ρ(D, two-session-mean TEX A_r), ≥ 100 participants. RAN S1–S2 ICC(A,1) on participants with ≥ 5 pairs per session, computed before transfer. Sensitivity: horizontal-component PV. No caliper widening.

## Change log
- (none)
