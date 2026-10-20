# Frozen Analysis Specification (B0)

**Frozen:** 2026-10-18, before any new (post-Gate-1) analysis. Code: `experiments/common.py`. Plan: `refine-logs/EXPERIMENT_PLAN.md`.
Any change after this date is logged at the bottom with reason and date; no setting is chosen by maximum ICC.

## Data
- GazeBaseVR TEX, `data/collected/gazebasevr/full/S_<round><pid>_S<session>_<n>_TEX.csv`, 250 Hz, columns lx, ly, rx, ry (deg).
- RAN stays inside `gazebasevr.zip` until the B5 parameters below are frozen (they are, as of this file) and M4 starts.

## Event extraction (shared by all passes)
- Raw NaN mask (any eye, any of lx/ly/rx/ry) recorded **before** interpolation; linear interpolation inside the recording with at most 12 consecutive filled samples per gap (pandas `interpolate(limit=12, limit_area='inside')`: longer gaps are partially filled from their start) [A1].
- Per-eye velocity: Savitzky–Golay, window 7, order 2, first derivative, dt = 4 ms, NaN → 0 before filtering; per-eye 2-D speed s_L, s_R.
- Detector (primary): I-VT on mean speed ½(s_L + s_R) > 40 °/s AND all four channels finite after interpolation; runs ≥ 3 samples.
- Event fields: onset s, run exclusive end e; e1 = min(e, n−1), i.e. the first sub-threshold sample [A1]. Per-eye amplitudes from sample s−1 to e1, ampx = mean horizontal amplitude, per-eye PV = max speed over [s, e1], duration = (e1 − s)·4 ms.
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
- Targets: R1 TEX recordings. Placement from the baseline pass. Injected eye uniform L/R. Seed per (condition, file) = crc32(condition + "|" + basename(file)) [A6].
- Per recording per condition: up to 3 RS + 3 progressive target events (F2, F3, N1) or up to 6 fixation intervals ≥ 300 ms with no raw NaN between consecutive detected events (F1) [A6], injections ≥ 500 ms apart across both classes [A6].
- F1 one-eye pseudo-saccade: at fixation start + 20 samples (fixation start = preceding event's e1, the first sub-threshold sample [A6]), horizontal step a ∈ {2, 4, 8}°, raised-cosine rise over 5 samples, then exponential return τ = 150 ms for 5τ (750 ms), then 0.
- F2a slow-down: injected eye's H and V on [s − 1, e1 + 1] time-stretched by f ∈ {1.5, 2, 3} (linear interpolation of the original segment, start-anchored; the stretched segment overwrites the following samples), so PV ≈ PV/f with amplitude preserved.
- F2b glitch: + b ∈ {3, 6}° on 3 consecutive horizontal samples of the injected eye centred on its speed peak within [s, e1].
- F3 dropout: injected eye set to NaN for k ~ U{3..12} samples starting at s + U{−10..10}.
- N1 coherent mild asymmetry: injected eye x' = x + (κ − 1)(x − x[s − 1]) on [s, e1] for H and V, then the end offset decays exponentially (τ = 150 ms) for 5τ, then 0; κ ∈ {0.85, 1.15}.
- Contamination: F2a f = 3 on a random fraction p ∈ {.05, .10, .20, .30} of a recording's eligible RS/prog events (no per-recording cap).
- Recovery: 5% of R1 participants (seeded), one random session each; F2a f = 3 on 30% of eligible events + F1 a = 4° on 10% of fixation intervals ≥ 300 ms.
- Matching: injected window ([s − 2, e1 + 2] for event faults; [c, c + 8] for F1) ↔ re-detected event with greatest sample overlap. Untouched events = events ≥ 25 samples from every injection window and every modified sample range (incl. F1/N1 exponential returns and F2a overwrites) [A3], matched to baseline by onset (±1 sample).
- Criteria (C2): abstention recall among admitted ≥ .80 with participant-cluster LB > .70 for F1 a ∈ {4, 8}, F2a f = 3, F2b b = 6°; N1 excess rejection ≤ 5 pp (UB < 10 pp), defined as the paired change in rejection on the same admitted targets (N1 pass vs baseline pass), participant-cluster bootstrap; admission loss reported separately [A4]; untouched decision flips ≤ 1%; recovery |ICC_centred(inj) − ICC_centred(base)| ≤ .05 and centred − unfiltered ΔICC CI > 0 on injected data. Floors (F1 2°, F2a f ∈ {1.5, 2}, F2b 3°) reported, not judged. Admitted rate < 30% ⇒ inconclusive for that fault.

## B4 per-eye matched detector
- Each eye: s_eye > 40 °/s & that eye's own channels finite after interpolation [A2], runs ≥ 3 samples. Candidate pair = overlapping intervals AND |onset gap| ≤ 3 samples. 1:1 assignment by greatest overlap, then smallest onset gap; exact ties discarded. Event interval = union; amplitudes over union; per-eye PV over own run [start, end) (exclusive) [A2]. Same exclusions/classes/abstention/floors.
- Independent agreement: a mean-speed event is "matched" iff some 1:1 per-eye pair overlaps it. Criterion: OR(unmatched | rejected vs retained) with participant-cluster LB > 2; Cohen κ reported.

## B6 gain grid
- Per eye about the recording median of that eye's channel: x' = med + g (x − med) for H and V; g_L = q e^{δ/2}, g_R = q e^{−δ/2}; q ∈ {.9, 1, 1.1}, δ ∈ {−.30, −.15, 0, .15, .30}; plus session-asymmetric (q = 1; S1 δ = +.15 / S2 δ = −.15, and reverse). Full re-extraction. Criteria on |δ| ≤ .15 (incl. asymmetric): median |ΔA| ≤ .01, P95 ≤ .03, |ΔICC| ≤ .05 on records/pairs eligible in both runs; coverage loss ≤ 5 pp = (eligible_base − eligible_cond) / all input recordings [A5]. C4 is reported INCOMPLETE unless all 11 judged conditions have results [A5].

## B5 matched-direction controls (frozen before RAN is unzipped)
- TEX: |amp_x| 1–10°, angle ≤ 15° from horizontal; within-recording greedy matching left↔right (Δ|amp| ≤ 1°, Δmid-x ≤ 2°, Δmid-y ≤ 1°, reversed-angle Δ ≤ 10°).
- RAN: first response 80–600 ms after target jump, |amp_x| 3–25°, angle ≤ 30°, direction agrees with target; matching Δ|amp| ≤ 2°, target-mid Δx ≤ 3°, Δy ≤ 2°, angle ≤ 10°; bands 3–10° and 15–25°.
- Same duration/missingness exclusions; no ratio filter in primary B5 (centred variant = sensitivity, centring within task pool).
- D_j = ½(l_left − l_right); session median; participant D = mean of the two R1 session medians; eligibility ≥ 5 pairs total. Transfer: Spearman ρ(D, two-session-mean TEX A_r), ≥ 100 participants. RAN S1–S2 ICC(A,1) on participants with ≥ 5 pairs per session, computed before transfer. Sensitivity: horizontal-component PV. No caliper widening.

## Change log
All amendments dated 2026-10-18, after external code review (Codex gpt-6.1-sol, trace `.aris/traces/code-review/2026-10-18_run01/`) and **before any full-data pass**; only a 1/8-participant smoke subset had been run. None was chosen by looking at an ICC. Pre-amendment text: `experiments/SPEC_20261018_frozen_v0.md`.
- A1 (extraction contract = pilot implementation): the frozen text described e1 = e − 1 and whole-gap interpolation, but the pilot that produced the reported exploratory numbers used e1 = e and pandas partial-fill interpolation. The contract is versioned to what the pilot ran, so the base pass is comparable to the pilot; text corrected.
- A2 (B4 bug fix): per-eye detector used the binocular finite mask and an inclusive end in own-eye PV; now own-eye mask and exclusive own-run PV, as the B4 text intends.
- A3 (B3a untouched controls): untouched specificity now excludes events near any modified sample, not only near the short matching windows.
- A4 (N1 excess rejection): previous code subtracted the untouched flip rate, which is not the targets' baseline rejection; now the paired before/after change on the same targets.
- A5 (B6 coverage / completeness): coverage loss in percentage points of all input recordings (per-pass input manifest, identical across conditions); C4 cannot pass with missing judged conditions.
- A6 (implementation details made explicit / corrected): seed serialization documented; F1 candidates exclude tracking-loss gaps (previously could log an injection into NaN); 500 ms spacing enforced across RS and progressive targets jointly.
- Quality-covariate CIs: participant-cluster bootstrap with 5000 shared draws (was recording bootstrap, 2000) — conforms to the existing statistics section.

- D1 (2026-10-18, disclosed after external review): N1 as implemented scales within-event displacement by κ and adds a decaying end offset (above). The plan described an amplitude-preserving time warp. The low N1 rejection rates therefore validate this mild displacement-scaling perturbation, not a time warp.
- D2 (2026-10-18, disclosed after external review): the ≥ 500 ms spacing rule is applied only to the k-per-class paths (F2, F3, N1) and F1. The fractional paths (C_p, REC) draw events without spacing. Also, F1/N1 modified tails last ≈ 750 ms, so spacing onsets does not guarantee non-overlapping modifications.
