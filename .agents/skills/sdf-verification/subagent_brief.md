# Subagent Brief Template — Parallel SDF Discovery

Copy this into the scratch directory, fill in the `<...>` fields, and point each subagent at it.
The parent agent owns integration: README rows, `sdf_gems.glsl`, independent re-verification, and one commit per primitive.

---

Repo: `<repo root>`. Read `AGENTS.md` and `.agents/skills/sdf-verification/SKILL.md` FIRST and follow them.
Goal: add ONE genuinely novel, exact (or provably conservative) SDF primitive/operator, fully verified and documented.

## Hard rules
- Only CREATE your own new files:
  - `glsl/<category>/<name>.glsl`: header comment block, `#include "../common/common.glsl"` if you use helpers,
    a unique `#ifndef SDF_GEMS_HAS_<NAME>` guard, and an `@author` line naming the model you are.
  - `docs/<category>/<name>.md`: title, images, Category / GLSL File / Attribution / Exactness bullets,
    Overview & Motivation, Mathematical Derivation, GLSL code, Verification table with measured numbers,
    before/after comparison, and a Novelty & Prior Art section.
  - `docs/images/<category>/<name>.png` (3D: 360x270 `render_3d_sdf`; 2D: 320x320 `render_2d_sdf`),
    `<name>_slice.png` (360x360) for 3D, and `docs/images/comparisons/<name>_comparison.png`.
- Do NOT edit `README.md`, `sdf_gems.glsl`, or any existing file. Do NOT run git add/commit/push.
- Never put personal names or usernames in any file, including repository URLs that contain a username.
- **Attribution**: credit the model YOU are, using its human-readable name. If you are not certain of your exact
  model name or version, ask the parent agent before writing it. Do not guess and do not copy another model's name.
- Run Python only via uv: `cd <repo root> && uv run python <script>`.
- Keep scratch scripts in `<scratch dir>/<your name>/`, never in the repo.
- Choose unique GLSL names; guard shared helpers with `#ifndef SDF_GEMS_HAS_<NAME>`.

## Tooling (all in `tools/`)
- `tools/gpu_eval.py`: `gpu_eval`, `gpu_fn`, `glsl_args`, `compile_umbrella(new=[...], call=...)`. Tests the LITERAL GLSL on
  the GPU (FP32). Format numbers with `glsl_args(...)`: interpolating numpy scalars yields `np.float64(...)` and breaks GLSL.
- `tools/ground_truth.py`: `polygon_sdf` / `polyline_distance` (2D curves), `convex_hull_sdf` (polyhedra, or hulls of dense
  curve samples), `implicit_boundary_cloud` + `cloud_sdf` (shapes defined by an inside predicate).
- `tools/sdf_checks.py`: `probe_points`, `check_finite`, `compare_gt`, `seam_pair_test`, `eikonal_fd`, `check_sdg`, `self_projection`,
  `lipschitz_pairs`, `surface_zero`, `exactness_certificate`.
- `tools/render_sdf.py`: `render_3d_sdf`, `render_slice_3d`, `render_2d_sdf`, `render_error_heatmap` (preferred before/after view),
  `raymarch_stats` (steps/pixel + overstep pixels), `render_eikonal_error_map`, `combine_panels`. Render from the literal GLSL via `gpu_fn`.
- View EVERY PNG with `read_file` and critique it (framing, clipping, artifacts, isolines, and whether the before/after
  makes the improvement obvious). Re-render until good.

## Method
1. Brute-force-test every lemma (candidate feature sets, "nearest point lies in X", unimodality, ...) on 10^5–10^6 points
   BEFORE writing GLSL. Record where it fails, not only where it holds.
2. Build an independent FP64 ground truth from the shape *definition* only (never from the formula under test).
3. Required checks (report numbers): NaN/Inf = 0 at probes and degenerate parameters; max |d - gt| on >= 20k points
   (uniform + near-surface + near sharp features); seam C0 on points placed ON each branch boundary; eikonal; `sdg` checks;
   self-projection (including points exactly on the surface); `exactness_certificate` with surface samples generated
   independently from the shape definition; a regression against an existing verified primitive in any special case
   (e.g. helix angle 0); and checks of the shape definition's own invariants (width, volume, chirality, symmetry).
   If a few points disagree with the reference, refine the reference before changing the shader.
4. Eikonal maps alone are NOT sufficient: max()/min() of exact distances has |grad| = 1 yet can be wrong. Use error vs ground truth.
5. FP32 pitfalls: avoid `sign(x + eps)` and `v / max(len, eps)` in gradients (they can return zero vectors on the surface). Use
   ternaries with unit-vector fallbacks. Write radicands as products, `(r - a) * (r + a)`, never `r*r - a*a`. Make branch
   selection rounding-tolerant so surface points cannot fall between an early exit and a region test. Watch for
   cancellation; cluster bracketing/sample nodes where roots bunch up.
6. Novelty: several `search_web` queries (Shadertoy, iq articles, hg_sdf, GitHub, papers). Report prior art honestly and
   state exactly what is new.

## Final report to parent
Files created; the novel idea in one paragraph; measured metrics; prior-art findings; the model name you used for
attribution; a proposed README row (3 columns); a proposed commit message `feat(<category>): ...`.
