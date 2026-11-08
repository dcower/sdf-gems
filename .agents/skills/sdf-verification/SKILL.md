---
name: sdf-verification
description: >-
  Workflow for verifying, benchmarking, and rendering 2D/3D Signed Distance Functions (SDFs)
  and analytical gradients in sdf-gems. Activate this skill whenever adding, optimizing,
  or auditing .glsl files or rendering documentation images in docs/images/.
---

# SDF Verification, Benchmarking & Rendering Workflow

Follow this checklist whenever creating, upgrading, or auditing any SDF or operator in `glsl/`.

## 0. Tooling (use these instead of writing new harnesses)
| Tool | Purpose |
| :--- | :--- |
| `tools/gpu_eval.py` | `gpu_eval` / `gpu_fn` run the **literal** `.glsl` on the GPU (OpenGL 4.1, FP32). `glsl_args(...)` formats numbers safely. `uv run tools/gpu_eval.py --new <file> --call "<vec4 expr of P>"` performs the single-translation-unit compile in forward, reverse and new-file-first order. |
| `tools/ground_truth.py` | Independent FP64 ground truth: `polygon_sdf`, `polyline_distance`, `convex_hull_sdf`, `point_triangle_distance`, `implicit_boundary_cloud` + `cloud_sdf`. |
| `tools/sdf_checks.py` | `probe_points`, `check_finite`, `compare_gt`, `seam_pair_test`, `eikonal_fd`, `check_sdg`, `self_projection`, `lipschitz_pairs`, `surface_zero`, `exactness_certificate`. |
| `tools/lint_math.py` | GitHub-safe math in Markdown: `--fix` converts to `` $`...`$ `` / ```` ```math ````, `--render` renders every formula with MathJax 3 + KaTeX, `--github` verifies against GitHub's own renderer (60 requests/hour unauthenticated). Also enforced by `tools/update_umbrella.py`. |
| `tools/render_sdf.py` | `render_3d_sdf`, `render_slice_3d`, `render_2d_sdf`, `render_error_heatmap`, `render_eikonal_error_map`, `raymarch_stats`, `combine_panels`. |

Wrap the literal shader with `gpu_fn(["glsl/<cat>/<name>.glsl"], "vec4(sdFoo(p, 1.0))")` and pass that callable to the
checks and renderers, so every number and image comes from the code being committed.

## 1. Preserve Geometric Intent & Evaluate True GPU Cost
1. **Verify Geometric Definition First**: Before optimizing an existing SDF or operator (e.g., `opUnionStairs`, `sdSteinmetz`), plot its zero-level set and cross-section against the reference implementation to confirm the exact intended geometry (sharp steps vs. smooth ripples, union vs. intersection).
2. **Count Total SFU Operations (`sqrt` + `rcp`/`div`), Not Just `sqrt`**:
   - On modern GPU architectures, `sqrt`, `inversesqrt`, and floating-point division (`1.0 / x`) all execute on the Special Function Unit (SFU) at roughly ~4x the latency of an FMA.
   - Never trade 1 `sqrt` for multiple floating-point divisions or deep divergent branching. Favor polynomial/dot-product Voronoi early exits (`0 SFU`) and fused `inversesqrt` normalizations.
3. **Brute-Force Every Lemma Before Coding It**: candidate feature sets ("only these faces/edges can be nearest"), "nearest point lies in the star of the nearest vertex", unimodality of a 1D objective, and so on. Test on 10^5–10^6 points first, and record where each lemma *fails* (e.g. exterior-only lemmas).

## 2. Direct-GLSL Numerical Audit (Never Trust Python Prototypes Alone)
Python prototypes can mask transcription bugs in constants, variable swaps, or integer/float promotion rules when ported to `.glsl`. Always test the literal `.glsl` file, using `tools/gpu_eval.py` on the GPU (FP32, the real target). An FP64 port is useful in addition to measure algorithmic error separately from FP32 rounding.

1. **Single-Translation-Unit Compilation Check**:
   - `uv run tools/gpu_eval.py --new <your files> --call "<expression calling your functions>"` compiles all `glsl/*/*.glsl` together in forward, reverse and new-file-first order with the real GLSL compiler.
2. **Mandatory Numerical Checks on Every Function** (`tools/sdf_checks.py`):
   - **Zero `NaN`/`Inf`s**: `probe_points` covers the origin, the axes, symmetry planes and diagonals; also test degenerate parameters (e.g., equal radii $r_1 = r_2$, zero bulge $b \to 0$, scale $\to 0$).
   - **Independent Brute-Force Ground Truth**: build it from the shape *definition* only (`tools/ground_truth.py`), never from the formula under test. Sample $N \ge 20,000$ points (uniform + near-surface + near sharp features) and verify that $|d(p) - d_{\text{exact}}(p)| < 10^{-5}$ (for exact SDFs; GPU FP32 typically reaches `~1e-6`) and $d(p) \le d_{\text{exact}}(p) + 10^{-5}$ (for conservative bounds). When integrating someone else's primitive, re-verify with a *second*, differently constructed ground truth.
   - **$C^0$ Branch Continuity**: `seam_pair_test` on points placed **on** each branch boundary (Voronoi seams, fold planes, solver bracket boundaries), not only random points: $|d(p+\epsilon u) - d(p-\epsilon u)| \le 2\epsilon + 10^{-4}$.
   - **Eikonal Accuracy**: `eikonal_fd`. **Caution**: an eikonal map cannot reveal every wrong SDF. A `max()`/`min()` of exact distances has $\|\nabla d\| = 1$ almost everywhere and can still be badly wrong, so always compare values against ground truth too.
   - **Analytical Gradient Verification (`sdg*`)** via `check_sdg`:
     - Verify `dg.x == sd*(p)` (`< 1e-6`).
     - Verify `length(dg.yzw) == 1.0` (`< 1e-4`), including finite unit-length fallbacks on medial axes and at the origin.
     - Verify `dg.yzw` matches central finite differences of `dg.x` (`> 0.999` cosine similarity away from medial seams).
   - **Self-Projection** (`self_projection`): $x = p - d(p)\,\nabla d(p)$ must land on the surface ($|d(x)| \approx 0$). This needs no ground truth and catches wrong distances and wrong gradients at once. Include points placed *exactly on* the surface: FP32 rounding shells are where branches disagree.
   - **Exactness Certificate** (`exactness_certificate`): (1) 1-Lipschitz on random pairs, (2) $d = 0$ on surface samples generated independently from the shape definition, and (3) self-projection lands on the zero set. Together these imply $d$ is the exact SDF, with no ground-truth distance computation. Run it on every exact primitive, in addition to the brute-force comparison.
   - **Special-Case Regression**: when a new primitive reduces to an existing verified one for some parameter (helix angle $0$, equal radii, $n \to$ a known polygon), compare against it numerically.
   - **Validate the Shape Definition Itself**: before trusting a ground truth built from a definition, check its invariants (e.g. constant width in thousands of directions, known volume, chirality vs. the mirror image, symmetry orbit invariance).
   - **Before Blaming the Shader**: when a handful of points disagree with the reference, refine the reference first (denser samples, local search, a second method). Several "outliers" in this repo's history were reference misses, not shader bugs.
3. **FP32 Pitfalls Seen in Practice**:
   - Format GLSL literals with `glsl_args(...)`: interpolating numpy scalars produces `np.float64(...)`, which does not compile.
   - Avoid `sign(x + 1e-30)` and `v / max(len, eps)` in gradients: they can yield zero vectors exactly on the surface. Use ternaries with unit-vector fallbacks.
   - Fixed-iteration solvers: cluster bracketing nodes (e.g. Chebyshev–Lobatto) where roots bunch up near symmetry planes; uniform nodes can lose two orders of magnitude of accuracy there.
   - Write radicands as products, not differences: `sqrt(max((r - a) * (r + a), 0.0))` instead of `sqrt(r*r - a*a)`. With constant-folded parameters and FMA contraction, `sqrt(rb*rb - rb*rb)` can become `sqrt(-eps)` = NaN.
   - Make branch selection rounding-tolerant: if an early-exit test (`dm <= 0`) and a later region test (`h > 0`) are computed along different arithmetic paths, a surface point can satisfy neither. Accept candidates within a tiny tolerance or provably on the far side of the boundary (each candidate is still a true distance to a boundary point, so the result stays an upper bound and exact up to rounding).
   - A finite-difference gradient in FP32 makes self-projection noise-dominated (noise × $d$). Use the analytical `sdg` or an FP64 port for self-projection.

## 3. Rendering & Mandatory Visual Inspection
Use `tools/render_sdf.py` via `uv run` to generate documentation visuals, preferably from the literal GLSL through `gpu_fn`:
1. **3D Preview (`docs/images/<category>/<name>.png`)**:
   - Render with `render_3d_sdf(fn, res=(360, 270))` showing the 3D zero-isosurface and ground-plane distance isolines.
   - Ensure the shape is centered, occupies 40%–70% of the frame, and rests cleanly above `plane_y`.
2. **2D Cross-Section Slice (`docs/images/<category>/<name>_slice.png`)**:
   - Render with `render_slice_3d` or `render_2d_sdf` at `res=(360, 360)`.
3. **Before/After Comparison (`docs/images/comparisons/`)**:
   - Prefer `render_error_heatmap(fn, gt)` (log10 error vs ground truth, with a colorbar) over eikonal maps; add 3D renders of offsets (`d - r`) or shells where the old approach visibly fails.
   - For bounds and overestimating SDFs, add `raymarch_stats(fn)`: mean steps/pixel and overstep pixels (in red), the cost and safety a GPU raymarcher actually sees.
4. **Inspect with `read_file`**:
   - Call `read_file` on **every** generated `.png` before staging or committing. Reject and re-render any image that is clipped, blank, overly zoomed, or uses the wrong aspect ratio, or whose before/after does not make the improvement obvious.

## 4. Documentation Math (GitHub)
- Inline math is `` $`...`$ ``, display math is a ```` ```math ```` block; inside a list item use `` $`\displaystyle ...`$ `` on its own line (GitHub skips ```math blocks nested in lists when the page has inline math).
- Never write bare `$...`` $` / `$ ``$...$$`: GitHub consumes backslash escapes before MathJax (`\%` becomes a comment that deletes the rest of the formula, `\\` collapses matrix rows, `\{`/`\|`/`\,` vanish) and pairs `*` / `[...]` across formulas.
- In table rows, no unescaped `|` in math (use `\Vert`, `\vert`) or in code spans (write `\|`).
- Before committing docs: `uv run tools/lint_math.py --render` (offline, MathJax + KaTeX); for large doc changes also `--github`.

## 5. Parallel Discovery with Subagents
For batches of independent primitives, start from [`subagent_brief.md`](subagent_brief.md). Subagents only create new files. The parent agent re-verifies each result independently, adds the README row and the `sdf_gems.glsl` include (relative to `HEAD`, so in-progress files from other subagents are not swept in), runs the compile check, views every image, and makes one commit per primitive. Subagents attribute work to the model they are, and ask the parent when unsure of their model name.
