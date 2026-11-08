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
| `tools/sdf_checks.py` | `probe_points`, `check_finite`, `compare_gt`, `seam_pair_test`, `eikonal_fd`, `check_sdg`, `self_projection`. |
| `tools/render_sdf.py` | `render_3d_sdf`, `render_slice_3d`, `render_2d_sdf`, `render_error_heatmap`, `render_eikonal_error_map`, `combine_panels`. |

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
   - **Self-Projection** (`self_projection`): $x = p - d(p)\,\nabla d(p)$ must land on the surface ($|d(x)| \approx 0$). This needs no ground truth and catches wrong distances and wrong gradients at once.
3. **FP32 Pitfalls Seen in Practice**:
   - Format GLSL literals with `glsl_args(...)`: interpolating numpy scalars produces `np.float64(...)`, which does not compile.
   - Avoid `sign(x + 1e-30)` and `v / max(len, eps)` in gradients: they can yield zero vectors exactly on the surface. Use ternaries with unit-vector fallbacks.
   - Fixed-iteration solvers: cluster bracketing nodes (e.g. Chebyshev–Lobatto) where roots bunch up near symmetry planes; uniform nodes can lose two orders of magnitude of accuracy there.

## 3. Rendering & Mandatory Visual Inspection
Use `tools/render_sdf.py` via `uv run` to generate documentation visuals, preferably from the literal GLSL through `gpu_fn`:
1. **3D Preview (`docs/images/<category>/<name>.png`)**:
   - Render with `render_3d_sdf(fn, res=(360, 270))` showing the 3D zero-isosurface and ground-plane distance isolines.
   - Ensure the shape is centered, occupies 40%–70% of the frame, and rests cleanly above `plane_y`.
2. **2D Cross-Section Slice (`docs/images/<category>/<name>_slice.png`)**:
   - Render with `render_slice_3d` or `render_2d_sdf` at `res=(360, 360)`.
3. **Before/After Comparison (`docs/images/comparisons/`)**:
   - Prefer `render_error_heatmap(fn, gt)` (log10 error vs ground truth, with a colorbar) over eikonal maps; add 3D renders of offsets (`d - r`) or shells where the old approach visibly fails.
4. **Inspect with `read_file`**:
   - Call `read_file` on **every** generated `.png` before staging or committing. Reject and re-render any image that is clipped, blank, overly zoomed, or uses the wrong aspect ratio, or whose before/after does not make the improvement obvious.

## 4. Parallel Discovery with Subagents
For batches of independent primitives, start from [`subagent_brief.md`](subagent_brief.md). Subagents only create new files. The parent agent re-verifies each result independently, adds the README row and the `sdf_gems.glsl` include (relative to `HEAD`, so in-progress files from other subagents are not swept in), runs the compile check, views every image, and makes one commit per primitive. Subagents attribute work to the model they are, and ask the parent when unsure of their model name.
