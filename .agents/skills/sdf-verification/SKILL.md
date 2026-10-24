---
name: sdf-verification
description: >-
  Workflow for verifying, benchmarking, and rendering 2D/3D Signed Distance Functions (SDFs)
  and analytical gradients in sdf-gems. Activate this skill whenever adding, optimizing,
  or auditing .glsl files or rendering documentation images in docs/images/.
---

# SDF Verification, Benchmarking & Rendering Workflow

Follow this checklist whenever creating, upgrading, or auditing any SDF or operator in `glsl/`.

## 1. Preserve Geometric Intent & Evaluate True GPU Cost
1. **Verify Geometric Definition First**: Before optimizing an existing SDF or operator (e.g., `opUnionStairs`, `sdSteinmetz`), plot its zero-level set and cross-section against the reference implementation to confirm the exact intended geometry (sharp steps vs. smooth ripples, union vs. intersection).
2. **Count Total SFU Operations (`sqrt` + `rcp`/`div`), Not Just `sqrt`**:
   - On modern GPU architectures, `sqrt`, `inversesqrt`, and floating-point division (`1.0 / x`) all execute on the Special Function Unit (SFU) at roughly ~4x the latency of an FMA.
   - Never trade 1 `sqrt` for multiple floating-point divisions or deep divergent branching. Favor polynomial/dot-product Voronoi early exits (`0 SFU`) and fused `inversesqrt` normalizations.

## 2. Direct-GLSL Numerical Audit (Never Trust Python Prototypes Alone)
Python prototypes can mask transcription bugs in constants, variable swaps, or integer/float promotion rules when ported to `.glsl`. Always compile and test the literal `.glsl` file in C++17:

1. **Single-Translation-Unit Compilation Check**:
   - Verify that all files in `glsl/*/*.glsl` compile together in a single namespace (both forward and reverse inclusion order) with `0` duplicate symbols.
2. **Five Mandatory Numerical Checks on Every Function**:
   - **Zero `NaN`/`Inf`s**: Test at `p = (0,0,0)`, along all coordinate axes `(x,0,0)`, `(0,y,0)`, `(0,0,z)`, on symmetry planes (`x = y`, `x = y = z`), and for degenerate parameters (e.g., equal radii $r_1 = r_2$, zero bulge $b \to 0$).
   - **Brute-Force Surface Ground Truth**: Sample $N \ge 20,000$ surface points in FP64 and verify that $|d(p) - d_{\text{exact}}(p)| < 10^{-5}$ (for exact SDFs) and $d(p) \le d_{\text{exact}}(p) + 10^{-5}$ (for conservative bounds).
   - **$C^0$ Branch Continuity**: Probe pairs of points $p \pm \epsilon u$ ($\epsilon = 10^{-5}$) straddling Voronoi region boundaries and verify $|d(p+\epsilon u) - d(p-\epsilon u)| \le 2\epsilon + 10^{-4}$.
   - **Eikonal Accuracy**: Verify $\|\nabla d(p)\| \approx 1$ via central finite differences across both interior and exterior regions.
   - **Analytical Gradient Verification (`sdg*`)**:
     - Verify `dg.x == sd*(p)` (`< 1e-6`).
     - Verify `length(dg.yzw) == 1.0` (`< 1e-4`), including finite unit-length fallbacks on medial axes and at the origin.
     - Verify `dg.yzw` matches central finite differences of `dg.x` (`> 0.999` cosine similarity away from medial seams).

## 3. Rendering & Mandatory Visual Inspection
Use `tools/render_sdf.py` via `uv run` to generate documentation visuals:
1. **3D Preview (`docs/images/<category>/<name>.png`)**:
   - Render with `render_3d_sdf(fn, res=(360, 270))` showing the 3D zero-isosurface and ground-plane distance isolines.
   - Ensure the shape is centered, occupies 40%–70% of the frame, and rests cleanly above `plane_y`.
2. **2D Cross-Section Slice (`docs/images/<category>/<name>_slice.png`)**:
   - Render with `render_slice_3d` or `render_2d_sdf` at `res=(360, 360)`.
3. **Inspect with `read_file`**:
   - Call `read_file` on **every** generated `.png` before staging or committing. Reject and re-render any image that is clipped, blank, overly zoomed, or uses the wrong aspect ratio.
