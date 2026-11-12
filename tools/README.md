# `sdf-gems` Developer Tooling (`tools/`)

This directory contains the repeatable Python scripts used to render 3D/2D SDF previews, visualize Eikonal gradient errors, regenerate the [`sdf_gems.glsl`](../sdf_gems.glsl) umbrella header, and audit documentation links.

## 1. Regenerate [`sdf_gems.glsl`](../sdf_gems.glsl) & Audit Docs / Images

```bash
uv run tools/update_umbrella.py
```

- Automatically writes [`sdf_gems.glsl`](../sdf_gems.glsl) with `#include` directives for all `.glsl` files in [`glsl/`](../glsl/) (`glsl/common/common.glsl` first, then alphabetical).
- Verifies that every `.glsl` primitive/operator has a matching Markdown article in [`docs/`](../docs/).
- Verifies that every PNG image linked in [`README.md`](../README.md) and [`docs/**/*.md`](../docs/) exists on disk.

## 2. Render 3D Surfaces, 2D Slices, & Comparison Heatmaps ([`tools/render_sdf.py`](render_sdf.py))

Import [`tools/render_sdf.py`](render_sdf.py) from any script or one-liner via `uv run`:

```bash
uv run python3 -c "
from tools.render_sdf import render_3d_sdf, render_slice_3d, render_2d_sdf
import numpy as np

def sd_jack_3d(p, w=0.80, r=0.14):
    q = np.abs(p)
    s = np.minimum(q[..., 0] + q[..., 1] + q[..., 2], w) / 3.0
    return np.linalg.norm(q - s[..., None], axis=-1) - r

img_3d = render_3d_sdf(sd_jack_3d, res=(360, 270))
img_slice = render_slice_3d(sd_jack_3d, res=320, label='sdJack3D')
"
```

## 3. Test the Literal GLSL on the GPU ([`tools/gpu_eval.py`](gpu_eval.py))

```bash
uv run tools/gpu_eval.py                                                    # 1-TU compile, forward + reverse order
uv run tools/gpu_eval.py --new glsl/primitives3d/sdFoo.glsl --call "vec4(sdFoo(P, 1.0))"   # + new-file-first
```

```python
from tools.gpu_eval import gpu_eval, gpu_fn, glsl_args
sd  = gpu_fn(["glsl/polyhedra/sdSnubCube.glsl"], f"vec4(sdSnubCube(p, {glsl_args(0.55)}))")  # numpy callable
out = gpu_eval(["glsl/polyhedra/sdSnubCube.glsl"], "sdgSnubCube(p, 0.55)", pts)             # (N,4) float32
```

- Compiles with the platform's real GLSL compiler (OpenGL 4.1 core via `moderngl`) and evaluates in FP32 via transform feedback.
- Always format numbers with `glsl_args(...)`: interpolating numpy scalars directly yields `np.float64(...)`, which does not compile.

## 4. Independent Ground Truth ([`tools/ground_truth.py`](ground_truth.py))

Builders that use only the shape *definition*, never the formula under test:
- `polygon_sdf` / `polyline_distance`: dense outlines of 2D shapes and open curves.
- `convex_hull_sdf` / `point_triangle_distance`: polyhedra, or convex hulls of dense curve samples (e.g. an oloid as the hull of two sampled circles).
- `implicit_boundary_cloud` + `cloud_sdf`: shapes given by an inside predicate (ray-cast bisection, KD-tree, tangent-plane refinement). Add exact samples of sharp edges and vertices to the cloud.

`uv run tools/ground_truth.py` self-tests the builders against analytic sphere, disk and box distances.

## 5. Numerical Checks ([`tools/sdf_checks.py`](sdf_checks.py))

`probe_points`, `check_finite`, `compare_gt`, `seam_pair_test`, `eikonal_fd`, `check_sdg`, `self_projection`, `lipschitz_pairs`, `surface_zero`, `exactness_certificate`.
`self_projection` (step `x = p - d(p) grad(p)` and require `d(x) = 0`) needs no ground truth and catches wrong distances and gradients, including FP32 rounding-shell bugs on the surface itself.
`exactness_certificate` combines 1-Lipschitz pairs, zero on independently generated surface samples, and self-projection: together they imply the SDF is exact, with no ground-truth distances.

## 6. Error Heatmaps (`render_error_heatmap` in [`tools/render_sdf.py`](render_sdf.py))

`render_error_heatmap(fn, gt, ...)` draws log10 |fn - gt| on a slice with a colorbar and the true surface in cyan. Prefer it to eikonal maps for before/after comparisons: a `max()` of exact distances has unit gradient yet can be badly wrong.

## 7. Raymarching Cost & Safety (`raymarch_stats` in [`tools/render_sdf.py`](render_sdf.py))

`stats, img = raymarch_stats(fn)` sphere-traces like `render_3d_sdf` and reports mean/95th-percentile steps per pixel, pixels needing more than 18 steps, and **overstep pixels** (rays ending inside the surface, drawn in red). Use it in before/after comparisons of bounds (slower marching) and overestimating SDFs (unsafe).

## 8. GitHub-Safe Math in Markdown ([`tools/lint_math.py`](lint_math.py))

```bash
uv run tools/lint_math.py                       # fail on bare $...$ / $$...$$ or `|` in table-row math/code
uv run tools/lint_math.py --fix                 # convert to $`...`$ (inline) and ```math (display)
uv run tools/lint_math.py --render              # render every formula with MathJax 3 + KaTeX (Node.js)
uv run tools/lint_math.py --github              # verify against the GitHub Markdown API (60 req/h; GITHUB_TOKEN)
```

GitHub processes Markdown before MathJax: backslash escapes before punctuation are consumed (`\%` becomes a TeX comment, `\\` collapses matrix rows, `\{` `\|` `\,` vanish) and emphasis/link syntax runs across formulas. GitHub's Markdown-proof syntax (`` $`...`$ `` and ```` ```math ````) passes TeX through verbatim. The offline lint also runs inside `tools/update_umbrella.py`.
