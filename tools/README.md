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
