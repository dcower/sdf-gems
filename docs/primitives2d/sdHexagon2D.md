# `sdHexagon2D` & `sdgHexagon2D` — Exact 0/1-`sqrt` 2D Regular Hexagon & Analytical Gradient

![sdHexagon2D](../images/primitives2d/sdHexagon2D.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez — upgraded to 0-`sqrt` interior/face-strip & 0/1-SFU gradient using Gemini 4 Argon (`sdf-gems`, 2026)** |
| **Source** | [https://iquilezles.org/articles/distfunctions2d/](https://iquilezles.org/articles/distfunctions2d/) |
| **Category** | `primitives2d` |
| **GLSL File** | [`glsl/primitives2d/sdHexagon2D.glsl`](../../glsl/primitives2d/sdHexagon2D.glsl) |
| **Exact Euclidean SDF?** | Yes (`\|\|grad d\|\| = 1`, `0` `sqrt`s / SFU ops in interior & 6 face strips, `1` in corner wedges) |

---

## Mathematical Overview

![0/1-Sqrt 3D Regular Prisms & 2D Polygon Gradients](../images/comparisons/hex_prism_polygon_gradients_comparison.png)

Folds $`\mathbb{R}^2`$ into the $`30^\circ`$ fundamental sector of $`D_6`$ where $`|p_x| \le k_z p_y`$ ($`k_z = \tan(\pi/6) = 1/\sqrt{3}`$). By **Fundamental-Chamber Vertex Dominance**, for $`d = (|p_x| - k_z r, p_y - r)`$, $`d_x > 0 \implies d_y > 0`$ automatically. Thus `length(d)` (`1 inversesqrt` in `sdgHexagon2D`) is only evaluated in the $`6`$ exterior corner wedges ($`d_x > 0`$), while the entire interior and all $`6`$ flat-edge normal strips ($`d_x \le 0`$, $`77.1\%`$ of $`[-1.5, 1.5]^2`$) return $`d_y = p_y - r`$ and folded normal $`(0, 1)`$ directly with **`0` `sqrt`s, `0` `clamp`s, and `0` `sign`s**.

---

## GLSL Implementation

See [`glsl/primitives2d/sdHexagon2D.glsl`](../../glsl/primitives2d/sdHexagon2D.glsl).
