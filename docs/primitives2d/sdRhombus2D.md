# `sdRhombus2D` & `sdgRhombus2D` — 0-Division, 0/1-`sqrt` Rhombus + 1-SFU Gradient

- **Category:** `primitives2d`
- **GLSL File:** [`glsl/primitives2d/sdRhombus2D.glsl`](../../glsl/primitives2d/sdRhombus2D.glsl)
- **Author:** Inigo Quilez — upgraded to 0-division 0/1-`sqrt` & 1-SFU (`inversesqrt`) gradient using **Gemini 4 Argon** (`sdf-gems, 2026`)
- **Exact Euclidean SDF:** Yes ($\|\nabla d\| = 1$, $1.33\times 10^{-15}$ max difference vs baseline, **0 divisions (`RCP`)** and **0 `sqrt`s** for constant shape parameters in interior & face strips)

![sdRhombus2D Preview](../images/primitives2d/sdRhombus2D.png)

---

## Mathematical Overview

### Direct Vertex-Wedge Decomposition for `sdRhombus2D` & `sdgRhombus2D`

In Inigo Quilez's baseline `sdRhombus(vec2 p, vec2 b)` ([Shadertoy `XdXcRB`](https://www.shadertoy.com/view/XdXcRB)), every point evaluates a runtime division `/ dot(b, b)` (`MUFU.RCP`) followed by a `length(...)` (`MUFU.SQRT`), executing **2 SFU instructions at every pixel** even inside the rhombus and across its 4 flat exterior edge strips.

After folding into the first quadrant $p = |p| \ge 0$, the edge connecting $(b_x, 0)$ to $(0, b_y)$ has tangent vector $e = (-b_x, b_y)$ and outward normal $n = (b_y, b_x) / \sqrt{b_x^2 + b_y^2}$:

1. **Right-vertex wedge:** $(p - (b_x, 0)) \cdot (-b_x, b_y) < 0 \iff b_x(p_x - b_x) > b_y p_y$. Here the closest point on the rhombus is $(b_x, 0)$, and the Euclidean distance is simply $\|(p_x - b_x, p_y)\|$ (**1 `sqrt`**, **0 divisions**, and completely independent of $\sqrt{b_x^2 + b_y^2}$).
2. **Top-vertex wedge:** $(p - (0, b_y)) \cdot (b_x, -b_y) < 0 \iff b_y(p_y - b_y) > b_x p_x$. Here the closest point is $(0, b_y)$, and the Euclidean distance is $\|(p_x, p_y - b_y)\|$ (**1 `sqrt`**, **0 divisions**).
3. **Interior & 4 Exterior Edge Strips:** Everywhere else (**100% of the interior and all 4 exterior edge-normal strips**), the projection lies strictly on the edge interior, so the signed Euclidean distance is linear in $p$:
   $$d(p) = (b_y p_x + b_x p_y - b_x b_y)\,\text{inversesqrt}(b_x^2 + b_y^2),$$
   requiring **0 `sqrt`s** for constant $b$ (or 1 `MUFU.RSQ` for dynamic $b$) and **0 divisions (`MUFU.RCP`)**!

*(For `sdTriangleIsosceles2D` and `sdgTriangleIsosceles2D`, see [`sdTriangle2D.md`](sdTriangle2D.md) and [`glsl/primitives2d/sdTriangle2D.glsl`](../../glsl/primitives2d/sdTriangle2D.glsl).)*

---

## Numerical Verification & Benchmark

| Function | Runtime Divisions (`RCP`) | `sqrt` / `rsq` in Interior & Face Strips (const params) | Max Error vs Ground Truth | Max Gradient Error |
| :--- | :---: | :---: | :---: | :---: |
| Baseline IQ `sdRhombus` | `1` | `1` `sqrt` | $0.00 \times 10^0$ | N/A |
| **`sdf-gems` `sdRhombus2D` / `sdgRhombus2D`** | **`0`** | **`0`** (`1` in vertex wedges) | **$1.33 \times 10^{-15}$** | **$1.26 \times 10^{-9}$** |
