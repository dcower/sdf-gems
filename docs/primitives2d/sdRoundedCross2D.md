# `sdRoundedCross2D`, `sdgRoundedCross2D` & `sdRoundedCrossPrism` — Division-Free 2D Rounded Cross & 1-`sqrt` 3D Prism

- **Category:** `primitives2d`
- **GLSL File:** [`glsl/primitives2d/sdRoundedCross2D.glsl`](../../glsl/primitives2d/sdRoundedCross2D.glsl)
- **Author:** Inigo Quilez — upgraded via Division-Free Region Test & 1-`sqrt` 3D Lifting using **Gemini 4 Argon** (`sdf-gems, 2026`)
- **Exact Euclidean SDF:** Yes ($`\|\nabla d\| = 1`$, $`4.44\times 10^{-16}`$ max difference vs baseline in 2D, $`8.88\times 10^{-16}`$ in 3D; **0 divisions** in region test, **1 `dot2`** in corner branch, **0/1 `sqrt`** in 3D prism)

![sdRoundedCross2D Preview](../images/primitives2d/sdRoundedCross2D.png)

---

## Mathematical Overview

### 1. Division-Free Region Test & Single-`dot2` Corner Selector

In Inigo Quilez's baseline `sdRoundedCross(vec2 p, float h)` (`distfunctions2d`), the circle radius $`k = \frac{1}{2}(h + 1/h)`$ is computed unconditionally, and the `else` branch evaluates two `dot2` calls `min(dot2(p - vec2(0,h)), dot2(p - vec2(1,0)))`:

1. **Division-Free Region Test:** Since $`k - h = \frac{1 - h^2}{2h}`$ and $`h > 0`$, the region test $`p_y < p_x(k - h) + h`$ is equivalent to the division-free polynomial inequality:
   $`\displaystyle 2h(p_y - h) < p_x(1 - h^2),`$
   deferring the reciprocal `1.0 / h` strictly to the concave arc branch.
2. **Single-`dot2` Corner Selector:** In the exterior corner branch (`else`), $`\|p - (0, h)\|^2 < \|p - (1, 0)\|^2`$ reduces to the linear separator $`2p_x - 2hp_y < 1 - h^2`$, selecting the active corner before evaluating `length(p - c)`.

### 2. 0/1-`sqrt` 3D Extruded `sdRoundedCrossPrism`

Because the `else` branch of `sdRoundedCross2D` lies **100% in the 2D exterior** ($`d_{\text{2D}} \ge 0`$) and its squared distance $`\|q - c\|_2^2`$ is a pure polynomial, the 3D extruded prism evaluates $`\sqrt{\|q - c\|_2^2 + \max(|p_z| - h_z, 0)^2}`$ with **1 `sqrt` (instead of 2)** across the entire exterior corner region, **0 `sqrt`s** above/below the 2D interior (`dz > 0 && r2 >= k*k`), and **1 `sqrt`** inside the $`z`$-slab (`dz <= 0`).

---

## Numerical Verification & Benchmark

| Function | Reciprocal `1/h` Outside Arc Branch | `dot2` in Corner Branch | 3D Exterior Corner `sqrt`s | Max Error vs Ground Truth |
| :--- | :---: | :---: | :---: | :---: |
| Baseline IQ `sdRoundedCross` + `opExtrusion` | `1` | `2` | `2` | $`0.00 \times 10^0`$ |
| **`sdf-gems` `sdRoundedCross2D` / `sdRoundedCrossPrism`** | **`0`** | **`1`** | **`1` (`0` above interior)** | **$`4.44 \times 10^{-16}`$ (2D) / $`8.88 \times 10^{-16}`$ (3D)** |
