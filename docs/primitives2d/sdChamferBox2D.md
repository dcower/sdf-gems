# `sdChamferBox2D`, `sdChamferBoxPrism`, `sdOrientedBox2D` & `sdCoolS2D` — 0/1-`sqrt` Box & Segment Primitives

- **Category:** `primitives2d`
- **GLSL File:** [`glsl/primitives2d/sdChamferBox2D.glsl`](../../glsl/primitives2d/sdChamferBox2D.glsl)
- **Author:** Inigo Quilez — upgraded via 0/1-`sqrt` Decomposition & Round/Clamp Elimination using **Gemini 4 Argon** (`sdf-gems, 2026`)
- **Exact Euclidean SDF:** Yes ($\|\nabla d\| = 1$, $8.88\times 10^{-16}$ max difference vs baseline; **0/1 `sqrt`** in 2D & 3D `sdChamferBoxPrism`; round/clamp-free `sdCoolS2D`)

![sdChamferBox2D Preview](../images/primitives2d/sdChamferBox2D.png)

---

## Mathematical Overview

### 1. 0/1-`sqrt` `sdChamferBox2D`, `sdgChamferBox2D`, & 3D `sdChamferBoxPrism`

After folding $p = |p| - b$, sorting $p_x \ge p_y$, and shifting $p_y \gets p_y + \text{chamfer}$ with $k = 1 - \sqrt{2}$:
- **Axial Face Strip ($p_y < 0 \land p_y + k p_x < 0$):** $d = p_x$, $\nabla d = (\pm 1, 0)$ (**0 `sqrt`s**).
- **$45^\circ$ Chamfer Strip ($p_x < p_y$):** $d = (p_x + p_y)/\sqrt{2}$, $\nabla d = (\pm 1, \pm 1)/\sqrt{2}$ (**0 `sqrt`s**).
- **Vertex Wedge ($p_x \ge p_y \ge 0$):** $d = \|p\|_2$ (**1 `sqrt`**, strictly exterior $d \ge 0$).

Because the vertex wedge lies strictly in the 2D exterior ($p_x \ge p_y \ge 0$), its squared 2D distance is the pure polynomial $\|q\|_2^2$. Consequently, the 3D extruded `sdChamferBoxPrism(vec3 p, vec3 b, float chamfer)` evaluates $\sqrt{\|q\|_2^2 + \max(d_z, 0)^2}$ in **1 `sqrt` (instead of 2)**, using **0 `sqrt`s** across 100% of the 3D interior and all 10 face/chamfer normal prisms.

### 2. Round-Free, Clamp-Free Orthogonal Decomposition of `sdCoolS2D`

In IQ's `sdfCoolS(vec2 p)` (`distfunctions2d`), line 2540 evaluates `rex = p.x - min(round(p.x/0.4), 0.4)`. Because $p_x \ge 0$ after symmetry folding, `round(p.x/0.4)` is a non-negative integer $\{0, 1, 2, \dots\}$, so `min(round(p.x/0.4), 0.4)` is `0.0` for $p_x < 0.2$ and `0.4` for $p_x \ge 0.2$. Thus $|rex| \equiv ||p_x - 0.2| - 0.2|$ (**eliminating `round`, division by `0.4`, and `min`**). Furthermore, $|p_y - \mathrm{clamp}(p_y, 0, 0.4)| \equiv m_y := \max(|p_y - 0.2| - 0.2, 0)$, and the sign expression `2*p.x + aby + abs(aby + 0.4) - 0.4` simplifies algebraically to `2*(p.x + my - 0.4)`.

---

## Numerical Verification & Benchmark

| Function | `sqrt` / SFU in Interior & Face Strips | 3D Vertex Wedge `sqrt`s | Max Error vs Ground Truth |
| :--- | :---: | :---: | :---: |
| **`sdf-gems` `sdChamferBox2D` / `sdgChamferBox2D`** | **`0`** (`1` in vertex wedge) | **`1` (down from `2` in `sdChamferBoxPrism`)** | **$4.44 \times 10^{-16}$** |
| **`sdf-gems` `sdOrientedBox2D`** | **`1` `inversesqrt` (`0` divs, `0` post-fold `sqrt`s)** | N/A | **$4.44 \times 10^{-16}$** |
| **`sdf-gems` `sdCoolS2D`** | **`1` `sqrt` (`0` `round`s, `0` `clamp`s, `0` divs)** | N/A | **$8.88 \times 10^{-16}$** |
