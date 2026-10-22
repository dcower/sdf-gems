# `sdgCappedCone` — Exact 0-Division, 0/1/2-SFU 3D Capped Cone Distance & Analytical Gradient

![sdgCappedCone Preview](../images/gradients3d/sdgCappedCone.png)

## Overview

Inigo Quilez's `sdgCappedCone` ([`distgradfunctions3d`](https://iquilezles.org/articles/distgradfunctions3d/), Shadertoy `3ctcW7`, 2025) returns:
```glsl
return (la<lb)?
    vec4(s*sqrt(la), 0.0, sign(p.y), 0.0) :
    vec4(s*sqrt(lb), vec3(k.y*p.xz/l,-k.x).xzy/sqrt(m));
```
We uncover two critical defects in this baseline implementation:
1. **Missing Corner-Rim Gradient Formula (Error up to $\sqrt{2} \approx 1.414$ across $>40\%$ of `[-2, 2]^3`)**: In both the exterior top-rim corner wedge ($l > r_2, \mathbf{q}\cdot\mathbf{k} \le 0$) and bottom-rim corner wedge ($l > r_1, (\mathbf{q}+\mathbf{k})\cdot\mathbf{k} \ge 0$), `a` and `b` clamp to the exact same rim point, so `la == lb` (`la < lb` is `false` or noisy to 1 ULP). IQ's code takes the second branch and returns the **constant slanted-wall normal** (or chatters between the cap normal and wall normal on the bottom rim) instead of the toroidal corner-rim normal!
2. **Symmetry-Axis `NaN` (`l = 0`)**: On the symmetry axis $p_x = p_z = 0$, whenever `la >= lb`, IQ divides by `l = 0.0` (`NaN`).

`sdgCappedCone` fixes both bugs while eliminating all clamped divisions and adding a **0-SFU Cap-Voronoi early exit**.

---

## Mathematical Formulation

![Capped Cone Gradient Comparison](../images/comparisons/capped_cone_gradient_comparison.png)

| Implementation | Max Distance Error | Max Smooth Gradient Error | Axis $p_{xz} = 0$ Safety | Clamped Divisions | Cap-Zone SFU Count |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **IQ `sdgCappedCone` (`3ctcW7`)** | `8.88e-16` | **`1.414e+00` (Broken Corner Rim)** | `NaN` (`0/0`) | `1` | `2` `sqrt`s |
| **Ours `sdgCappedCone`** | **`8.88e-16`** | **`4.31e-10`** | **100% `NaN`-Free** | **`0`** | **`0` SFUs** |

---

## GLSL Implementation

See [`glsl/gradients3d/sdgCappedCone.glsl`](../../glsl/gradients3d/sdgCappedCone.glsl).
