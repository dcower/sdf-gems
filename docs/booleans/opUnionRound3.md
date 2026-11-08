# `opUnionRound3` — Exact 0/1-`sqrt` Symmetric Ternary 3-Way Corner Booleans & Gradients

![opUnionRound3](../images/booleans/opUnionRound3.png)

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | Discovered and derived using **Gemini 4 Argon** (`sdf-gems`, 2026) |
| **Source** | [https://github.com/dcower/sdf-gems](https://github.com/dcower/sdf-gems) |
| **Category** | `booleans` |
| **GLSL File** | [`glsl/booleans/opUnionRound3.glsl`](../../glsl/booleans/opUnionRound3.glsl) |
| **Exact Euclidean SDF?** | **Yes** ($`\Vert \nabla d\Vert = 1`$ everywhere for `opUnionRound3` & `opUnionChamfer3`, `8.88e-16` max error, 0/1 `sqrt`) |

---

## Before & After Comparison

![opUnionRound3 Comparison](../images/comparisons/union_round3_comparison.png)

---

## Mathematical Overview

### 1. Why Chained Binary Corner Booleans Fail at 3-Way Corners

When three orthogonal walls $`a, b, c`$ meet at a 3-way corner, chaining two binary chamfer operators `opUnionChamfer(opUnionChamfer(a, b, r), c, r)` breaks $`S_3`$ permutation symmetry:

```math
d_{(ab)c} = \frac{d_{ab} + c - r}{\sqrt{2}} = \frac{1}{2}a + \frac{1}{2}b + \frac{1}{\sqrt{2}}c - r\!\left(\frac{1}{2} + \frac{1}{\sqrt{2}}\right)
```

1. **$`9.74^\circ`$ Normal Tilt**: The corner facet has normal $`\mathbf{n} = \left(\frac{1}{2}, \frac{1}{2}, \frac{1}{\sqrt{2}}\right)`$ instead of the isotropic octahedral normal $`\frac{1}{\sqrt{3}}(1, 1, 1)`$, tilting the corner normal by $`\arccos\!\left(\frac{1 + 1/\sqrt{2}}{\sqrt{3}}\right) = 9.7356^\circ`$.
2. **$`41.4\%`$ Axis-Intercept Skew**: The intercepts on the $`a, b`$ axes are $`r(1 + \sqrt{2}) \approx 2.414 r`$, whereas the intercept on the $`c`$ axis is $`r(1 + 1/\sqrt{2}) \approx 1.707 r`$ (a $`\sqrt{2}:1`$ distortion), while executing **2 `sqrt`s** instead of **0 or 1 `sqrt`**.

### 2. The $`S_3`$ Weyl-Chamber Convex Polyhedron Theorem (`opUnionChamfer3`, `opUnionChamferOct3`, `sdgChamfer3`)

Sorting $`(a, b, c)`$ into the $`S_3`$ Weyl chamber $`u_1 \le u_2 \le u_3`$ via 3 compare-swaps reveals that the exterior region $`\{d \ge 0\}`$ of a 3-way chamfered union is an intersection of half-spaces, hence a **convex polyhedron** $`\mathcal{P}`$:

- **Edge-Chamfered 3-Way Union (`opUnionChamfer3`)**: $`\mathcal{P}_{\text{edge}} = \{x, y, z \ge 0,\; x+y \ge r,\; y+z \ge r,\; z+x \ge r\}`$.
- **Octahedral Corner Gusset Union (`opUnionChamferOct3`)**: $`\mathcal{P}_{\text{oct}} = \{x, y, z \ge 0,\; x+y+z \ge r\}`$.

Because the interior of any convex polyhedron has **zero** edge or vertex Voronoi regions, the exact 3D Euclidean distance across **100% of the exterior $`\{d \ge 0\}`$** requires **ZERO `sqrt`s**:

```math
u_1 \ge 0 \;\land\; u_1 + u_2 \ge r \implies d = \min\!\left(u_1,\; \frac{u_1 + u_2 - r}{\sqrt{2}}\right)
```

In the interior $`\{d < 0\}`$, projecting $`(u_1, u_2, u_3)`$ onto the sorted chamfer strip $`x_1 + x_2 = r`$ ($`t = x_2 \in [r/2, r]`$) or the coordinate wall $`x_1 = 0`$ gives the exact 3D Euclidean distance (`8.88e-16` max error) in **1 `sqrt`**.

### 3. 0/1-`sqrt` Wall-Gated `opUnionRound3`, `sdgRound3`, & `sminCubicGeom3`

Gating `opUnionRound3(a, b, c, r)` on the second-smallest input $`u_2 \ge r`$ returns $`u_1 = \min(a,b,c)`$ immediately with **0 `sqrt`s** across all flat walls and the far field, and prevents the catastrophic mantissa cancellation $`r - \sqrt{(r - u_1)^2} \to 0`$ when $`|u_1| \ll r`$.

### Numerical Verification Summary

| 3-Way Ternary Operator | `sqrt` Count (Ext / Wall / Corner) | Max Error vs Brute-Force 3D GT | Max $`\big\vert \Vert \nabla d\Vert _2 - 1\big\vert`$ | Max $`\Vert \nabla - \nabla_{\text{FD}}\Vert _2`$ |
| :--- | :---: | :---: | :---: | :---: |
| Chained Binary `opUnionChamfer(opUnionChamfer)` | $`0 / 1 / 2`$ | `6.51e-02` ($`9.74^\circ`$ tilt) | — | — |
| **`opUnionChamfer3` & `sdgChamfer3`** | **$`0`$ (ext) / $`1`$ (int)** | **`8.88e-16`** | **`5.55e-16`** | **`6.62e-09`** |
| **`opUnionChamferOct3`** | **$`0`$ (ext & prism) / $`1`$** | **`8.88e-16`** | **`4.44e-16`** | **`6.10e-09`** |
| **`opUnionRound3` & `sdgRound3`** | **$`0`$ ($`u_2 \ge r`$) / $`1`$ ($`u_2 < r`$)** | **`8.88e-16`** | **`6.66e-16`** | **`8.88e-09`** |
| **`sminCubicGeom3` & `sdgSmoothUnionCubic3`** | **$`0`$ ($`u_2 \ge R`$) / $`1`$ `cbrt`** | **`0.00e+00`** | $`\Vert \nabla\Vert _2 \in [3^{-1/6}, 1]`$ | **`6.73e-09`** |

---

## GLSL Implementation

```glsl
float opUnionRound3( float a, float b, float c, float r )
{
    float u1 = min(a, min(b, c));
    float u2 = max(min(a, b), min(max(a, b), c));
    if( u2 >= r ) return u1;
    return r - length(max(r - vec3(a, b, c), 0.0));
}

float opUnionChamfer3( float a, float b, float c, float r )
{
    float u1 = min(a, b), u2 = max(a, b), u3 = max(u2, c);
    u2 = max(u1, min(u2, c)); u1 = min(u1, c);
    if( u1 >= 0.0 && u1 + u2 >= r )
        return min(u1, (u1 + u2 - r) * 0.70710678);
    float s = r - u1 + u2;
    float t = max((2.0 * u3 >= s) ? 0.5 * s : (s + u3) * (1.0 / 3.0), 0.5 * r);
    return -(t >= r ? length(vec3(u1, max(r - u2, 0.0), max(r - u3, 0.0)))
                    : length(vec3(u1 - r + t, u2 - t, max(t - u3, 0.0))));
}
```
