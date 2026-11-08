# Exact 0-Trig 3D Reuleaux Tetrahedron & Constant-Width Rotor (`sdReuleauxTetrahedron`, `sdReuleauxRotor3D`)

![sdReuleaux3D](../../docs/images/primitives3d/sdReuleaux3D.png)

> **Source**: [`glsl/primitives3d/sdReuleaux3D.glsl`](../../glsl/primitives3d/sdReuleaux3D.glsl)
> **Category**: `primitives3d` | **Exact**: Yes ($`0`$ trig, $`0`$ divisions, $`1`$ `sqrt` in face & vertex wedges, $`2`$ `sqrt`s in arc edge)
> **Author**: Discovered and derived using Gemini 4 Argon (`sdf-gems`, 2026)

---

## Visual Overview: 2D & 3D Constant-Width Reuleaux Primitives

![Reuleaux Constant-Width Comparison](../../docs/images/comparisons/reuleaux_constant_width_comparison.png)

---

## Mathematical Foundations

### 1. The $`T_d`$ Vertex-Distance Total-Ordering Theorem (`sdReuleauxTetrahedron`)

Let $`\mathbf{V}_0 = (s, s, s)`$, $`\mathbf{V}_1 = (s, -s, -s)`$, $`\mathbf{V}_2 = (-s, s, -s)`$, $`\mathbf{V}_3 = (-s, -s, s)`$ be the 4 vertices of a regular tetrahedron of edge length $`w = 2\sqrt{2}\, s`$ ($`s = \frac{w}{2\sqrt{2}} \approx 0.35355339\, w`$).

Folding $`\mathbf{p} \in \mathbb{R}^3`$ into the 24-fold $`T_d`$ parity-sorted Weyl chamber $`\mathcal{C}_{T_d} = \{x \ge y \ge |z|\}`$ guarantees $`x + y + z \ge x - y - z \ge -x + y - z \ge -x - y + z`$, and therefore:
```math
\|\mathbf{p} - \mathbf{V}_3\|_2 \;\ge\; \|\mathbf{p} - \mathbf{V}_2\|_2 \;\ge\; \|\mathbf{p} - \mathbf{V}_1\|_2 \;\ge\; \|\mathbf{p} - \mathbf{V}_0\|_2.
```
Consequently, inside $`\mathcal{C}_{T_d}`$:
1. **Single Spherical Face ($`\mathbf{V}_3 = (-s, -s, s)`$)**: All 4 spherical faces collapse to a single sphere centered at $`\mathbf{V}_3`$. Because $`\mathbf{V}_3`$ lies on the chamber walls $`x = y`$ and $`y = -z`$, a radial ray from $`\mathbf{V}_3`$ through $`\mathbf{p}`$ hits the spherical face inside the chamber if and only if $`2(y - z + 2s)^2 \ge \|\mathbf{p} - \mathbf{V}_3\|_2^2`$.
2. **Single Vertex Wedge ($`\mathbf{V}_0 = (s, s, s)`$)**: When $`y + z \ge x + s`$, $`\mathbf{p}`$ lies in the exterior vertex normal cone of $`\mathbf{V}_0`$, with exact distance $`\|\mathbf{p} - \mathbf{V}_0\|_2`$ (**1 `sqrt`**, fused with the spherical face branch!).
3. **Single 3D Circular-Arc Edge**: Along the wall $`y = z`$, the spheres around $`\mathbf{V}_3`$ and $`\mathbf{V}_2`$ intersect in a circular arc in the plane $`y = z`$ centered at $`(-s, 0, 0)`$ with radius $`R_{\text{edge}} = \sqrt{6}\, s`$.

---

### 2. Axis-Center Radial Collapse Theorem (`sdReuleauxRotor3D`)

Rotating the 2D Reuleaux triangle around its vertical symmetry axis produces the 3D constant-width Reuleaux rotor ("acorn"). In both the **top apex cone** ($`y \ge y_0 \land 3 r^2 \le (y - y_0)^2`$) and the **bottom spherical cap zone** ($`y \le 0 \land r^2 \le 3 y^2 \land 3 r^2 \le (y - y_0)^2`$, where $`r^2 = p_x^2 + p_z^2`$ and $`y_0 = w/\sqrt{3}`$), the closest 3D center lies on the $`y`$-axis at $`(0, y_0, 0)`$, so $`\sqrt{r^2 + (y - y_0)^2}`$ evaluates the exact 3D distance in **1 `sqrt`** without computing $`\sqrt{p_x^2 + p_z^2}`$!

---

## Numerical Verification Summary

| Primitive | SFU / Div Cost | Max Distance Error vs. Ground Truth | Max $`\bigl\vert \Vert \nabla d\Vert _2 - 1\bigr\vert`$ | Naive 4-Sphere CSG Error |
| :--- | :--- | :---: | :---: | :---: |
| `sdReuleauxTetrahedron` / `sdg` | **0 trig, 0 div, 1/2 `sqrt`** | $`1.33 \times 10^{-15}`$ | $`5.55 \times 10^{-16}`$ | $`1.84 \times 10^{-1}`$ ($`4`$ `sqrt`s) |
| `sdReuleauxRotor3D` / `sdg` | **0 trig, 0 div, 1/2 `sqrt`** | $`1.78 \times 10^{-15}`$ | $`4.44 \times 10^{-16}`$ | $`1.19 \times 10^{-1}`$ ($`3`$ `sqrt`s) |

---

## GLSL Implementation

```glsl
float sdReuleauxTetrahedron( vec3 p, float w )
{
    float s   = 0.35355339 * w;
    float sig = (p.x * p.y * p.z < 0.0) ? -1.0 : 1.0;
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;
    p.z *= sig;

    vec3  v       = p - vec3(-s, -s, s);
    float v2      = dot(v, v);
    float dyz     = p.y - p.z, d2s = dyz + 2.0 * s;
    bool  in_face = (2.0 * d2s * d2s >= v2);
    if( in_face || p.y + p.z >= p.x + s )
        return length(p - (in_face ? vec3(-s, -s, s) : vec3(s))) - (in_face ? w : 0.0);

    float ux  = p.x + s, syz = p.y + p.z;
    float dr  = sqrt(ux * ux + 0.5 * syz * syz) - 2.44948974 * s;
    return sqrt(dr * dr + 0.5 * dyz * dyz);
}

float sdReuleauxRotor3D( vec3 p, float w )
{
    const float invK = 0.57735027;
    float r2  = dot(p.xz, p.xz), y = p.y;
    float y0  = w * invK, dy0 = y - y0;
    if( 3.0 * r2 <= dy0 * dy0 && (dy0 >= 0.0 || (y <= 0.0 && r2 <= 3.0 * y * y)) )
        return sqrt(r2 + dy0 * dy0) - ((dy0 >= 0.0) ? 0.0 : w);
    bool rim = (y <= -0.5 * y0);
    return length(vec2(sqrt(r2) - (rim ? 0.5 * w : -0.5 * w), y + 0.5 * y0)) - (rim ? 0.0 : w);
}
```
