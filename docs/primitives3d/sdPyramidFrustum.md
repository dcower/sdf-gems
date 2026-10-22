# `sdPyramidFrustum` & `sdgPyramidFrustum` — Exact 0/1-Sqrt 3D Square Frustum & Analytical Gradient

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdPyramidFrustum](../images/primitives3d/sdPyramidFrustum.png) | ![sdPyramidFrustum slice 1](../images/primitives3d/sdPyramidFrustum_slice.png) | ![sdPyramidFrustum slice 2](../images/primitives3d/sdPyramidFrustum_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://github.com/dcower/sdf-gems](https://github.com/dcower/sdf-gems) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdPyramidFrustum.glsl`](../../glsl/primitives3d/sdPyramidFrustum.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1` everywhere in $\mathbb{R}^3$; `0` per-pixel `sqrt`s in interior, caps, and slanted faces) |

---

## Mathematical Overview & Tangential Hip-Edge Pruning Theorem

![0/1-Sqrt Frustums & Barrels with Tangential Hip-Edge Pruning](../images/comparisons/frustum_barrel_solid_angle_comparison.png)

Computes the exact 3D Euclidean SDF and analytical unit gradient `vec4(d, grad.xyz)` of a square frustum of bottom half-width $b_1$, top half-width $b_2$, and half-height $h$:
1. **0-`sqrt` Three-Zone Euclidean Decomposition**:
   After folding $p_{xz} = |p_{xz}|$ into $p_x \ge p_z \ge 0$, let $b_m = \frac{1}{2}(b_1 + b_2)$, $k = \frac{1}{2}(b_1 - b_2)$, $q = p - (b_m, 0, b_m)$, and $m_2 = h^2 + k^2$. In the entire interior ($\max(d_s, d_c) \le 0$), the exterior top/bottom cap prisms ($d_c > 0 \;\land\; p_x \le b_{\text{cap}}$), and the exterior slanted-face normal prism ($d_s > 0 \;\land\; |h q_y - k q_x| \le m_2 \;\land\; q_z + k u_1 \le 0$), the exact signed distance and unit gradient are purely linear (**0 per-pixel `sqrt`s**; `35.6%` of $[-2, 2]^3$).
2. **Tangential Hip-Edge Pruning ($q_z + k u_1 \le 0$)**:
   Let $u_1 = \operatorname{clamp}((h q_y - k q_x)/m_2, -1, 1)$ minimize $(q_x + k u)^2 + (q_y - h u)^2$ on $[-1, 1]$. Whenever $q_z + k u_1 \le 0$, $\max(q_z + k u_1, 0) = 0$, so for every $u_2 \in [-1, 1]$:
   $$d_1 = (q_x + k u_1)^2 + (q_y - h u_1)^2 \le (q_x + k u_2)^2 + (q_y - h u_2)^2 + (q_z + k u_2)^2 = d_2.$$
   Thus the corner hip edge $u_2$ is **100% dominated ($d_1 \le d_2$ identically, `0.0` difference)** whenever $q_z + k u_1 \le 0$ and never needs to be evaluated.

---

## GLSL Implementation

```glsl
float sdPyramidFrustum( vec3 p, float b1, float b2, float h )
{
    p.xz = abs(p.xz);
    p.xz = (p.z > p.x) ? p.zx : p.xz;

    float bm = 0.5 * (b1 + b2), k = 0.5 * (b1 - b2);
    vec3  q  = p - vec3(bm, 0.0, bm);
    float m2 = h * h + k * k;
    float ds = (h * q.x + k * q.y) * inversesqrt(m2);
    float dc = abs(q.y) - h;
    float u1 = clamp((h * q.y - k * q.x) / m2, -1.0, 1.0);
    if( max(ds, dc) <= 0.0 || (dc > 0.0 && p.x <= ((p.y < 0.0) ? b1 : b2)) ||
        (ds > 0.0 && abs(h * q.y - k * q.x) <= m2 && q.z + k * u1 <= 0.0) )
        return max(ds, dc);

    vec3 d1 = vec3(q.x + k * u1, q.y - h * u1, max(q.z + k * u1, 0.0));
    if( q.z + k * u1 <= 0.0 )
        return length(d1.xy);

    float u2 = clamp((h * q.y - k * (q.x + q.z)) / (m2 + k * k), -1.0, 1.0);
    vec3  d2 = vec3(q.x + k * u2, q.y - h * u2, q.z + k * u2);
    return sqrt(min(dot(d1, d1), dot(d2, d2)));
}
```
