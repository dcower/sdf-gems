# `sdRectFrustum` & `sdgRectFrustum` — Exact 0/1-Sqrt 3D Rectangular & Hip-Roof Frustum & Gradient

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdRectFrustum](../images/primitives3d/sdRectFrustum.png) | ![sdRectFrustum slice 1](../images/primitives3d/sdRectFrustum_slice.png) | ![sdRectFrustum slice 2](../images/primitives3d/sdRectFrustum_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://github.com/dcower/sdf-gems](https://github.com/dcower/sdf-gems) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdRectFrustum.glsl`](../../glsl/primitives3d/sdRectFrustum.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1` everywhere in $\mathbb{R}^3$; `0` per-pixel `sqrt`s in interior, caps, and slanted faces) |

---

## Mathematical Overview

![0/1-Sqrt Frustums & Barrels with Tangential Hip-Edge Pruning](../images/comparisons/frustum_barrel_solid_angle_comparison.png)

Unifies rectangular boxes, rectangular pyramids, keycap frustums, and hip roofs into a **0/1-`sqrt` exact 3D Euclidean SDF (`sdRectFrustum`) and analytical gradient (`sdgRectFrustum`)**, using the **Tangential Hip-Edge Pruning Theorem** so the corner hip edge $u_e$ is evaluated only in the exterior corner wedge where both $q_z + k_y u_x > 0$ and $q_x + k_x u_z > 0$.

---

## GLSL Implementation

```glsl
float sdRectFrustum( vec3 p, vec2 b1, vec2 b2, float h )
{
    p.xz = abs(p.xz);
    vec2  bm  = 0.5 * (b1 + b2), k = 0.5 * (b1 - b2);
    vec3  q   = p - vec3(bm.x, 0.0, bm.y);
    float mx2 = h * h + k.x * k.x, mz2 = h * h + k.y * k.y;
    float dsx = (h * q.x + k.x * q.y) * inversesqrt(mx2);
    float dsz = (h * q.z + k.y * q.y) * inversesqrt(mz2);
    float dc  = abs(q.y) - h;
    float df  = max(max(dsx, dsz), dc);

    float ux  = clamp((h * q.y - k.x * q.x) / mx2, -1.0, 1.0);
    float uz  = clamp((h * q.y - k.y * q.z) / mz2, -1.0, 1.0);
    vec2  bc  = (p.y < 0.0) ? b1 : b2;
    if( df <= 0.0 || (dc > 0.0 && p.x <= bc.x && p.z <= bc.y) ||
        (dsx > 0.0 && abs(h * q.y - k.x * q.x) <= mx2 && q.z + k.y * ux <= 0.0) ||
        (dsz > 0.0 && abs(h * q.y - k.y * q.z) <= mz2 && q.x + k.x * uz <= 0.0) )
        return df;

    vec3  vx = vec3(q.x + k.x * ux, q.y - h * ux, max(q.z + k.y * ux, 0.0));
    vec3  vz = vec3(max(q.x + k.x * uz, 0.0), q.y - h * uz, q.z + k.y * uz);
    float d  = min(dot(vx, vx), dot(vz, vz));
    if( q.z + k.y * ux > 0.0 && q.x + k.x * uz > 0.0 )
    {
        float ue = clamp((h * q.y - dot(k, q.xz)) / (h * h + dot(k, k)), -1.0, 1.0);
        vec3  ve = vec3(q.x + k.x * ue, q.y - h * ue, q.z + k.y * ue);
        d = min(d, dot(ve, ve));
    }
    return sqrt(d);
}
```
