# `sdSphericalSegmentSym` — Exact 3D Two-Plane Cut Sphere (Conjugate-Square Reduction)

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdSphericalSegmentSym](../images/primitives3d/sdSphericalSegmentSym.png) | ![sdSphericalSegmentSym slice](../images/primitives3d/sdSphericalSegmentSym_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Discovered using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://github.com/dcower/sdf-gems](https://github.com/dcower/sdf-gems) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdSphericalSegmentSym.glsl`](../../glsl/primitives3d/sdSphericalSegmentSym.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`, `0` `sqrt`s in caps, `1` in sphere, `2` in rim) |

---

## Mathematical Overview

Computes the exact 3D Euclidean SDF of a sphere of radius $r$ truncated symmetrically by two parallel planes $y = \pm h$. By the **Conjugate-Square Rationalization Theorem** (see [`sdCutSphere`](sdCutSphere.md)), since $y = |p_y| \ge 0$ and $h > 0$, the cone test $w y < h q_x$ reduces to $w_2 y^2 < h^2 x_2$ with zero square roots, requiring **0 `sqrt`s** in the flat disk caps, **1 `sqrt`** in the spherical barrel region, and **2 `sqrt`s** in the circular rims.

---

## GLSL Implementation

```glsl
float sdSphericalSegmentSym( vec3 p, float r, float h )
{
    float x2 = dot(p.xz, p.xz), y = abs(p.y), w2 = r*r - h*h, d = x2 - w2, dy = y - h;
    if( w2*y*y < h*h*x2 && (r - h)*(r + h + 2.0*dy) < x2 )
        return sqrt(x2 + y*y) - r;
    return (d < 0.0) ? dy : sqrt(d*d / (x2 + w2 + 2.0*sqrt(w2*x2)) + dy*dy);
}
```
