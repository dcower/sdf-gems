# `sdSphericalSegmentSym` — Exact 3D Two-Plane Cut Sphere (Spherical Segment)

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdSphericalSegmentSym](../images/primitives3d/sdSphericalSegmentSym.png) | ![sdSphericalSegmentSym slice](../images/primitives3d/sdSphericalSegmentSym_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Discovered using Gemini 4 Argon (sdf-gems, 2026)** |
| **Source** | [https://github.com/dcower/sdf-gems](https://github.com/dcower/sdf-gems) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdSphericalSegmentSym.glsl`](../../glsl/primitives3d/sdSphericalSegmentSym.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Exact Euclidean SDF of a sphere of radius $r$ truncated symmetrically by parallel planes $y = \pm h$.

---

## GLSL Implementation

```glsl
float sdSphericalSegmentSym( vec3 p, float r, float h )
{
    float w = sqrt(r*r - h*h);
    vec2  q = vec2( length(p.xz), abs(p.y) );
    float s = max( w*q.y - h*q.x, (r - h)*(r - h + 2.0*q.y) - q.x*q.x );
    return (s < 0.0) ? length(q) - r :
           (q.x < w) ? q.y - h       :
                       length(q - vec2(w, h));
}
```
