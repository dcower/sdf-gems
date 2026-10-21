# `sdCutSphere` — Exact 3D Single-Plane Cut Sphere

| 3D Surface ($f(p) = 0$) | 2D Cross-Section Slice |
| :---: | :---: |
| ![sdCutSphere](../images/primitives3d/sdCutSphere.png) | ![sdCutSphere slice](../images/primitives3d/sdCutSphere_slice.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdCutSphere.glsl`](../../glsl/primitives3d/sdCutSphere.glsl) |
| **Exact Euclidean SDF?** | Yes (`||grad d|| = 1`) |

---

## Mathematical Overview

Sphere of radius $r$ truncated by the horizontal plane $y = h$.

---

## GLSL Implementation

```glsl
float sdCutSphere( vec3 p, float r, float h )
{
    float w = sqrt(r*r - h*h);
    vec2 q = vec2( length(p.xz), p.y );
    float s = max( (h - r)*q.x*q.x + w*w*(h + r - 2.0*q.y), h*q.x - w*q.y );
    return (s < 0.0) ? length(q) - r :
           (q.x < w) ? h - q.y       :
                       length(q - vec2(w,h));
}
```
