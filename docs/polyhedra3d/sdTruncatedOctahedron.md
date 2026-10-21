# `sdTruncatedOctahedron` — 3D Truncated Octahedron / Kelvin Cell (GDF Plane Bound)

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdTruncatedOctahedron](../images/polyhedra3d/sdTruncatedOctahedron.png) | ![sdTruncatedOctahedron slice 1](../images/polyhedra3d/sdTruncatedOctahedron_slice.png) | ![sdTruncatedOctahedron slice 2](../images/polyhedra3d/sdTruncatedOctahedron_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Mercury Demogroup** |
| **Source** | [https://mercury.sexy/hg_sdf/](https://mercury.sexy/hg_sdf/) |
| **Category** | `polyhedra3d` |
| **GLSL File** | [`glsl/polyhedra3d/sdTruncatedOctahedron.glsl`](../../glsl/polyhedra3d/sdTruncatedOctahedron.glsl) |
| **Exact Euclidean SDF?** | No (Plane Bound) |

---

## Mathematical Overview

Evaluates the maximum of the 14 bounding planes of the Truncated Octahedron (Kelvin cell).

---

## GLSL Implementation

```glsl
float sdTruncatedOctahedron( vec3 p, float r )
{
    p = abs(p);
    return max(max(p.x, max(p.y, p.z)) - r, (p.x + p.y + p.z - 1.5*r) * 0.57735027);
}
```
