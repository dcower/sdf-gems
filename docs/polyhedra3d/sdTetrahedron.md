# `sdTetrahedron` — 3D Regular Tetrahedron (GDF Plane Bound)

| 3D Surface ($f(p) = 0$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdTetrahedron](../images/polyhedra3d/sdTetrahedron.png) | ![sdTetrahedron slice 1](../images/polyhedra3d/sdTetrahedron_slice.png) | ![sdTetrahedron slice 2](../images/polyhedra3d/sdTetrahedron_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Mercury Demogroup** |
| **Source** | [https://mercury.sexy/hg_sdf/](https://mercury.sexy/hg_sdf/) |
| **Category** | `polyhedra3d` |
| **GLSL File** | [`glsl/polyhedra3d/sdTetrahedron.glsl`](../../glsl/polyhedra3d/sdTetrahedron.glsl) |
| **Exact Euclidean SDF?** | No (Plane Bound) |

---

## Mathematical Overview

Evaluates the maximum of the 4 bounding planes of a regular tetrahedron.

---

## GLSL Implementation

```glsl
float sdTetrahedron( vec3 p, float r )
{
    return (max(abs(p.x + p.y) - p.z, abs(p.x - p.y) + p.z) - r) * 0.57735027;
}
```
