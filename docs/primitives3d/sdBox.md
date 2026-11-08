# `sdBox` — Exact 0/1-`sqrt` 3D Axis-Aligned Box (0-`sqrt` Interior & 6 Face Prisms)

| 3D Surface ($`f(p) = 0`$) | Primary 2D Cross-Section | Secondary 2D Cross-Section |
| :---: | :---: | :---: |
| ![sdBox](../images/primitives3d/sdBox.png) | ![sdBox slice 1](../images/primitives3d/sdBox_slice.png) | ![sdBox slice 2](../images/primitives3d/sdBox_slice2.png) |

| Property | Value |
| :--- | :--- |
| **Author / Discoverer** | **Inigo Quilez** — upgraded to 0-`sqrt` median face-prism check using **Gemini 4 Argon** (`sdf-gems`, 2026) |
| **Source** | [https://iquilezles.org/articles/distfunctions/](https://iquilezles.org/articles/distfunctions/) |
| **Category** | `primitives3d` |
| **GLSL File** | [`glsl/primitives3d/sdBox.glsl`](../../glsl/primitives3d/sdBox.glsl) |
| **Exact Euclidean SDF?** | **Yes** ($`\Vert \nabla d\Vert = 1`$, `0` `sqrt`s in interior + all 6 exterior face prisms [$`64.9\%`$ of space], `1` in edges/corners) |

---

## Before & After Comparison

![0-SFU Median Face-Prism Upgrade](../images/comparisons/box_0sfu_face_prism_comparison.png)

---

## Mathematical Overview

For $`\mathbf{q} = |\mathbf{p}| - \mathbf{b}`$, let $`q_{\text{med}} = \max(\min(q_x, q_y), \min(\max(q_x, q_y), q_z))`$ be the second-largest component of $`\mathbf{q}`$. Whenever $`q_{\text{med}} \le 0`$ (meaning **at most ONE** component of $`\mathbf{q}`$ is positive, covering **100% of the box interior AND all 6 infinite exterior face prisms**), the exact Euclidean distance is identically:

```math
\|\max(\mathbf{q}, \mathbf{0})\| + \min(q_{\max}, 0) \equiv \max(q_x, \max(q_y, q_z))
```

with **ZERO `sqrt`s**, evaluating `length(max(q, 0.0))` only in the 12 exterior edge and 8 corner wedges ($`q_{\text{med}} > 0`$).

---

## GLSL Implementation

```glsl
float sdBox( vec3 p, vec3 b )
{
    vec3 q = abs(p) - b;
    if( max(min(q.x, q.y), min(max(q.x, q.y), q.z)) <= 0.0 )
        return max(q.x, max(q.y, q.z));
    return length(max(q, 0.0));
}
```
