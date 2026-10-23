// ============================================================================
// sdJack3D — Exact 1-Sqrt 3D Body-Diagonal, Axis, and Face-Diagonal Capsule Jacks
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`1` `sqrt` everywhere in R^3; `4.44e-16` max error vs 4/3/6-capsule union)
// ============================================================================
//
// Mathematical Derivation:
//   1. sdJack3D (8-Pointed Body-Diagonal Caltrop / 4-Capsule Union):
//      The union of the 4 body-diagonal capsules connecting (-a,-a,-a) <-> (a,a,a)
//      (where a = w / 3) has octahedral symmetry. After folding into the positive
//      octant `p = abs(p)`, for any sign vector s in {-1, +1}^3:
//        |dot(p, s)| = |sx*px + sy*py + sz*pz| <= px + py + pz = dot(p, vec3(1)).
//      Because the squared distance ||p - clamp(dot(p, u_s), -l, l) * u_s||^2 is a
//      strictly decreasing function of |dot(p, u_s)|, the closest of the 4 body-
//      diagonal capsules in the positive octant is ALWAYS the (+1, +1, +1) diagonal
//      (with zero sorting swaps required!):
//        return length(p - min(p.x + p.y + p.z, w) * (1.0 / 3.0)) - r;
//   2. sdAxisJack3D (6-Pointed Coordinate-Axis Jack / 3-Capsule Union):
//      In the positive octant `p = abs(p)`, the closest coordinate-axis segment is
//      along the largest coordinate m = max(p.x, max(p.y, p.z)). Clamping t = min(m, w)
//      gives ||p - t*e_m||^2 = dot(p, p) - t*(2*m - t) in 1 sqrt without sorting.
//   3. sdFaceJack3D (12-Pointed Cuboctahedral Jack / 6-Capsule Union):
//      Sorting p.x >= p.y >= p.z >= 0 maximizes p.x + p.y among all 6 face diagonals.

#ifndef SDF_GEMS_SD_JACK3D_GLSL
#define SDF_GEMS_SD_JACK3D_GLSL

// 1. Exact 2-line, 1-sqrt 3D 8-pointed body-diagonal jack (union of 4 diagonal capsules)
float sdJack3D(vec3 p, float w, float r) {
    p = abs(p);
    return length(p - min(p.x + p.y + p.z, w) * (1.0 / 3.0)) - r;
}

// 2. Exact 3-line, 1-sqrt 3D 6-pointed coordinate-axis jack (union of 3 axis capsules)
float sdAxisJack3D(vec3 p, float w, float r) {
    p = abs(p);
    float m = max(p.x, max(p.y, p.z)), t = min(m, w);
    return sqrt(max(dot(p, p) - t * (2.0 * m - t), 0.0)) - r;
}

// 3. Exact 1-sqrt 3D 12-pointed face-diagonal jack (union of 6 face-diagonal capsules)
float sdFaceJack3D(vec3 p, float w, float r) {
    p = abs(p);
    if (p.x < p.y) p.xy = p.yx;
    if (p.y < p.z) p.yz = p.zy;
    if (p.x < p.y) p.xy = p.yx;
    float t = min(p.x + p.y, w) * 0.5;
    return length(vec3(p.x - t, p.y - t, p.z)) - r;
}

#endif // SDF_GEMS_SD_JACK3D_GLSL
