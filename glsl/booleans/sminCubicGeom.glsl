// Corner-Local Associative C2 (G2 Curvature-Continuous) Smooth Minimum (sminCubicGeom)
// Discovered using Gemini 4 Argon (sdf-gems, 2026)
// Reference: Inigo Quilez, https://iquilezles.org/articles/smin/
//
// Solves the 5-property Smooth-Minimum problem:
//   - Polynomial sminCubic (CD family) is C2-continuous, but FAILS Locality (its
//     blend strip |a - b| < 6k extends to infinity along the medial bisector a = b,
//     underestimating distance by 29.3% forever) and FAILS 3-way Associativity.
//   - Geometrical sminCircularGeom (L2 quarter-circle) is Rigid, Corner-Local,
//     Conservative, and Associative, but FAILS C2 / G2 Curvature Continuity (its
//     second derivative jumps by 1/R at the blend boundaries a = R and b = R,
//     causing visible Mach-band creases in specular reflections).
//
// Mathematical Formulation (L3 Lame Corner Fillet):
//   Using the L3 norm ||u||_3 = (u.x^3 + u.y^3)^(1/3) on u = max(R - vec2(a, b), 0)
//   with corner support radius R = k / (1 - 2^(-1/3)) ≈ 4.84732210 * k:
//     1. Strictly Corner-Local & Rigid: Identically equals min(a, b) outside a < R, b < R.
//     2. C2 / G2 Curvature-Continuous: Near u.x -> 0+, (u.x^3 + u.y^3)^(1/3) =
//        u.y + u.x^3 / (3 * u.y^2) + O(u.x^6); both linear and quadratic terms vanish,
//        so the gradient and full Hessian match min(a, b) with zero curvature jump at
//        a = R and b = R, and the profile is C^infinity across the diagonal a = b.
//     3. Holder-Guaranteed Conservative (for grad(a) perp grad(b)):
//        ||grad f||_2 = (||u||_4 / ||u||_3)^2 in [2^(-1/6), 1] = [0.8908987, 1.0],
//        never overestimating distance and cutting max Eikonal error by 2.7x vs sminCubic.
//     4. Exact 3-Way Associativity: ((u.x^3 + u.y^3) + u.z^3)^(1/3) is order-independent.

#ifndef SDF_GEMS_SMIN_CUBIC_GEOM_GLSL
#define SDF_GEMS_SMIN_CUBIC_GEOM_GLSL

// 1. Corner-Local C2 Smooth Union (normalized so sminCubicGeom(k, k, k) == 0)
float sminCubicGeom(float a, float b, float k) {
    k *= 4.84732210; // 1.0 / (1.0 - pow(2.0, -1.0 / 3.0))
    vec2 u = max(k - vec2(a, b), 0.0);
    return max(k, min(a, b)) - pow(u.x * u.x * u.x + u.y * u.y * u.y, 1.0 / 3.0);
}

// 2. Corner-Local C2 Smooth Intersection (De Morgan dual)
float smaxCubicGeom(float a, float b, float k) {
    return -sminCubicGeom(-a, -b, k);
}

// 3. Corner-Local C2 Smooth Subtraction (subtracts d1 from d2 with C2 corner fillet)
float opSmoothSubtractionCubicGeom(float d1, float d2, float k) {
    return -sminCubicGeom(d1, -d2, k);
}

// 4. Corner-Local C2 Smooth Union with Analytical 3D Gradient (vec4(d, grad.xyz))
vec4 sdgSmoothUnionCubic(vec4 a, vec4 b, float k) {
    k *= 4.84732210;
    vec2 u = max(k - vec2(a.x, b.x), 0.0);
    float s3 = u.x * u.x * u.x + u.y * u.y * u.y;
    if (s3 <= 0.0) return (a.x < b.x) ? a : b;
    float l3 = pow(s3, 1.0 / 3.0);
    vec2 w = (u * u) / (l3 * l3);
    return vec4(k - l3, w.x * a.yzw + w.y * b.yzw);
}

#endif // SDF_GEMS_SMIN_CUBIC_GEOM_GLSL
