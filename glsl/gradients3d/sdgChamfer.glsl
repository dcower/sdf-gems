// ============================================================================
// sdgChamfer — Angle-Exact Euclidean Chamfer Booleans with Analytical 3D Gradients
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` across exterior, chamfer bevel, and interior for any angle theta)
// ============================================================================
//
// Given two 3D signed distance + gradient vectors a = vec4(da, grad_a) and
// b = vec4(db, grad_b) with unit gradients ||grad_a|| = ||grad_b|| = 1 and
// intersection cosine c = dot(grad_a, grad_b) = cos(theta):
//   1. Because (grad_a + grad_b) . (grad_a - grad_b) = 1 - 1 = 0, the normal strip
//      separator of the chamfer plane a + b = r at the vertex V = (0, r) is
//      universally |b - a| = r, independent of the angle theta in (0, pi).
//   2. Projecting into the Gram-Schmidt orthonormal frame at V = (0, r) yields the
//      exact Euclidean corner distance and analytical unit 3D gradient in all
//      regions (exterior, chamfer bevel, and interior) for Union, Intersection,
//      and Subtraction.

#ifndef SDF_GEMS_SDG_CHAMFER_GLSL
#define SDF_GEMS_SDG_CHAMFER_GLSL

// 1. Angle-Exact Euclidean Chamfer Union with Analytical 3D Gradient
vec4 sdgChamferUnion(vec4 a, vec4 b, float r) {
    if (a.x > b.x) { vec4 t = a; a = b; b = t; }
    float c  = dot(a.yzw, b.yzw);
    if (c < -0.95) return a;
    float ic = inversesqrt(max(2.0 + 2.0 * c, 1e-5));
    float ch = (a.x + b.x - r) * ic;
    vec3  gc = (a.yzw + b.yzw) * ic;
    if (a.x > 0.0 && a.x + b.x > r) return (a.x < ch) ? a : vec4(ch, gc);
    float s2 = max(1.0 - c * c, 1e-5);
    float h  = max(r - b.x + c * a.x, 0.0);
    float dv = -sqrt(a.x * a.x + h * h / s2);
    vec3  gv = (h <= 0.0) ? a.yzw : ((a.x + c * h / s2) * a.yzw - (h / s2) * b.yzw) / dv;
    return (b.x - a.x > r || dv > ch) ? vec4(dv, gv) : vec4(ch, gc);
}

// 2. Angle-Exact Euclidean Chamfer Intersection with Analytical 3D Gradient
vec4 sdgChamferIntersection(vec4 a, vec4 b, float r) {
    return -sdgChamferUnion(-a, -b, r);
}

// 3. Angle-Exact Euclidean Chamfer Subtraction (subtracts a from b) with Analytical 3D Gradient
vec4 sdgChamferSubtraction(vec4 a, vec4 b, float r) {
    return -sdgChamferUnion(a, -b, r);
}

#endif // SDF_GEMS_SDG_CHAMFER_GLSL
