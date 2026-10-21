// ============================================================================
// sdgRoundCone & sdgTorus — Exact 1-Sqrt 3D Round Cone & 2-Sqrt Torus Gradients
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distgradfunctions3d/
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt` everywhere for sdgRoundCone/Seg, 2 `sqrt`s for sdgTorus)
// ============================================================================

#ifndef SDF_GEMS_SDG_ROUND_CONE_GLSL
#define SDF_GEMS_SDG_ROUND_CONE_GLSL

// 1. Exact 1-Sqrt Vertical Round Cone with Analytical 3D Gradient
//    Uses signed-square tangent-cone separators and w = sqrt((1 - b^2)*(px^2 + pz^2))
//    to compute both distance and unit gradient in 1 single sqrt everywhere in R^3!
vec4 sdgRoundCone(vec3 p, float r1, float r2, float h) {
    float b  = (r1 - r2) / h, a2 = 1.0 - b * b;
    float x2 = dot(p.xz, p.xz), y = p.y, z = y - h, s = b * abs(b) * x2;
    if (a2 * y * abs(y) < s) { float w = sqrt(x2 + y * y); return vec4(w - r1, p / max(w, 1e-12)); }
    if (a2 * z * abs(z) > s) { float w = sqrt(x2 + z * z); return vec4(w - r2, vec3(p.x, z, p.z) / max(w, 1e-12)); }
    float w = max(sqrt(a2 * x2), 1e-12);
    return vec4(w + b * y - r1, a2 * p.x / w, b, a2 * p.z / w);
}

// 2. Exact 1-Sqrt O(L^2) Scale-Safe Arbitrary-Segment Round Cone with Analytical 3D Gradient
//    Eliminates IQ's O(L^8) FP16 overflow and uses 1 sqrt everywhere in R^3.
vec4 sdgRoundConeSeg(vec3 p, vec3 a, vec3 b, float r1, float r2) {
    vec3  ba = b - a, pa = p - a, pb = p - b;
    float l2 = dot(ba, ba), il2 = 1.0 / l2;
    float rr = r1 - r2, b2 = rr * rr * il2, a2 = 1.0 - b2;
    float y  = dot(pa, ba) * il2, z = y - 1.0;
    vec3  u  = pa - ba * y;
    float x2 = dot(u, u), s = sign(rr) * b2 * x2;
    if (sign(y) * a2 * (y * y * l2) < s) { float w = max(length(pa), 1e-12); return vec4(w - r1, pa / w); }
    if (sign(z) * a2 * (z * z * l2) > s) { float w = max(length(pb), 1e-12); return vec4(w - r2, pb / w); }
    float w = max(sqrt(a2 * x2), 1e-12);
    return vec4(w + y * rr - r1, (rr * il2) * ba + (a2 / w) * u);
}

// 3. Exact 2-Sqrt 3D Torus with Analytical Gradient (down from 3 sqrts in IQ's sdgTorus
//    by factoring ||(px*(h - ra), py*h, pz*(h - ra))||_2 = h * l).
vec4 sdgTorus(vec3 p, float ra, float rb) {
    float h = max(length(p.xz), 1e-12);
    float u = h - ra;
    float l = max(length(vec2(u, p.y)), 1e-12);
    return vec4(l - rb, vec3(p.x * u / h, p.y, p.z * u / h) / l);
}

#endif // SDF_GEMS_SDG_ROUND_CONE_GLSL
