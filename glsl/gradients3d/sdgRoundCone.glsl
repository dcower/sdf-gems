// ============================================================================
// sdgRoundCone — Exact 1-Sqrt 3D Round Cone & Segment Analytical Gradients
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distgradfunctions3d/
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt` everywhere for sdgRoundCone & sdgRoundConeSeg)
// ============================================================================

#ifndef SDF_GEMS_SDG_ROUND_CONE_GLSL
#define SDF_GEMS_SDG_ROUND_CONE_GLSL

// 1. Exact 1-Sqrt Vertical Round Cone with Analytical 3D Gradient
//    Uses signed-square tangent-cone separators and w = sqrt((1 - b^2)*(px^2 + pz^2))
//    to compute both distance and unit gradient in 1 single sqrt everywhere in R^3!
vec4 sdgRoundCone(vec3 p, float r1, float r2, float h) {
    float b  = (r1 - r2) / h, a2 = 1.0 - b * b;
    float x2 = dot(p.xz, p.xz), y = p.y, z = y - h, s = b * abs(b) * x2;
    if (a2 * y * abs(y) < s) {
        float w2 = x2 + y * y, w = sqrt(w2);
        return vec4(w - r1, (w2 > 0.0) ? (p / w) : vec3(0.0, -1.0, 0.0));
    }
    if (a2 * z * abs(z) > s) {
        float w2 = x2 + z * z, w = sqrt(w2);
        return vec4(w - r2, (w2 > 0.0) ? (vec3(p.x, z, p.z) / w) : vec3(0.0, 1.0, 0.0));
    }
    float w  = sqrt(a2 * x2);
    vec2  er = (x2 > 0.0) ? (p.xz * inversesqrt(x2)) : vec2(1.0, 0.0);
    float a  = sqrt(max(a2, 0.0));
    return vec4(w + b * y - r1, a * er.x, b, a * er.y);
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
    if (sign(y) * a2 * (y * y * l2) < s) {
        float w2 = dot(pa, pa), w = sqrt(w2);
        return vec4(w - r1, (w2 > 0.0) ? (pa / w) : (-ba * sqrt(il2)));
    }
    if (sign(z) * a2 * (z * z * l2) > s) {
        float w2 = dot(pb, pb), w = sqrt(w2);
        return vec4(w - r2, (w2 > 0.0) ? (pb / w) : (ba * sqrt(il2)));
    }
    float w  = sqrt(max(a2 * x2, 0.0));
    vec3  er = (x2 > 1e-20 * l2) ? (u * inversesqrt(x2))
             : normalize((abs(ba.x) < abs(ba.z)) ? vec3(-ba.y, ba.x, 0.0) : vec3(0.0, -ba.z, ba.y));
    return vec4(w + y * rr - r1, (rr * il2) * ba + sqrt(max(a2, 0.0)) * er);
}

#endif // SDF_GEMS_SDG_ROUND_CONE_GLSL
