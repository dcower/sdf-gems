// ============================================================================
// sdgArchimedeanSolid — 0-Loop, 0-Trig, 1-Sqrt Analytical Unit Gradients & Exact
// 3D Euclidean SDFs for Space-Filling Archimedean & Catalan Polyhedra
// (Truncated Octahedron, Cuboctahedron, Rhombic Dodecahedron)
// @author Discovered using Gemini 4 Argon & Gemini 3.8 Flash (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1`, 1-sqrt Voronoi cone decomposition)
// ============================================================================

#ifndef SDF_GEMS_SDG_ARCHIMEDEAN_SOLID_GLSL
#define SDF_GEMS_SDG_ARCHIMEDEAN_SOLID_GLSL

#include "../common/common.glsl"

vec4 sdgTruncatedOctahedron(vec3 p, float s) {
    vec3 sg = vec3(p.x < 0.0 ? -1.0 : 1.0, p.y < 0.0 ? -1.0 : 1.0, p.z < 0.0 ? -1.0 : 1.0);
    vec3 q = abs(p);
    bool s1 = (q.x < q.y); if (s1) q.xy = q.yx;
    bool s2 = (q.y < q.z); if (s2) q.yz = q.zy;
    bool s3 = (q.x < q.y); if (s3) q.xy = q.yx;

    const float invSqrt3 = 0.5773502691896258;
    float dSq  = q.x - s;
    float dHex = (q.x + q.y + q.z - 1.5 * s) * invSqrt3;

    bool inSq   = (q.y + q.z <= 0.5 * s) && (dSq >= 0.0);
    bool inHex  = (max(2.0 * q.x - q.y - q.z, q.x + q.y - 2.0 * q.z) <= 1.5 * s) && (dHex >= 0.0);
    bool inside = (dSq <= 0.0) && (dHex <= 0.0);

    float d;
    vec3 g;
    if (inside || inSq || inHex) {
        bool useSq = inside ? (dSq >= dHex) : inSq;
        d = useSq ? dSq : dHex;
        g = useSq ? vec3(1.0, 0.0, 0.0) : vec3(invSqrt3);
    } else {
        float t1 = clamp(0.5 * (0.5 * s - q.y + q.z), 0.0, 0.25 * s);
        vec3  d1 = vec3(max(q.x - s, 0.0), q.y - (0.5 * s - t1), q.z - t1);
        float t2 = clamp(0.5 * (0.5 * s - q.x + q.y), 0.0, 0.25 * s);
        vec3  d2 = vec3(q.x - (s - t2), q.y - (0.5 * s + t2), q.z);
        vec3  v  = (dot(d1, d1) <= dot(d2, d2)) ? d1 : d2;
        float v2 = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV;
        g = v * invV;
    }

    if (s3) g.xy = g.yx;
    if (s2) g.yz = g.zy;
    if (s1) g.xy = g.yx;
    return vec4(d, sg * g);
}

vec4 sdgCuboctahedron(vec3 p, float s) {
    vec3 sg = vec3(p.x < 0.0 ? -1.0 : 1.0, p.y < 0.0 ? -1.0 : 1.0, p.z < 0.0 ? -1.0 : 1.0);
    vec3 q = abs(p);
    bool s1 = (q.x < q.y); if (s1) q.xy = q.yx;
    bool s2 = (q.y < q.z); if (s2) q.yz = q.zy;
    bool s3 = (q.x < q.y); if (s3) q.xy = q.yx;

    const float invSqrt3 = 0.5773502691896258;
    float dSq  = q.x - s;
    float dTri = (q.x + q.y + q.z - 2.0 * s) * invSqrt3;

    bool inSq   = (q.y + q.z <= s) && (dSq >= 0.0);
    bool inTri  = (2.0 * q.x - q.y - q.z <= s) && (dTri >= 0.0);
    bool inside = (dSq <= 0.0) && (dTri <= 0.0);

    float d;
    vec3 g;
    if (inside || inSq || inTri) {
        bool useSq = inside ? (dSq >= dTri) : inSq;
        d = useSq ? dSq : dTri;
        g = useSq ? vec3(1.0, 0.0, 0.0) : vec3(invSqrt3);
    } else {
        float t = clamp(0.5 * (s - q.y + q.z), 0.0, 0.5 * s);
        vec3  v = vec3(q.x - s, q.y - (s - t), q.z - t);
        float v2 = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV;
        g = v * invV;
    }

    if (s3) g.xy = g.yx;
    if (s2) g.yz = g.zy;
    if (s1) g.xy = g.yx;
    return vec4(d, sg * g);
}

#ifndef SDF_GEMS_HAS_SDG_RHOMBIC_DODECAHEDRON
#define SDF_GEMS_HAS_SDG_RHOMBIC_DODECAHEDRON
vec4 sdgRhombicDodecahedron(vec3 p, float s) {
    vec3 sg = vec3(p.x < 0.0 ? -1.0 : 1.0, p.y < 0.0 ? -1.0 : 1.0, p.z < 0.0 ? -1.0 : 1.0);
    vec3 q = abs(p);
    bool s1 = (q.x < q.y); if (s1) q.xy = q.yx;
    bool s2 = (q.y < q.z); if (s2) q.yz = q.zy;
    bool s3 = (q.x < q.y); if (s3) q.xy = q.yx;

    const float invSqrt2 = 0.7071067811865475;
    float m = q.x + q.y - s;
    float d;
    vec3 g;
    if (m <= 0.0 || (q.x - q.y + 2.0 * q.z <= s)) {
        d = m * invSqrt2;
        g = vec3(invSqrt2, invSqrt2, 0.0);
    } else {
        float t = clamp((s - q.x + q.y + q.z) / 3.0, 0.0, 0.5 * s);
        vec3  v = vec3(q.x - (s - t), q.y - t, q.z - t);
        float v2 = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV;
        g = v * invV;
    }

    if (s3) g.xy = g.yx;
    if (s2) g.yz = g.zy;
    if (s1) g.xy = g.yx;
    return vec4(d, sg * g);
}
#endif // SDF_GEMS_HAS_SDG_RHOMBIC_DODECAHEDRON

#endif // SDF_GEMS_SDG_ARCHIMEDEAN_SOLID_GLSL
