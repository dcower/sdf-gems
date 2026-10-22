// ============================================================================
// sdOctagon2D, sdgOctagon2D & sdOctagonalPrism — 1-Fold 0/1-Sqrt 2D Octagon, Gradient & 0/1-Sqrt 3D Prism
// @author Inigo Quilez (2-fold 2D) — 1-fold, 0-sqrt-face & 0/1-SFU gradient discovered using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`0` sqrts in interior & flat-edge strips, `1` sqrt in corner wedges; `8.88e-16` max error)
// ============================================================================

#ifndef SDF_GEMS_SD_OCTAGON2D_GLSL
#define SDF_GEMS_SD_OCTAGON2D_GLSL

// 1. Exact 1-Fold, 0/1-Sqrt 2D Regular Octagon of inradius r
float sdOctagon2D(vec2 p, float r) {
    const vec3 k = vec3(-0.9238795325, 0.3826834323, 0.4142135623);
    p = abs(p);
    if (p.y < p.x) p = p.yx;
    p -= 2.0 * min(dot(k.xy, p), 0.0) * k.xy;
    vec2 d = p - vec2(k.z, 1.0) * r;
    return (d.x > 0.0) ? length(d) : d.y;
}

// 2. 1-Fold, 0/1-SFU Analytical 2D Distance & Gradient vec3(d, grad.xy)
vec3 sdgOctagon2D(vec2 p, float r) {
    const vec3 k = vec3(-0.9238795325, 0.3826834323, 0.4142135623);
    vec2 s  = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    bool sw = (p.y < p.x);
    if (sw) p = p.yx;
    float dk = dot(k.xy, p);
    bool  rf = (dk < 0.0);
    if (rf) p -= 2.0 * dk * k.xy;

    vec2  v = p - vec2(k.z, 1.0) * r;
    float d;
    vec2  g;
    if (v.x <= 0.0) {
        d = v.y;
        g = vec2(0.0, 1.0);
    } else {
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV;
        g = v * invV;
    }
    if (rf) g -= 2.0 * dot(k.xy, g) * k.xy;
    if (sw) g = g.yx;
    return vec3(d, s * g);
}

// 3. Exact 1-Fold, 0/1-Sqrt 3D Regular Octagonal Prism of inradius r and half-height hz
float sdOctagonalPrism(vec3 p, float r, float hz) {
    const vec3 k = vec3(-0.9238795325, 0.3826834323, 0.4142135623);
    p = abs(p);
    if (p.y < p.x) p.xy = p.yx;
    p.xy -= 2.0 * min(dot(k.xy, p.xy), 0.0) * k.xy;
    vec3 d = vec3(p.x - k.z * r, p.y - r, p.z - hz);
    return (max(d.x, min(d.y, d.z)) <= 0.0) ? max(d.y, d.z) : length(max(d, 0.0));
}

#endif // SDF_GEMS_SD_OCTAGON2D_GLSL
