// ============================================================================
// sdParallelogram2D & sdParallelogramPrism — Exact 1-Sqrt 2D Parallelogram & 3D Prism
// @author Inigo Quilez (2D) — 1-sqrt 3D prism extrusion discovered using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`1` `sqrt` everywhere in R^2 and R^3 across all skews sk < 0, sk == 0, sk > 0)
// ============================================================================

#ifndef SDF_GEMS_SD_PARALLELOGRAM2D_GLSL
#define SDF_GEMS_SD_PARALLELOGRAM2D_GLSL

// 1. Exact 1-sqrt 2D Parallelogram of width wi, half-height he, and horizontal skew sk
float sdParallelogram2D(vec2 p, float wi, float he, float sk) {
    vec2  e = vec2(sk, he);
    p = (p.y < 0.0) ? -p : p;
    vec2  w = p - e; w.x -= clamp(w.x, -wi, wi);
    vec2  d = vec2(dot(w, w), -w.y);
    float s = p.x * e.y - p.y * e.x;
    p = (s < 0.0) ? -p : p;
    vec2  v = p - vec2(wi, 0.0); v -= e * clamp(dot(v, e) / dot(e, e), -1.0, 1.0);
    d = min(d, vec2(dot(v, v), wi * he - abs(s)));
    return sqrt(d.x) * sign(-d.y);
}

// 2. Exact 1-sqrt 3D Extruded Parallelogram Prism (Universal 1-Sqrt Extrusion Theorem:
//    fuses the 2D boundary sqrt and 3D cap-corner sqrt into a single sqrt everywhere in R^3)
float sdParallelogramPrism(vec3 p, float wi, float he, float sk, float hz) {
    vec2  e = vec2(sk, he);
    vec2  q = (p.y < 0.0) ? -p.xy : p.xy;
    vec2  w = q - e; w.x -= clamp(w.x, -wi, wi);
    vec2  d = vec2(dot(w, w), -w.y);
    float s = q.x * e.y - q.y * e.x;
    q = (s < 0.0) ? -q : q;
    vec2  v = q - vec2(wi, 0.0); v -= e * clamp(dot(v, e) / dot(e, e), -1.0, 1.0);
    d = min(d, vec2(dot(v, v), wi * he - abs(s)));
    float qz = abs(p.z) - hz, m = max(qz, 0.0);
    float r = sqrt((d.y < 0.0) ? d.x + m * m : d.x);
    return (d.y < 0.0) ? r : max(-r, qz);
}

#endif // SDF_GEMS_SD_PARALLELOGRAM2D_GLSL
