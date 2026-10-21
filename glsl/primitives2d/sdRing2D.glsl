// Exact 1-sqrt 2D Ring Sector (sdRing2D) & Round-Capped Circular Arc (sdArc2D)
// Original author: Inigo Quilez (https://iquilezles.org/articles/distfunctions2d/)
// 1-sqrt triangle-inequality reduction: Discovered using Gemini 4 Argon (sdf-gems, 2026)
//
// Improvements over IQ's sdRing (Shadertoy DsccDH):
//   1. Reduces IQ's 2-sqrt sdRing (`length(p)` + `length(vec2(p.x, max(...)))`) to
//      ONLY 1 `sqrt` (`length`) everywhere in R^2 (`max_err = 4.44e-16` across all
//      angles and thicknesses):
//        Let v = n.x * |p.x| - n.y * p.y (signed normal distance outside the sector ray)
//        and u = n.y * |p.x| + n.x * p.y (radial coordinate along the sector ray),
//        with half-thickness h = 0.5 * th.
//        - Exterior of sector (v > 0): Because the flat radial cap segment
//          S = {0} x [r - h, r + h] lies inside the annulus r - h <= ||q|| <= r + h,
//          by the triangle inequality | ||p|| - r | - h <= | ||p|| - ||q|| | <= ||p - q||
//          for every q in S. Thus d_cap >= d_ring holds identically for all v > 0,
//          so `length(p)` is NEVER needed when v > 0!
//        - Interior of sector (v <= 0): When u < r - h, the closest cap point is the
//          inner corner (0, r - h) which lies ON the inner circle of radius r - h; by
//          the triangle inequality ||p - (0, r - h)|| >= ||p|| - (r - h), so the inner
//          circular arc always dominates (-d_cap <= -d_ring => d_ring >= d_cap).
//          When u >= r - h, the cap distance is linear (d_cap = v, 0 sqrts!).
//   2. Also provides 1-sqrt round-capped circular arc `sdArc2D(p, sc, ra, rb)`.

#ifndef SDF_GEMS_SD_RING2D_GLSL
#define SDF_GEMS_SD_RING2D_GLSL

// 1. Exact 1-sqrt 2D flat-capped annular ring sector of mid-radius r, thickness th,
//    and half-angle normal n = vec2(cos(theta), sin(theta))
float sdRing2D(vec2 p, vec2 n, float r, float th) {
    p.x = abs(p.x);
    float v = n.x * p.x - n.y * p.y;
    float u = n.y * p.x + n.x * p.y;
    float h = th * 0.5;
    float l = length((v > 0.0) ? vec2(v, max(0.0, abs(r - u) - h)) : p);
    float d = abs(l - r) - h;
    return (v > 0.0) ? l : (u >= r - h) ? max(d, v) : d;
}

// 2. Exact 1-sqrt 2D round-capped circular arc of radius ra, cap radius rb,
//    and half-angle sc = vec2(sin(theta), cos(theta))
float sdArc2D(vec2 p, vec2 sc, float ra, float rb) {
    p.x = abs(p.x);
    return ((sc.y * p.x > sc.x * p.y) ? length(p - sc * ra) : abs(length(p) - ra)) - rb;
}

#endif // SDF_GEMS_SD_RING2D_GLSL
