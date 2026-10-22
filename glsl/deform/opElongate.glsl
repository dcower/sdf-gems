// ============================================================================
// opElongate, opElongate1D & sdgElongateFinalize — Exact 3D Domain Elongation & Analytical Gradient
// @author Inigo Quilez (opElongate) — upgraded with sdgElongateFinalize using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1` in exterior, elongated slabs, and interior core)
// ============================================================================

// 1. 1D/Direct Domain Clamping (exact for 1D elongation along a single axis; zero-gradient core in 2D/3D)
vec3 opElongate1D( in vec3 p, in vec3 h )
{
    return p - clamp(p, -h, h);
}

// 2. 3D Interior-Exact Domain Elongation: returns w = vec4(q_pos, core_d); evaluate as w.w + sdPrimitive(w.xyz)
vec4 opElongate( in vec3 p, in vec3 h )
{
    vec3 q = abs(p) - h;
    return vec4( max(q, 0.0), min(max(q.x, max(q.y, q.z)), 0.0) );
}

// 3. 3D Analytical Distance + Gradient Elongation Finalizer for dg = sdgPrimitive(max(abs(p) - h, 0.0)):
//    Propagates exact unit gradients through the exterior, clamped slabs (masking coordinate-plane creases),
//    and the interior core box B_h (where q_+ = 0 is constant and grad d = sign(p) * e_argmax(q)).
vec4 sdgElongateFinalize( in vec3 p, in vec3 h, in vec4 dg )
{
    vec3  s = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    vec3  q = abs(p) - h;
    float g = max(q.x, max(q.y, q.z));
    vec3  c = (q.x >= max(q.y, q.z)) ? vec3(1.0, 0.0, 0.0)
            : ((q.y >= q.z)          ? vec3(0.0, 1.0, 0.0) : vec3(0.0, 0.0, 1.0));
    if( g <= 0.0 ) return vec4(dg.x + g, s * c);
    vec3  m  = vec3((q.x > 0.0) ? dg.y : 0.0, (q.y > 0.0) ? dg.z : 0.0, (q.z > 0.0) ? dg.w : 0.0);
    float l2 = dot(m, m);
    return vec4(dg.x, s * ((l2 > 0.0) ? (m * inversesqrt(l2)) : c));
}
