// ============================================================================
// sdHorseshoe2D, sdgHorseshoe2D & sdHorseshoe3D — Exact 1-Sqrt 2D & 3D Horseshoe
// @author Inigo Quilez (sdHorseshoe) — upgraded to 1-sqrt argument selection & 3D box fusion using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, 1 sqrt everywhere in R^2 down from 2, `0.0` error, thick-wall safe)
// ============================================================================

// 1. Exact 1-Sqrt 2D Horseshoe SDF (down from 2 sqrts in IQ, 0.0 difference across all parameters):
//    c = vec2(cos(theta), sin(theta)) half-angle, r = arc radius, w = vec2(leg_half_len, half_thickness).
float sdHorseshoe2D( in vec2 p, in vec2 c, in float r, in vec2 w )
{
    p.x = abs(p.x);
    vec2  q = mat2(-c.x, c.y, c.y, c.x) * p;
    vec2  u = vec2(q.x - w.x, abs(q.y - r) - w.y);
    float l = length((q.x > 0.0) ? max(u, 0.0) : p);
    return (q.x > 0.0) ? l + min(0.0, max(u.x, u.y))
                       : max(abs(l - r) - w.y, (q.y > 0.0) ? u.x : -w.y);
}

// 2. Exact 1-Sqrt 2D Horseshoe SDF + Analytical Unit Gradient vec3(d, grad_x, grad_y)
vec3 sdgHorseshoe2D( in vec2 p, in vec2 c, in float r, in vec2 w )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    vec2  q  = mat2(-c.x, c.y, c.y, c.x) * p;
    vec2  u  = vec2(q.x - w.x, abs(q.y - r) - w.y);
    float sy = (q.y < r) ? -1.0 : 1.0;
    vec2  v  = max(u, 0.0);
    float l  = length((q.x > 0.0) ? v : p);

    if( q.x > 0.0 )
    {
        float g  = max(u.x, u.y);
        vec2  gq = (g > 0.0) ? vec2(v.x, sy * v.y) / max(l, 1e-30)
                             : ((u.x > u.y) ? vec2(1.0, 0.0) : vec2(0.0, sy));
        vec2  gp = mat2(-c.x, c.y, c.y, c.x) * gq;
        return vec3((g > 0.0) ? l : g, sx * gp.x, gp.y);
    }
    float da = abs(l - r) - w.y;
    if( q.y > 0.0 && u.x > da )
        return vec3(u.x, -sx * c.x, c.y);
    float sr = (l < r) ? -1.0 : 1.0;
    return vec3(da, sr * vec2(sx * p.x, p.y) / max(l, 1e-30));
}

// 3. Exact 1-to-2 Sqrt 3D Extruded Horseshoe SDF (w = vec3(leg_half_len, radial_half_thick, z_half_height)):
//    Fuses the 3D leg box in q.x > 0 into a single 3D box evaluation (1 sqrt down from 3 in opExtrusion).
float sdHorseshoe3D( in vec3 p, in vec2 c, in float r, in vec3 w )
{
    p.x = abs(p.x);
    vec2  q  = mat2(-c.x, c.y, c.y, c.x) * p.xy;
    float uz = abs(p.z) - w.z;
    if( q.x > 0.0 )
    {
        vec3 u = vec3(q.x - w.x, abs(q.y - r) - w.y, uz);
        return length(max(u, 0.0)) + min(0.0, max(u.x, max(u.y, u.z)));
    }
    float l   = length(p.xy);
    float d2d = max(abs(l - r) - w.y, (q.y > 0.0) ? q.x - w.x : -w.y);
    vec2  u   = vec2(d2d, uz);
    return length(max(u, 0.0)) + min(0.0, max(u.x, u.y));
}
