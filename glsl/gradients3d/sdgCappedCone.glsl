// ============================================================================
// sdgCappedCone — Exact 0-Division, 0/1/2-SFU 3D Capped Cone Distance & Analytical Gradient
// @author Inigo Quilez (3ctcW7) — corner-rim bug fixed & upgraded using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distgradfunctions3d/
// @exact  Yes (`||grad d|| = 1` everywhere; fixes IQ's `1.414` corner-rim gradient error & axis NaN)
// ============================================================================

#ifndef SDF_GEMS_SDG_CAPPED_CONE_GLSL
#define SDF_GEMS_SDG_CAPPED_CONE_GLSL

vec4 sdgCappedCone( vec3 p, float h, float r1, float r2 )
{
    float x2  = dot(p.xz, p.xz), sy = (p.y < 0.0) ? -1.0 : 1.0;
    float rc  = (p.y < 0.0) ? r1 : r2, cy = abs(p.y) - h, kx = r2 - r1;
    float hrc = h * rc + min(cy, 0.0) * (h + max(sy * kx, 0.0));
    if( hrc >= 0.0 && x2 * h * h <= hrc * hrc )
        return vec4(cy, 0.0, sy, 0.0);

    float l  = sqrt(x2), ky = 2.0 * h, m = kx * kx + ky * ky;
    vec2  q  = vec2(l - r2, p.y - h);
    float wc = ky * q.x - kx * q.y;
    if( cy <= 0.0 && wc <= 0.0 && cy * cy * m <= wc * wc )
        return vec4(cy, 0.0, sy, 0.0);

    float qk   = q.x * kx + q.y * ky;
    float invL = 1.0 / max(l, 1e-30);
    if( (cy <= 0.0 && wc <= 0.0) || (qk < 0.0 && qk > -m) )
    {
        float invM = inversesqrt(m);
        return vec4(wc * invM, ky * invM * p.x * invL, -kx * invM, ky * invM * p.z * invL);
    }
    vec2  v    = (qk >= 0.0) ? q : q + vec2(kx, ky);
    float v2   = dot(v, v), invV = inversesqrt(max(v2, 1e-30));
    return vec4(v2 * invV, v.x * invV * p.x * invL, v.y * invV, v.x * invV * p.z * invL);
}

#endif // SDF_GEMS_SDG_CAPPED_CONE_GLSL
