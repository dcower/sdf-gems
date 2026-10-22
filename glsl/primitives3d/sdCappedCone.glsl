// ============================================================================
// sdCappedCone — Exact 0-Division, 0/1/2-Sqrt 3D Capped Cone (Cap-Voronoi & Cross-Product Wall)
// @author Inigo Quilez — upgraded to 0-div, 0-clamp, 0/1/2-sqrt using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in flat-cap Voronoi cylinder/cone, `1` coord sqrt on slanted wall, `0` divs)
// ============================================================================

#ifndef SDF_GEMS_SD_CAPPED_CONE_GLSL
#define SDF_GEMS_SD_CAPPED_CONE_GLSL

float sdCappedCone( vec3 p, float h, float r1, float r2 )
{
    float x2  = dot(p.xz, p.xz), sy = (p.y < 0.0) ? -1.0 : 1.0;
    float rc  = (p.y < 0.0) ? r1 : r2, cy = abs(p.y) - h, kx = r2 - r1;
    float hrc = h * rc + min(cy, 0.0) * (h + max(sy * kx, 0.0));
    if( hrc >= 0.0 && x2 * h * h <= hrc * hrc ) return cy;

    float l  = sqrt(x2), ky = 2.0 * h, m = kx * kx + ky * ky;
    vec2  q  = vec2(l - r2, p.y - h);
    float wc = ky * q.x - kx * q.y;
    if( cy <= 0.0 && wc <= 0.0 )
        return (cy * cy * m <= wc * wc) ? cy : wc * inversesqrt(m);
    float qk = q.x * kx + q.y * ky;
    if( qk < 0.0 && qk > -m ) return wc * inversesqrt(m);
    return length((qk >= 0.0) ? q : q + vec2(kx, ky));
}

#endif // SDF_GEMS_SD_CAPPED_CONE_GLSL
