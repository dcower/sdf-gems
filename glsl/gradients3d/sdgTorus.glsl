// ============================================================================
// sdgTorus — Singularity-Free All-Regime 3D Torus Distance + Analytical Gradient
// @author Inigo Quilez — upgraded to All-Regime & NaN-Free using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distgradfunctions3d/
// @exact  Yes (`||grad d|| = 1` everywhere including symmetry axis & core circle, 0 NaNs)
// ============================================================================

#ifndef SDF_GEMS_SDG_TORUS_GLSL
#define SDF_GEMS_SDG_TORUS_GLSL

// Returns vec4(d, grad.xyz) for a 3D Torus of major radius ra and minor radius rb
// across all geometric regimes (Ring ra > rb, Horn ra = rb, and Apple Torus ra < rb):
//   - Eliminates IQ's 3rd `sqrt` (`normalize(...)`) via factored radial unit vector er.
//   - Eliminates `0.0 / 0.0 = NaN` on both the symmetry axis (p.xz = 0) and core circle (u = 0, p.y = 0).
//   - Evaluates the exact interior dimple-cusp cone for self-intersecting Apple Tori (ra < rb)
//     in 1 per-pixel `inversesqrt`!
vec4 sdgTorus( vec3 p, float ra, float rb )
{
    float x2 = dot(p.xz, p.xz), ay = abs(p.y), sy = (p.y < 0.0) ? -1.0 : 1.0;
    float b2 = rb * rb - ra * ra;
    if( b2 > 0.0 )
    {
        float b = sqrt(b2), dy = ay - b;
        if( dy < 0.0 && b2 * x2 < ra * ra * dy * dy )
        {
            float l2   = x2 + dy * dy;
            float invL = inversesqrt(max(l2, 1e-30));
            return vec4(-l2 * invL, (l2 > 0.0) ? (-vec3(p.x, sy * dy, p.z) * invL) : vec3(0.0, sy, 0.0));
        }
    }
    float invH = inversesqrt(max(x2, 1e-30));
    float h    = (x2 > 0.0) ? (x2 * invH) : 0.0;
    vec2  er   = (x2 > 0.0) ? (p.xz * invH) : vec2(1.0, 0.0);
    float u    = h - ra, l2 = u * u + p.y * p.y;
    float invL = inversesqrt(max(l2, 1e-30));
    return vec4(l2 * invL - rb, (l2 > 0.0) ? (vec3(u * er.x, p.y, u * er.y) * invL) : vec3(er.x, 0.0, er.y));
}

#endif // SDF_GEMS_SDG_TORUS_GLSL
