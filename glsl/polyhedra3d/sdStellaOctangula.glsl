// ============================================================================
// sdStellaOctangula & sdgStellaOctangula — Exact 0/1-Sqrt 3D Stella Octangula
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in face strip, 1 sqrt in ridge/valley)
// @novelty Exact 0/1-sqrt Euclidean SDF and analytical gradient for Kepler's non-convex
//          Stella Octangula (compound of two dual tetrahedra, 24 visible outer triangles,
//          36 edges, 14 vertices) via O_h Weyl chamber folding, where T_- subset T_+ collapses
//          the entire non-convex solid to a single half-space x + y - z <= s with 1 outer
//          ridge and 1 unclamped inward valley seam.
// ============================================================================

#ifndef SDF_GEMS_SD_STELLA_OCTANGULA_GLSL
#define SDF_GEMS_SD_STELLA_OCTANGULA_GLSL

float sdStellaOctangula( vec3 p, float s )
{
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;

    float m = p.x + p.y - p.z - s;
    if( m >= -3.0 * p.z && m <= 1.5 * (p.y - p.z) )
        return m * 0.57735027; // 0 sqrts in interior & exterior face prism!

    float t = min(0.5 * (p.y + p.z), s);
    float h = 0.5 * (p.x + p.y - s);
    vec3  v = (m < 0.0) ? vec3(h, h, p.z) : vec3(p.x - s, p.y - t, p.z - t);
    return ((m < 0.0) ? -1.0 : 1.0) * length(v);
}

vec4 sdgStellaOctangula( vec3 p, float s )
{
    vec3 sg = vec3((p.x < 0.0) ? -1.0 : 1.0,
                   (p.y < 0.0) ? -1.0 : 1.0,
                   (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    bool s1 = (p.x < p.y); if( s1 ) p.xy = p.yx;
    bool s2 = (p.y < p.z); if( s2 ) p.yz = p.zy;
    bool s3 = (p.x < p.y); if( s3 ) p.xy = p.yx;

    float m = p.x + p.y - p.z - s;
    float d;
    vec3  g;
    if( m >= -3.0 * p.z && m <= 1.5 * (p.y - p.z) )
    {
        d = m * 0.57735027;
        g = vec3(0.57735027, 0.57735027, -0.57735027);
    }
    else
    {
        float t    = min(0.5 * (p.y + p.z), s);
        float h    = 0.5 * (p.x + p.y - s);
        vec3  v    = (m < 0.0) ? vec3(h, h, p.z) : vec3(p.x - s, p.y - t, p.z - t);
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        float sm   = (m < 0.0) ? -1.0 : 1.0;
        d = sm * v2 * invV;
        g = (sm * invV) * v;
    }

    if( s3 ) g.xy = g.yx;
    if( s2 ) g.yz = g.zy;
    if( s1 ) g.xy = g.yx;
    return vec4(d, sg * g);
}

#endif // SDF_GEMS_SD_STELLA_OCTANGULA_GLSL
