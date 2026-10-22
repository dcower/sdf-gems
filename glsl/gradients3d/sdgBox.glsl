// ============================================================================
// sdgBox & sdgRoundBox — Exact 0/1-SFU 3D Box & Rounded Box with Analytical Gradient
// @author Inigo Quilez — upgraded to 0-SFU face-prism & division-free rsq using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distgradfunctions3d/
// @exact  Yes (`||grad d|| = 1`, 0 SFUs & 0 divs in interior + all 6 exterior face prisms)
// ============================================================================

#ifndef SDF_GEMS_SDG_BOX_GLSL
#define SDF_GEMS_SDG_BOX_GLSL

vec4 sdgBox( in vec3 p, in vec3 b )
{
    vec3 w = abs(p) - b;
    vec3 s = vec3(p.x < 0.0 ? -1.0 : 1.0,
                  p.y < 0.0 ? -1.0 : 1.0,
                  p.z < 0.0 ? -1.0 : 1.0);
    if( max(min(w.x, w.y), min(max(w.x, w.y), w.z)) <= 0.0 )
    {
        return vec4( max(w.x, max(w.y, w.z)),
                     s * ((w.x > w.y) ? ((w.x > w.z) ? vec3(1,0,0) : vec3(0,0,1))
                                      : ((w.y > w.z) ? vec3(0,1,0) : vec3(0,0,1))) );
    }
    vec3  q    = max(w, 0.0);
    float l2   = dot(q, q), invL = inversesqrt(l2);
    return vec4( l2 * invL, s * q * invL );
}

vec4 sdgRoundBox( in vec3 p, in vec3 b, in float r )
{
    vec4 res = sdgBox( p, b - r );
    return vec4( res.x - r, res.yzw );
}

#endif // SDF_GEMS_SDG_BOX_GLSL
