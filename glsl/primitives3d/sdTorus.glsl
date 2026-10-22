// ============================================================================
// sdTorus — Exact All-Regime 3D Torus (Ring R > r, Horn R = r, & Apple R < r)
// @author Inigo Quilez — upgraded to All-Regime Apple Torus (R < r) using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1` for all R > 0, r > 0; 1 per-pixel sqrt in polar cusp cone!)
// ============================================================================

#ifndef SDF_GEMS_SD_TORUS_GLSL
#define SDF_GEMS_SD_TORUS_GLSL

// Exact 3D Euclidean SDF of a Torus of major radius t.x = R and minor radius t.y = r:
//   - Ring Torus (R > r) & Horn Torus (R = r): standard 2-sqrt radial revolution.
//   - Self-Intersecting Apple Torus (R < r): eliminates the internal lemon/spindle
//     self-intersection sheet inside the polar dimple-cusp cone
//       |y| < b  and  (r^2 - R^2) (p_x^2 + p_z^2) < R^2 (b - |y|)^2,
//     where b = sqrt(r^2 - R^2) is the polar dimple height, and evaluates the exact
//     interior distance -sqrt(p_x^2 + p_z^2 + (b - |p_y|)^2) in a single per-pixel sqrt!
float sdTorus( vec3 p, vec2 t )
{
    float x2 = dot(p.xz, p.xz), ay = abs(p.y), b2 = t.y * t.y - t.x * t.x;
    if( b2 > 0.0 )
    {
        float b = sqrt(b2), dy = b - ay;
        if( dy > 0.0 && b2 * x2 < t.x * t.x * dy * dy )
            return -sqrt(x2 + dy * dy);
    }
    return length(vec2(sqrt(x2) - t.x, p.y)) - t.y;
}

#endif // SDF_GEMS_SD_TORUS_GLSL
