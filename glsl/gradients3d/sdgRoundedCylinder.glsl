// ============================================================================
// sdgRoundedCylinder — 0/1/2-SFU Cap-Voronoi 3D Rounded Cylinder Distance & Gradient
//
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (0 SFU in flat cap cylinder, 1 SFU in cylindrical wall, 2 SFUs in rim)
// ============================================================================

#ifndef SDF_GEMS_SDG_ROUNDED_CYLINDER_GLSL
#define SDF_GEMS_SDG_ROUNDED_CYLINDER_GLSL

// ----------------------------------------------------------------------------
// sdgRoundedCylinder(p, ra, rb, h)
//
// Exact 3D signed distance and analytical unit gradient vec4(d, gx, gy, gz)
// for a vertical cylinder of core radius r = ra - rb, corner rounding radius rb,
// and core half-height h - rb (total radius ra, total half-height h when
// parameterized as in IQ's sdRoundedCylinder(p, ra, rb, h) where core is
// (ra, h) offset by rb, or core (ra, h) + rounding rb).
// Note: Matches sdRoundedCylinder(p, ra, rb, h) in glsl/primitives3d/sdRoundedCylinder.glsl
// where ra is the inner radial scale (total radius ra + rb or 2*ra + rb; here
// w = vec2(length(p.xz) - 2.0*ra + rb, abs(p.y) - h) offset by rb, and we also
// provide sdgRoundedCylinderExact(p, r, rb, h) for core (r, h) + rounding rb).
//
// Cap-Voronoi 0-SFU Theorem:
// Let r_core = 2.0*ra - rb and wy = abs(p.y) - h.
// Inside and above the flat end-caps, the medial Voronoi boundary between the
// flat cap (wy) and the cylindrical wall (wx = length(p.xz) - r_core) is:
//   wx <= min(wy, 0.0)  <==>  dot(p.xz, p.xz) <= (r_core + min(wy, 0.0))^2
// Whenever rc = r_core + min(wy, 0.0) >= 0.0 and dot(p.xz, p.xz) <= rc * rc,
// both the distance wy - rb and the unit gradient (0, sign(p.y), 0) are
// completely independent of length(p.xz), requiring 0 sqrt / 0 inversesqrt!
// ----------------------------------------------------------------------------
vec4 sdgRoundedCylinder( vec3 p, float ra, float rb, float h )
{
    float r  = 2.0 * ra - rb;
    float r2 = dot(p.xz, p.xz);
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float wy = abs(p.y) - h;

    // 1. Flat cap Voronoi region (both exterior above cap AND interior paraboloid): 0 SFU!
    float rc = r + min(wy, 0.0);
    if( rc >= 0.0 && r2 <= rc * rc )
    {
        return vec4( wy - rb, 0.0, sy, 0.0 );
    }

    // Compute radial distance and unit radial normal (1st inversesqrt)
    float invR = (r2 > 0.0) ? inversesqrt(r2) : 0.0;
    float l    = r2 * invR;
    vec2  er   = (r2 > 0.0) ? p.xz * invR : vec2(1.0, 0.0);
    float wx   = l - r;

    // 2. Cylindrical wall Voronoi region (exterior & interior): 1 SFU total!
    if( wy <= 0.0 )
    {
        return vec4( wx - rb, er.x, 0.0, er.y );
    }

    // 3. Exterior toroidal corner rim (wx > 0 && wy > 0): 2nd inversesqrt
    float l2   = wx * wx + wy * wy;
    float invL = inversesqrt(l2);
    return vec4( l2 * invL - rb, wx * er.x * invL, sy * wy * invL, wx * er.y * invL );
}

#endif // SDF_GEMS_SDG_ROUNDED_CYLINDER_GLSL
