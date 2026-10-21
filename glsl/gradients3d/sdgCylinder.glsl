// ============================================================================
// sdgCylinder — 0/1/2-Sqrt Vertical & 1/2-Sqrt Segment Cylinders with 3D Gradients
// @author Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distgradfunctions3d/
// @exact  Yes (`||grad d|| = 1`, fixes 100% exterior axis `0.0 * NaN` bug in IQ's `3fcyz2`)
// ============================================================================

// 1. Vertical Capped Cylinder (radius r, half-height h) via Unified Cap-Voronoi Radius Theorem:
//    - 0 sqrts in the entire exterior & interior flat-cap Voronoi zones (and on the exterior y-axis!)
//    - 1 sqrt in the interior & exterior cylindrical wall Voronoi zones
//    - 2 sqrts only in the exterior circular rim wedge (wy > 0 && r2 > r^2)
vec4 sdgCylinder( vec3 p, float r, float h )
{
    float r2 = dot(p.xz, p.xz), wy = abs(p.y) - h, sy = (p.y < 0.0) ? -1.0 : 1.0;
    float rc = r + min(wy, 0.0);
    if( rc >= 0.0 && r2 <= rc * rc ) return vec4(wy, 0.0, sy, 0.0);
    float l = sqrt(r2), wx = l - r;
    vec2  dir = (l > 0.0) ? (p.xz / l) : vec2(1.0, 0.0);
    if( wy <= 0.0 ) return vec4(wx, dir.x, 0.0, dir.y);
    float d = sqrt(wx * wx + wy * wy);
    return vec4(d, vec3(wx * dir.x, sy * wy, wx * dir.y) / d);
}

// 2. Arbitrary-Segment Capped Cylinder via L-Cancellation Identity:
//    wy = uy * L and dv = sy * ba / L  ==>  wy * dv = (sy * uy) * ba (L cancels identically!).
//    Requires 1 sqrt in exterior cap & wall zones, 2 sqrts in rim & interior cap zones (vs 3 in IQ).
vec4 sdgCylinderSeg( vec3 p, vec3 a, vec3 b, float r )
{
    vec3  ba = b - a, pa = p - a;
    float l2 = dot(ba, ba), t = dot(pa, ba) / l2;
    vec3  u  = pa - t * ba;
    float x2 = dot(u, u), s = t - 0.5, sy = (s < 0.0) ? -1.0 : 1.0, uy = abs(s) - 0.5;
    if( uy >= 0.0 ) {
        if( x2 <= r * r ) {
            float invL = inversesqrt(l2);
            return vec4((uy * l2) * invL, (sy * invL) * ba);
        }
        float x = sqrt(x2), wx = x - r, d = sqrt(wx * wx + uy * uy * l2);
        return vec4(d, ((wx / x) * u + (sy * uy) * ba) / d);
    }
    if( x2 > r * r ) {
        float x = sqrt(x2);
        return vec4(x - r, u / x);
    }
    float x = sqrt(x2), wx = x - r, wy2 = uy * uy * l2;
    if( wx * wx <= wy2 ) {
        vec3 du = (x > 0.0) ? (u / x) : normalize(vec3(ba.y - ba.z, ba.z - ba.x, ba.x - ba.y));
        return vec4(wx, du);
    }
    float invL = inversesqrt(l2);
    return vec4((uy * l2) * invL, (sy * invL) * ba);
}
