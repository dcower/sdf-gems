// ============================================================================
// sdReuleauxTetrahedron, sdgReuleauxTetrahedron, sdReuleauxRotor3D & sdgReuleauxRotor3D
// — Exact 1/2-Sqrt 3D Reuleaux Tetrahedron & Constant-Width Reuleaux Rotor
//
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (0 trig, 0 divs, 1 sqrt in face & vertex wedges, 2 sqrts in arc edge)
// ============================================================================

#ifndef SDF_GEMS_SD_REULEAUX_3D_GLSL
#define SDF_GEMS_SD_REULEAUX_3D_GLSL

// ----------------------------------------------------------------------------
// sdReuleauxTetrahedron(p, w) — Exact 3D Reuleaux Tetrahedron (Meissner Core)
//
// The 3D Reuleaux tetrahedron of sphere radius / edge length w is the
// intersection of 4 spheres of radius w centered at the 4 vertices of a
// regular tetrahedron V0=(s,s,s), V1=(s,-s,-s), V2=(-s,s,-s), V3=(-s,-s,s)
// where s = w / (2 * sqrt(2)) = 0.35355339 * w.
//
// Td Vertex-Distance Total-Ordering Theorem:
// Folding p into the 24-fold Td parity-sorted Weyl chamber x >= y >= |z|
// guarantees x+y+z >= x-y-z >= -x+y-z >= -x-y+z, so:
//   ||p - V3|| >= ||p - V2|| >= ||p - V1|| >= ||p - V0||
// Thus all 4 spherical faces collapse to a SINGLE sphere centered at V3=(-s,-s,s),
// all 4 vertices collapse to V0=(s,s,s), and all 6 circular-arc edges collapse
// to the single arc in the plane y = z centered at (-s, 0, 0) with radius
// sqrt(3)/2 * w = sqrt(6) * s!
// ----------------------------------------------------------------------------
float sdReuleauxTetrahedron( vec3 p, float w )
{
    float s   = 0.35355339 * w; // w / (2.0 * sqrt(2.0))
    float sig = (p.x * p.y * p.z < 0.0) ? -1.0 : 1.0;
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;
    p.z *= sig;

    vec3  v       = p - vec3(-s, -s, s);
    float v2      = dot(v, v);
    float dyz     = p.y - p.z, d2s = dyz + 2.0 * s;
    bool  in_face = (2.0 * d2s * d2s >= v2);
    if( in_face || p.y + p.z >= p.x + s )
        return length(p - (in_face ? vec3(-s, -s, s) : vec3(s))) - (in_face ? w : 0.0);

    float ux  = p.x + s, syz = p.y + p.z;
    float dr  = sqrt(ux * ux + 0.5 * syz * syz) - 2.44948974 * s;
    return sqrt(dr * dr + 0.5 * dyz * dyz);
}

// ----------------------------------------------------------------------------
// sdgReuleauxTetrahedron(p, w) — Exact 1/2-SFU 3D Distance + Unit Gradient
// ----------------------------------------------------------------------------
vec4 sdgReuleauxTetrahedron( vec3 p, float w )
{
    float s  = 0.35355339 * w;
    vec3  sg = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    float sig = sg.x * sg.y * sg.z;
    p = abs(p);
    bool s1 = (p.x < p.y); if( s1 ) p.xy = p.yx;
    bool s2 = (p.y < p.z); if( s2 ) p.yz = p.zy;
    bool s3 = (p.x < p.y); if( s3 ) p.xy = p.yx;
    p.z *= sig;

    vec3  vf      = p - vec3(-s, -s, s);
    float vf2     = dot(vf, vf);
    float dyz     = p.y - p.z, d2s = dyz + 2.0 * s;
    bool  in_face = (2.0 * d2s * d2s >= vf2);
    float d;
    vec3  g;
    if( in_face || p.y + p.z >= p.x + s )
    {
        vec3  v    = in_face ? vf : (p - vec3(s));
        float v2   = in_face ? vf2 : dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV - (in_face ? w : 0.0);
        g = v * invV;
    }
    else
    {
        float ux   = p.x + s, syz = p.y + p.z;
        float r2   = ux * ux + 0.5 * syz * syz;
        float invR = inversesqrt(max(r2, 1e-30));
        float dr   = r2 * invR - 2.44948974 * s;
        float d2   = dr * dr + 0.5 * dyz * dyz;
        float invD = inversesqrt(max(d2, 1e-30));
        float sr   = dr * invR;
        d = d2 * invD;
        g = vec3(sr * ux, 0.5 * (sr * syz + dyz), 0.5 * (sr * syz - dyz)) * invD;
    }

    g.z *= sig;
    if( s3 ) g.xy = g.yx;
    if( s2 ) g.yz = g.zy;
    if( s1 ) g.xy = g.yx;
    return vec4(d, sg * g);
}

// ----------------------------------------------------------------------------
// sdReuleauxRotor3D(p, w) — Exact 1/2-Sqrt 3D Constant-Width Reuleaux Rotor
//
// Solid of revolution of a 2D Reuleaux triangle of constant width w around its
// vertical symmetry axis (the 3D constant-width "acorn" body).
// By the Axis-Center Radial Collapse Theorem, both the top apex cone and the
// entire bottom spherical cap share the on-axis center (0, w / sqrt(3), 0),
// evaluating in a SINGLE 3D sqrt without computing length(p.xz)!
// ----------------------------------------------------------------------------
float sdReuleauxRotor3D( vec3 p, float w )
{
    const float invK = 0.57735027;
    float r2  = dot(p.xz, p.xz), y = p.y;
    float y0  = w * invK, dy0 = y - y0;
    if( 3.0 * r2 <= dy0 * dy0 && (dy0 >= 0.0 || (y <= 0.0 && r2 <= 3.0 * y * y)) )
        return sqrt(r2 + dy0 * dy0) - ((dy0 >= 0.0) ? 0.0 : w);
    bool rim = (y <= -0.5 * y0);
    return length(vec2(sqrt(r2) - (rim ? 0.5 * w : -0.5 * w), y + 0.5 * y0)) - (rim ? 0.0 : w);
}

// ----------------------------------------------------------------------------
// sdgReuleauxRotor3D(p, w) — Exact 1/2-SFU 3D Distance + Unit Gradient
// ----------------------------------------------------------------------------
vec4 sdgReuleauxRotor3D( vec3 p, float w )
{
    const float invK = 0.57735027;
    float r2  = dot(p.xz, p.xz), y = p.y;
    float y0  = w * invK, dy0 = y - y0;
    if( 3.0 * r2 <= dy0 * dy0 && (dy0 >= 0.0 || (y <= 0.0 && r2 <= 3.0 * y * y)) )
    {
        float l2   = r2 + dy0 * dy0;
        float invL = inversesqrt(max(l2, 1e-30));
        return vec4(l2 * invL - ((dy0 >= 0.0) ? 0.0 : w), vec3(p.x, dy0, p.z) * invL);
    }
    float invR = inversesqrt(max(r2, 1e-30));
    float rho  = (r2 > 0.0) ? (r2 * invR) : 0.0;
    vec2  er   = (r2 > 0.0) ? (p.xz * invR) : vec2(1.0, 0.0);
    bool  rim  = (y <= -0.5 * y0);
    vec2  u    = vec2(rho - (rim ? 0.5 * w : -0.5 * w), y + 0.5 * y0);
    float l2   = dot(u, u);
    float invL = inversesqrt(max(l2, 1e-30));
    return vec4(l2 * invL - (rim ? 0.0 : w), vec3(u.x * er.x, u.y, u.x * er.y) * invL);
}

#endif // SDF_GEMS_SD_REULEAUX_3D_GLSL
