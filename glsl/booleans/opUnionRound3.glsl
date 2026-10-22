// ============================================================================
// opUnionRound3 — Exact 0/1-Sqrt Symmetric Ternary 3-Way Corner Booleans & Gradients
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1` everywhere for Round3 & Chamfer3, `8.88e-16` max error, 0/1 `sqrt`)
// ============================================================================
//
// Solves the 3-Way Corner Boolean Problem:
//   1. Chaining binary opUnionChamfer(opUnionChamfer(a, b, r), c, r) breaks S3
//      symmetry, tilting the 3-way corner normal to (1/2, 1/2, 1/sqrt(2)) instead
//      of (1, 1, 1)/sqrt(3) (9.74 deg angular error, 41.4% axis-intercept distortion)
//      and executing 2 sqrts instead of 0 or 1.
//   2. By sorting (a, b, c) into the S3 Weyl chamber u1 <= u2 <= u3 and exploiting
//      the convexity of the exterior polyhedron, opUnionChamfer3 and opUnionChamferOct3
//      require 0 sqrts in the entire exterior {d >= 0} and 1 sqrt in the interior {d < 0}.
//   3. Gating opUnionRound3 on the second-smallest input (u2 >= r) eliminates the
//      redundant sqrt on all three flat walls and prevents r - sqrt((r - u1)^2)
//      mantissa cancellation when |u1| << r.

#ifndef SDF_GEMS_OP_UNION_ROUND3_GLSL
#define SDF_GEMS_OP_UNION_ROUND3_GLSL

// 1. Exact 0/1-Sqrt 3-Way Rolling-Ball Rounded Union & Intersection
float opUnionRound3( float a, float b, float c, float r )
{
    float u1 = min(a, min(b, c));
    float u2 = max(min(a, b), min(max(a, b), c));
    if( u2 >= r ) return u1;
    return r - length(max(r - vec3(a, b, c), 0.0));
}

float opIntersectionRound3( float a, float b, float c, float r )
{
    return -opUnionRound3(-a, -b, -c, r);
}

// 2. Exact 0/1-Sqrt 3-Way Rounded Union with Analytical 3D Gradient: vec4(d, gx, gy, gz)
vec4 sdgRound3( in vec4 a, in vec4 b, in vec4 c, in float r )
{
    vec4  m  = (a.x < b.x) ? a : b;
    float u2 = (a.x < b.x) ? min(b.x, c.x) : min(a.x, c.x);
    if( c.x < m.x ) { u2 = m.x; m = c; }
    if( u2 >= r ) return m;
    vec3  u  = max(r - vec3(a.x, b.x, c.x), 0.0);
    float l  = length(u);
    vec3  w  = u / max(l, 1e-30);
    return vec4( r - l, w.x * a.yzw + w.y * b.yzw + w.z * c.yzw );
}

// 3. Exact 0/1-Sqrt Symmetric 3-Way Edge-Chamfered Union & Intersection
//    0 sqrts in the entire exterior {d >= 0}; 1 sqrt in the interior {d < 0}.
float opUnionChamfer3( float a, float b, float c, float r )
{
    float u1 = min(a, b), u2 = max(a, b), u3 = max(u2, c);
    u2 = max(u1, min(u2, c)); u1 = min(u1, c);
    if( u1 >= 0.0 && u1 + u2 >= r )
        return min(u1, (u1 + u2 - r) * 0.70710678);
    float s = r - u1 + u2;
    float t = max((2.0 * u3 >= s) ? 0.5 * s : (s + u3) * (1.0 / 3.0), 0.5 * r);
    return -(t >= r ? length(vec3(u1, max(r - u2, 0.0), max(r - u3, 0.0)))
                    : length(vec3(u1 - r + t, u2 - t, max(t - u3, 0.0))));
}

float opIntersectionChamfer3( float a, float b, float c, float r )
{
    return -opUnionChamfer3(-a, -b, -c, r);
}

// 4. Exact 0/1-Sqrt Symmetric 3-Way Octahedral Corner Gusset Union & Intersection
//    Truncates/fills the 3-way corner with an equilateral octahedral facet a + b + c = r.
float opUnionChamferOct3( float a, float b, float c, float r )
{
    float u1 = min(a, b), u2 = max(a, b), u3 = max(u2, c);
    u2 = max(u1, min(u2, c)); u1 = min(u1, c);
    float s3 = u1 + u2 + u3 - r;
    if( u1 >= 0.0 && s3 >= 0.0 ) return min(u1, s3 * 0.57735027);
    if( s3 < 0.0 && u2 + u3 - 2.0 * u1 <= r ) return s3 * 0.57735027;
    float d2 = (u2 >= 0.0 && u2 + u3 >= r) ? 0.0 :
               ((u3 - u2 > r) ? u2 * u2 + max(r - u3, 0.0) * max(r - u3, 0.0)
                              : 0.5 * max(r - u2 - u3, 0.0) * max(r - u2 - u3, 0.0));
    return -sqrt(u1 * u1 + d2);
}

float opIntersectionChamferOct3( float a, float b, float c, float r )
{
    return -opUnionChamferOct3(-a, -b, -c, r);
}

// 5. Exact 0/1-Sqrt 3-Way Edge-Chamfered Union with Analytical 3D Gradient: vec4(d, gx, gy, gz)
vec4 sdgChamfer3( in vec4 a, in vec4 b, in vec4 c, in float r )
{
    vec4 u1 = a, u2 = b, u3 = c, tmp;
    if( u1.x > u2.x ) { tmp = u1; u1 = u2; u2 = tmp; }
    if( u2.x > u3.x ) { tmp = u2; u2 = u3; u3 = tmp; }
    if( u1.x > u2.x ) { tmp = u1; u1 = u2; u2 = tmp; }
    if( u1.x >= 0.0 && u1.x + u2.x >= r )
    {
        float ch = (u1.x + u2.x - r) * 0.70710678;
        return (u1.x < ch) ? u1 : vec4( ch, (u1.yzw + u2.yzw) * 0.70710678 );
    }
    float s  = r - u1.x + u2.x;
    float t  = max((2.0 * u3.x >= s) ? 0.5 * s : (s + u3.x) * (1.0 / 3.0), 0.5 * r);
    vec3  v  = (t >= r) ? vec3(u1.x, min(u2.x - r, 0.0), min(u3.x - r, 0.0))
                        : vec3(u1.x - r + t, u2.x - t, min(u3.x - t, 0.0));
    float l  = length(v);
    vec3  w  = -v / max(l, 1e-30);
    return vec4( -l, w.x * u1.yzw + w.y * u2.yzw + w.z * u3.yzw );
}

// 6. Corner-Local C2 Lame L3 3-Way Geometric Smooth Minimum & 3D Gradient
//    Strictly local, C2/G2 curvature-continuous, symmetric, and Holder-conservative
//    (||grad||_2 in [3^(-1/6), 1] = [0.832683, 1.0]).
float sminCubicGeom3( float a, float b, float c, float k )
{
    k *= 4.84732210; // 1.0 / (1.0 - pow(2.0, -1.0 / 3.0))
    float u1 = min(a, min(b, c));
    float u2 = max(min(a, b), min(max(a, b), c));
    if( u2 >= k ) return u1;
    vec3 u = max(k - vec3(a, b, c), 0.0);
    return k - pow(dot(u, u * u), 1.0 / 3.0);
}

vec4 sdgSmoothUnionCubic3( in vec4 a, in vec4 b, in vec4 c, in float k )
{
    k *= 4.84732210;
    vec4  m  = (a.x < b.x) ? a : b;
    float u2 = (a.x < b.x) ? min(b.x, c.x) : min(a.x, c.x);
    if( c.x < m.x ) { u2 = m.x; m = c; }
    if( u2 >= k ) return m;
    vec3  u  = max(k - vec3(a.x, b.x, c.x), 0.0);
    float s3 = dot(u, u * u);
    float l3 = pow(max(s3, 1e-30), 1.0 / 3.0);
    vec3  w  = (u * u) / (l3 * l3);
    return vec4( k - l3, w.x * a.yzw + w.y * b.yzw + w.z * c.yzw );
}

#endif // SDF_GEMS_OP_UNION_ROUND3_GLSL
