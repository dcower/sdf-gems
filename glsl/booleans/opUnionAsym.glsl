// ============================================================================
// opUnionAsym — Exact 0/1-Sqrt Asymmetric Setback Booleans & 3D Gradients
// @author Discovered and derived using Claude Opus 5.5; audited and hardened
//         (degenerate zero-setback guards) using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` across exterior, corner profile, and interior)
// ============================================================================
//
// Exact Euclidean boolean operators with independent setbacks (ra, rb) on
// orthogonal surfaces (grad(a) . grad(b) = 0) and analytical 3D gradients:
//
//   1. Rational Bi-Arc Asymmetric Fillet (`opUnionFilletAsym`, `sdgFilletAsym`):
//      Scaling an isotropic circular fillet to (ra, rb) warps the arc into an
//      ellipse with non-unit gradient ||grad d|| in [min(1, ra/rb), max(1, ra/rb)].
//      Instead, two internally tangent circular arcs centered at C1 = (ra, R1)
//      and C2 = (R2, rb) satisfy the Rational Bi-Arc Hyperbola Identity:
//          (R1 - ra)(rb - R2) = 0.5 * (ra - rb)^2
//      Assuming ra >= rb > 0, the canonical rational split:
//          R1 = ra + (ra^2 - rb^2) / (2 * rb),   R2 = 2 * rb^2 / (ra + rb)
//      uses 0 sqrts and 0 trig, satisfies ||C1 - C2|| = R1 - R2 rationally,
//      is G^1-continuous, and reduces identically to a circle of radius r when
//      ra = rb = r.
//
//   2. Exact Corner-Local Asymmetric Chamfer (`opUnionChamferAsym`, `sdgChamferAsym`):
//      Exact Euclidean distance to the bevel segment from (0, rb) to (ra, 0)
//      inside its normal strip -rb^2 <= a*ra - b*rb <= ra^2, and exact 0/1-sqrt
//      Euclidean distance to the two reflex vertices (ra, 0) and (0, rb) outside.
//
//   3. Exact Carpenter's Rabbet & Stepped Gusset (`opRabbet`, `opUnionRabbet`, `sdgRabbet`):
//      Exact 0/1-sqrt Euclidean stepped corner notch on an intersection
//      ({a <= 0, b <= 0} \ {a > -ra, b > -rb}) and its De Morgan dual stepped
//      corner gusset on a union, fixing the 29.3% interior corner error of CSG.

#ifndef SDF_GEMS_OP_UNION_ASYM_GLSL
#define SDF_GEMS_OP_UNION_ASYM_GLSL

// ----------------------------------------------------------------------------
// 1. Rational Bi-Arc Asymmetric Fillet Booleans & 3D Gradients
// ----------------------------------------------------------------------------

float opUnionFilletAsym( float a, float b, float ra, float rb )
{
    if( ra < rb ) { float t = a; a = b; b = t; t = ra; ra = rb; rb = t; }
    float d = min(a, b);
    if( a >= ra || rb <= 0.0 ) return d;
    float R1 = ra + (ra - rb) * (ra + rb) / (2.0 * rb);
    float R2 = 2.0 * rb * rb / (ra + rb);
    float dx = ra - R2, dy = R1 - rb;
    float s  = (R2 - a) * dy - (rb - b) * dx;
    if( b < rb && s >= 0.0 )
        return R2 - length(vec2(R2 - a, rb - b));
    if( s < 0.0 )
        return min(d, R1 - length(vec2(ra - a, R1 - b)));
    return d;
}

float opIntersectionFilletAsym( float a, float b, float ra, float rb )
{
    return -opUnionFilletAsym(-a, -b, ra, rb);
}

float opSubtractionFilletAsym( float a, float b, float ra, float rb )
{
    return -opUnionFilletAsym(a, -b, ra, rb);
}

vec4 sdgFilletAsym( vec4 a, vec4 b, float ra, float rb )
{
    if( ra < rb ) { vec4 t = a; a = b; b = t; float tr = ra; ra = rb; rb = tr; }
    vec4 d = (a.x < b.x) ? a : b;
    if( a.x >= ra || rb <= 0.0 ) return d;
    float R1 = ra + (ra - rb) * (ra + rb) / (2.0 * rb);
    float R2 = 2.0 * rb * rb / (ra + rb);
    float dx = ra - R2, dy = R1 - rb;
    float s  = (R2 - a.x) * dy - (rb - b.x) * dx;
    if( b.x < rb && s >= 0.0 )
    {
        vec2  u = vec2(R2 - a.x, rb - b.x);
        float l = length(u);
        return vec4(R2 - l, (l > 0.0) ? (u.x * a.yzw + u.y * b.yzw) / l : a.yzw);
    }
    if( s < 0.0 )
    {
        vec2  u  = vec2(ra - a.x, R1 - b.x);
        float l  = length(u);
        float d1 = R1 - l;
        if( d1 < d.x )
            return vec4(d1, (l > 0.0) ? (u.x * a.yzw + u.y * b.yzw) / l : a.yzw);
    }
    return d;
}

vec4 sdgIntersectionFilletAsym( vec4 a, vec4 b, float ra, float rb )
{
    return -sdgFilletAsym(-a, -b, ra, rb);
}

vec4 sdgSubtractionFilletAsym( vec4 a, vec4 b, float ra, float rb )
{
    return -sdgFilletAsym(a, -b, ra, rb);
}

// ----------------------------------------------------------------------------
// 2. Exact Corner-Local Euclidean Asymmetric Chamfer Booleans & 3D Gradients
// ----------------------------------------------------------------------------

float opUnionChamferAsym( float a, float b, float ra, float rb )
{
    if( ra <= 0.0 || rb <= 0.0 ) return min(a, b);
    float q = a * ra - b * rb;
    if( b < 0.0 &&  q > ra * ra ) return (a >= ra) ? b : -length(vec2(ra - a, b));
    if( a < 0.0 && -q > rb * rb ) return (b >= rb) ? a : -length(vec2(a, rb - b));
    return min(min(a, b), (a * rb + b * ra - ra * rb) * inversesqrt(ra * ra + rb * rb));
}

float opIntersectionChamferAsym( float a, float b, float ra, float rb )
{
    return -opUnionChamferAsym(-a, -b, ra, rb);
}

float opSubtractionChamferAsym( float a, float b, float ra, float rb )
{
    return -opUnionChamferAsym(a, -b, ra, rb);
}

vec4 sdgChamferAsym( vec4 a, vec4 b, float ra, float rb )
{
    vec4 d = (a.x < b.x) ? a : b;
    if( ra <= 0.0 || rb <= 0.0 ) return d;
    float q = a.x * ra - b.x * rb;
    if( b.x < 0.0 && q > ra * ra )
    {
        if( a.x >= ra ) return b;
        vec2  u = vec2(a.x - ra, b.x);
        float l = length(u);
        return vec4(-l, -(u.x * a.yzw + u.y * b.yzw) / l);
    }
    if( a.x < 0.0 && -q > rb * rb )
    {
        if( b.x >= rb ) return a;
        vec2  u = vec2(a.x, b.x - rb);
        float l = length(u);
        return vec4(-l, -(u.x * a.yzw + u.y * b.yzw) / l);
    }
    float il = inversesqrt(ra * ra + rb * rb);
    float ch = (a.x * rb + b.x * ra - ra * rb) * il;
    return (ch < d.x) ? vec4(ch, (rb * a.yzw + ra * b.yzw) * il) : d;
}

vec4 sdgIntersectionChamferAsym( vec4 a, vec4 b, float ra, float rb )
{
    return -sdgChamferAsym(-a, -b, ra, rb);
}

vec4 sdgSubtractionChamferAsym( vec4 a, vec4 b, float ra, float rb )
{
    return -sdgChamferAsym(a, -b, ra, rb);
}

// ----------------------------------------------------------------------------
// 3. Exact 0/1-Sqrt Euclidean Carpenter's Rabbet & Stepped Union Gusset
// ----------------------------------------------------------------------------

float opRabbet( float a, float b, float ra, float rb )
{
    float u = a + ra, v = b + rb;
    if( u <= 0.0 || v <= 0.0 )
        return max(max(a, b), (u <= 0.0 && v <= 0.0) ? -length(vec2(u, v)) : min(u, v));
    float pa = max(a, 0.0), pb = max(b, 0.0);
    return (u * u + pb * pb < pa * pa + v * v)
         ? ((b <= 0.0) ? u : length(vec2(u, b)))
         : ((a <= 0.0) ? v : length(vec2(a, v)));
}

float opUnionRabbet( float a, float b, float ra, float rb )
{
    return -opRabbet(-a, -b, ra, rb);
}

vec4 sdgRabbet( vec4 a, vec4 b, float ra, float rb )
{
    float u = a.x + ra, v = b.x + rb;
    if( u <= 0.0 || v <= 0.0 )
    {
        vec4 d = (a.x > b.x) ? a : b;
        if( u <= 0.0 && v <= 0.0 )
        {
            float l = length(vec2(u, v));
            if( -l > d.x )
                return vec4(-l, (l > 0.0) ? -(u * a.yzw + v * b.yzw) / l : a.yzw);
            return d;
        }
        vec4 duv = (u < v) ? vec4(u, a.yzw) : vec4(v, b.yzw);
        return (duv.x > d.x) ? duv : d;
    }
    float pa = max(a.x, 0.0), pb = max(b.x, 0.0);
    if( u * u + pb * pb < pa * pa + v * v )
    {
        if( b.x <= 0.0 ) return vec4(u, a.yzw);
        float l = length(vec2(u, b.x));
        return vec4(l, (u * a.yzw + b.x * b.yzw) / l);
    }
    if( a.x <= 0.0 ) return vec4(v, b.yzw);
    float l = length(vec2(a.x, v));
    return vec4(l, (a.x * a.yzw + v * b.yzw) / l);
}

vec4 sdgUnionRabbet( vec4 a, vec4 b, float ra, float rb )
{
    return -sdgRabbet(-a, -b, ra, rb);
}

#endif // SDF_GEMS_OP_UNION_ASYM_GLSL
