// ============================================================================
// opRepPolar — 0-Sqrt Direct-Rotation Infinite/Limited Polar Repetition & 0-Trig Dihedral/Cyclic Folding
// @author Mercury Demogroup — upgraded to 0-sqrt Direct Rotation, Voronoi-Centered Limited Wrapping, & 0-Trig D_n/C_n Folding using Gemini 4 Argon (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (0 sqrts; eliminates theta = +-pi branch cut, asymmetric [lima, limb] Voronoi tear, and boundary neighbor drop)
// ============================================================================

#ifndef SDF_GEMS_OP_REP_POLAR_GLSL
#define SDF_GEMS_OP_REP_POLAR_GLSL

// 1. 0-Sqrt 1-Tap Infinite Polar Repetition for radially symmetric primitives
vec2 opRepPolar( vec2 p, float n, out float id )
{
    float sp  = 6.28318530718 / n;
    float a   = atan(p.y, p.x);
    id = round(a / sp);
    float ang = id * sp;
    float c   = cos(ang), s = sin(ang);
    return vec2(c * p.x + s * p.y, c * p.y - s * p.x);
}

vec2 opRepPolar( vec2 p, float n )
{
    float id;
    return opRepPolar(p, n, id);
}

// 2. 0-Sqrt 1-Tap Limited Polar Repetition on [lima, limb] (Voronoi-centered branch cut)
vec2 opRepPolarLim( vec2 p, float n, float lima, float limb, out float id )
{
    float sp     = 6.28318530718 / n;
    float a      = atan(p.y, p.x) / sp;
    float mid    = 0.5 * (lima + limb);
    float sector = mid + mod(a - mid + 0.5 * n, n) - 0.5 * n;
    id = clamp(round(sector), lima, limb);
    float ang    = id * sp;
    float c      = cos(ang), s = sin(ang);
    return vec2(c * p.x + s * p.y, c * p.y - s * p.x);
}

vec2 opRepPolarLim( vec2 p, float n, float lima, float limb )
{
    float id;
    return opRepPolarLim(p, n, lima, limb, id);
}

// 3. 2-Tap Limited Polar Repetition on [lima, limb] for asymmetric primitives (outputs angles id0, id1 in radians)
void opLimitedPolarRep2Tap( in vec2 p, in float n, in float lima, in float limb, out float id0, out float id1 )
{
    float sp     = 6.28318530718 / n;
    float a      = atan(p.y, p.x) / sp;
    float mid    = 0.5 * (lima + limb);
    float span   = limb - lima;
    float sector = mid + mod(a - mid + 0.5 * n, n) - 0.5 * n;
    float margin = 0.25 * max(n - span - 1.0, 0.0);
    float i0, i1;
    if( abs(sector - mid) <= 0.5 * span + margin )
    {
        i0 = clamp(floor(sector), lima, max(lima, limb - 1.0));
        i1 = min(i0 + 1.0, limb);
    }
    else
    {
        i0 = lima;
        i1 = limb;
    }
    id0 = i0 * sp;
    id1 = i1 * sp;
}

// Symmetric [-max_idx, max_idx] overload (outputs angles id0, id1 in radians)
void opLimitedPolarRep2Tap( in vec2 p, in float n, in float max_idx, out float id0, out float id1 )
{
    opLimitedPolarRep2Tap(p, n, -max_idx, max_idx, id0, id1);
}

// 4. 0-Trig, 0-Sqrt Dihedral Sector Folds (D_n into [0, pi/n]) & Proper Cyclic Rotations (C_n into [-pi/n, pi/n])
vec2 opFoldPolar3( vec2 p )
{
    p.y = abs(p.y);
    const vec2 n = vec2(-0.86602540378, 0.5);
    p -= 2.0 * max(dot(p, n), 0.0) * n;
    p.y = abs(p.y);
    return p;
}

vec2 opRepPolar3( vec2 p )
{
    float s = (p.y < 0.0) ? -1.0 : 1.0;
    p.y = abs(p.y);
    const vec2 n = vec2(-0.86602540378, 0.5);
    float h = dot(p, n);
    if( h > 0.0 ) { p -= 2.0 * h * n; s = -s; }
    p.y *= s;
    return p;
}

vec2 opFoldPolar4( vec2 p )
{
    p = abs(p);
    return (p.y > p.x) ? p.yx : p;
}

vec2 opRepPolar4( vec2 p )
{
    if( abs(p.y) > abs(p.x) ) return (p.y > 0.0) ? vec2(p.y, -p.x) : vec2(-p.y, p.x);
    return (p.x < 0.0) ? -p : p;
}

vec2 opFoldPolar5( vec2 p )
{
    p.y = abs(p.y);
    const vec2 n1 = vec2(-0.95105651630, 0.30901699437);
    const vec2 n2 = vec2(-0.58778525229, 0.80901699437);
    p -= 2.0 * max(dot(p, n1), 0.0) * n1;
    p -= 2.0 * max(dot(p, n2), 0.0) * n2;
    p.y = abs(p.y);
    return p;
}

vec2 opRepPolar5( vec2 p )
{
    float s = (p.y < 0.0) ? -1.0 : 1.0;
    p.y = abs(p.y);
    const vec2 n1 = vec2(-0.95105651630, 0.30901699437);
    const vec2 n2 = vec2(-0.58778525229, 0.80901699437);
    float h1 = dot(p, n1);
    if( h1 > 0.0 ) { p -= 2.0 * h1 * n1; s = -s; }
    float h2 = dot(p, n2);
    if( h2 > 0.0 ) { p -= 2.0 * h2 * n2; s = -s; }
    p.y *= s;
    return p;
}

vec2 opFoldPolar6( vec2 p )
{
    p = abs(p);
    const vec2 n = vec2(-0.5, 0.86602540378);
    p -= 2.0 * max(dot(p, n), 0.0) * n;
    p.y = abs(p.y);
    return p;
}

vec2 opRepPolar6( vec2 p )
{
    float s = 1.0;
    if( p.y < 0.0 ) { p.y = -p.y; s = -s; }
    if( p.x < 0.0 ) { p.x = -p.x; s = -s; }
    const vec2 n = vec2(-0.5, 0.86602540378);
    float h = dot(p, n);
    if( h > 0.0 ) { p -= 2.0 * h * n; s = -s; }
    p.y *= s;
    return p;
}

vec2 opFoldPolar8( vec2 p )
{
    p = abs(p);
    if( p.y > p.x ) p.xy = p.yx;
    const vec2 n = vec2(-0.38268343237, 0.92387953251);
    p -= 2.0 * max(dot(p, n), 0.0) * n;
    return p;
}

vec2 opRepPolar8( vec2 p )
{
    float s = 1.0;
    if( p.y < 0.0 ) { p.y = -p.y; s = -s; }
    if( p.x < 0.0 ) { p.x = -p.x; s = -s; }
    if( p.y > p.x ) { p.xy = p.yx; s = -s; }
    const vec2 n = vec2(-0.38268343237, 0.92387953251);
    float h = dot(p, n);
    if( h > 0.0 ) { p -= 2.0 * h * n; s = -s; }
    p.y *= s;
    return p;
}

vec2 opFoldPolar12( vec2 p )
{
    p = abs(p);
    if( p.y > p.x ) p.xy = p.yx;
    const vec2 n = vec2(-0.25881904510, 0.96592582629);
    p -= 2.0 * max(dot(p, n), 0.0) * n;
    p.y = abs(p.y);
    return p;
}

vec2 opRepPolar12( vec2 p )
{
    float s = 1.0;
    if( p.y < 0.0 ) { p.y = -p.y; s = -s; }
    if( p.x < 0.0 ) { p.x = -p.x; s = -s; }
    if( p.y > p.x ) { p.xy = p.yx; s = -s; }
    const vec2 n = vec2(-0.25881904510, 0.96592582629);
    float h = dot(p, n);
    if( h > 0.0 ) { p -= 2.0 * h * n; s = -s; }
    p.y *= s;
    return p;
}

#endif // SDF_GEMS_OP_REP_POLAR_GLSL
