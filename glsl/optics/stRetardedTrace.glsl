// ============================================================================
// stRetardedTrace — retarded-time sphere tracing: SDFs seen at a finite speed of light (c = 1)
//
//   A camera pixel at event (ro, t0) looks along the PAST null ray
//       x(l) = ro + l rd,   t(l) = t0 - l,   l >= 0,  |rd| = 1.
//   A rigidly moving body with rest shape d0 has the spacetime field d(x, t) = d0(q(x, t)), q = body-frame point.
//   Along the null ray q(l) is a curve with speed |dq/dl| = |rd + u(x, t)| (u = velocity of the body point under x),
//   so any step s with  integral_0^s |dq/dl| <= d0(q)  cannot cross the surface.
//
//   stRetardedStep        first order:  s = d / (1 + speedMax)
//   stRetardedStep2       second order: s = 2 d / (S + sqrt(S^2 + 2 A d)),  S = current |dq/dl|, A >= |d^2 q / dl^2|
//   stGalileanAim         uniform Galilean translation: the body-frame ray is the straight line along rd + v
//   stLorentzBoost / stLorentzToRest / stDoppler / stNormalToCamera / stSnapshotToRest
//                         uniform special-relativistic motion: trace the static SDF in the rest frame
//   stRotatingLocalStep / stRotatingStep
//                         rigid rotation: local step from the axis distance / from the exact body-frame ray speed
//   stTraceMovingBox*, stTraceRotatingBox, stTraceOscillatingBox
//                         concrete marchers for a rounded box (shape = vec4(half extents, rounding radius);
//                         vec4(0, 0, 0, R) is a sphere). They return vec2(l, steps) or vec4(l, steps, D, l_rest);
//                         l < 0 means no hit.
//
// @author Claude Fable 5.2 (sdf-gems, 2026)
// @exact  Uniform motion (Galilean re-aim, Lorentz rest-frame trace): same hit as a static sphere trace of the exact
//         rest SDF along the transformed ray (measured against closed forms, see docs/optics/stRetardedTrace.md).
//         Rotation / oscillation / generic steps: conservative bounds (0 overstep measured), not exact.
//         Physics: only the Lorentz path is physically exact; Galilean translation, rigid rotation and rigid
//         oscillation are kinematic light-delay models (and unphysical where a surface speed reaches 1).
// ============================================================================

#ifndef SDF_GEMS_ST_RETARDED_TRACE_GLSL
#define SDF_GEMS_ST_RETARDED_TRACE_GLSL

#ifndef ST_MAX_STEPS
#define ST_MAX_STEPS 512
#endif
#ifndef ST_HIT_EPS
#define ST_HIT_EPS 2e-5
#endif

// Exact SDF of a rounded box, shape = (half extents, rounding radius). shape = (0, 0, 0, R) is a sphere.
float stSdRoundBox( in vec3 q, in vec4 shape )
{
    vec3 a = abs(q) - shape.xyz;
    return length(max(a, 0.0)) + min(max(a.x, max(a.y, a.z)), 0.0) - shape.w;
}

// Rodrigues rotation of p about the unit axis n by angle ang.
vec3 stRotateAxis( in vec3 p, in vec3 n, in float ang )
{
    float c = cos(ang), s = sin(ang);
    return p * c + cross(n, p) * s + n * (dot(n, p) * (1.0 - c));
}

// ---------------------------------------------------------------------------------------------------------------
// Generic conservative steps along the past null ray
// ---------------------------------------------------------------------------------------------------------------

// First order: every body point under the ray moves slower than speedMax, so |dq/dl| <= 1 + speedMax.
float stRetardedStep( in float d, in float speedMax )
{
    return max(d, 0.0) / (1.0 + speedMax);
}

// Second order: `speed` = |dq/dl| at the current sample (exact, or an upper bound), `accel` >= |d^2 q / dl^2| ahead.
// Largest s with speed * s + accel * s^2 / 2 <= d, written without cancellation. accel = 0 gives d / speed.
float stRetardedStep2( in float d, in float speed, in float accel )
{
    float h = max(d, 0.0);
    return 2.0 * h / max(speed + sqrt(speed * speed + 2.0 * accel * h), 1e-30);
}

// ---------------------------------------------------------------------------------------------------------------
// (a) Uniform Galilean translation, body centre c0 + v t: d(x, t) = d0(x - c0 - v t)
// ---------------------------------------------------------------------------------------------------------------

// Body-frame ray. Returns k = |rd + v|: body-frame arc length = k * l (classical aberration).
float stGalileanAim( in vec3 ro, in float t0, in vec3 rd, in vec3 c0, in vec3 v, out vec3 roB, out vec3 rdB )
{
    vec3  u = rd + v;
    float k = length(u);
    roB = ro - c0 - v * t0;
    rdB = (k > 0.0) ? u / k : rd;
    return k;
}

// Exact: ordinary unit-step sphere tracing along the re-aimed ray.
vec2 stTraceMovingBoxGalilean( in vec3 ro, in float t0, in vec3 rd, in vec3 c0, in vec3 v, in vec4 shape, in float lMax )
{
    vec3 roB, rdB;
    float k = stGalileanAim(ro, t0, rd, c0, v, roB, rdB);
    float sMax = lMax * k, s = 0.0;
    for( int i = 0; i < ST_MAX_STEPS; i++ )
    {
        float h = stSdRoundBox(roB + s * rdB, shape);
        if( h < ST_HIT_EPS ) return vec2(s / max(k, 1e-30), float(i));
        s += h;
        if( s > sMax ) return vec2(-1.0, float(i + 1));
    }
    return vec2(-1.0, float(ST_MAX_STEPS));
}

// Conservative: camera-frame march of d(x, t) with the first-order step d / (1 + |v|).
vec2 stTraceMovingBoxRetarded( in vec3 ro, in float t0, in vec3 rd, in vec3 c0, in vec3 v, in vec4 shape, in float lMax )
{
    float sp = length(v), l = 0.0;
    for( int i = 0; i < ST_MAX_STEPS; i++ )
    {
        float h = stSdRoundBox(ro + l * rd - c0 - v * (t0 - l), shape);
        if( h < ST_HIT_EPS ) return vec2(l, float(i));
        l += stRetardedStep(h, sp);
        if( l > lMax ) return vec2(-1.0, float(i + 1));
    }
    return vec2(-1.0, float(ST_MAX_STEPS));
}

// ---------------------------------------------------------------------------------------------------------------
// (b) Uniform special-relativistic motion, velocity v (|v| < 1). Rest frame S': x' = x + g^2/(g+1) (x.v) v - g v t,
//     t' = g (t - v.x). The body is static in S' (its origin passes the camera-frame origin at t = 0).
// ---------------------------------------------------------------------------------------------------------------

float stGamma( in vec3 v )
{
    float b = length(v);
    return inversesqrt(max((1.0 - b) * (1.0 + b), 1e-12));
}

// Lorentz boost of the event e = (x, t) into the frame moving with velocity v. Inverse: stLorentzBoost(e', -v).
vec4 stLorentzBoost( in vec4 e, in vec3 v )
{
    float g  = stGamma(v);
    float xv = dot(e.xyz, v);
    return vec4(e.xyz + (g * g / (g + 1.0) * xv - g * e.w) * v, g * (e.w - xv));
}

// Doppler factor D = (observed frequency) / (rest-frame frequency) for a camera looking along rd at a body with
// velocity v: D = 1 / (gamma (1 + v.rd)). Radiance scales as D^3 (D^4 bolometric): relativistic beaming.
float stDoppler( in vec3 rd, in vec3 v )
{
    return 1.0 / (stGamma(v) * (1.0 + dot(v, rd)));
}

// Past null ray (ro, t0, rd) -> rest frame (roR, tR, rdR), rdR unit (relativistic aberration).
// Returns D; rest-frame ray parameter l' = l / D.
float stLorentzToRest( in vec3 ro, in float t0, in vec3 rd, in vec3 v, out vec3 roR, out float tR, out vec3 rdR )
{
    float g  = stGamma(v);
    vec4  e  = stLorentzBoost(vec4(ro, t0), v);
    float D  = 1.0 / (g * (1.0 + dot(v, rd)));
    roR = e.xyz;
    tR  = e.w;
    rdR = normalize(rd + (g * g / (g + 1.0) * dot(rd, v) + g) * v);
    return D;
}

// Rest-frame surface normal -> normal of the (Lorentz-contracted) constant-t surface in the camera frame.
vec3 stNormalToCamera( in vec3 nR, in vec3 v )
{
    float g = stGamma(v);
    return normalize(nR + (g * g / (g + 1.0) * dot(nR, v)) * v);
}

// "Snapshot" (infinite light speed) ray: the contracted body frozen at camera time t0, mapped to rest coordinates.
// Returns k = |d x'/dl|: rest-coordinate arc length = k * l. Not what a camera sees; provided for comparison.
float stSnapshotToRest( in vec3 ro, in float t0, in vec3 rd, in vec3 v, out vec3 roR, out vec3 rdR )
{
    float g = stGamma(v);
    float a = g * g / (g + 1.0);
    vec3  u = rd + a * dot(rd, v) * v;
    float k = length(u);
    roR = ro + (a * dot(ro, v) - g * t0) * v;
    rdR = u / k;
    return k;
}

// Exact rest-frame trace of a rounded box centred at cR in the rest frame (no rotation).
// Returns vec4(l, steps, D, l'): camera-frame hit ro + l rd at time t0 - l; rest-frame hit roR + l' rdR.
vec4 stTraceMovingBoxLorentz( in vec3 ro, in float t0, in vec3 rd, in vec3 v, in vec3 cR, in vec4 shape, in float lMax )
{
    vec3 roR, rdR; float tR;
    float D = stLorentzToRest(ro, t0, rd, v, roR, tR, rdR);
    float sMax = lMax / D, s = 0.0;
    roR -= cR;
    for( int i = 0; i < ST_MAX_STEPS; i++ )
    {
        float h = stSdRoundBox(roR + s * rdR, shape);
        if( h < ST_HIT_EPS ) return vec4(s * D, float(i), D, s);
        s += h;
        if( s > sMax ) return vec4(-1.0, float(i + 1), D, -1.0);
    }
    return vec4(-1.0, float(ST_MAX_STEPS), D, -1.0);
}

// Snapshot (infinite c) trace of the same body: the contracted shape at camera time t0. Returns vec2(l, steps).
vec2 stTraceMovingBoxSnapshot( in vec3 ro, in float t0, in vec3 rd, in vec3 v, in vec3 cR, in vec4 shape, in float lMax )
{
    vec3 roR, rdR;
    float k = stSnapshotToRest(ro, t0, rd, v, roR, rdR);
    float sMax = lMax * k, s = 0.0;
    roR -= cR;
    for( int i = 0; i < ST_MAX_STEPS; i++ )
    {
        float h = stSdRoundBox(roR + s * rdR, shape);
        if( h < ST_HIT_EPS ) return vec2(s / k, float(i));
        s += h;
        if( s > sMax ) return vec2(-1.0, float(i + 1));
    }
    return vec2(-1.0, float(ST_MAX_STEPS));
}

// ---------------------------------------------------------------------------------------------------------------
// (c) Rigid rotation about the unit axis n through the origin, angular velocity W: d(x, t) = d0(R(n, -W t) x - cB)
//     dq/dl = R (rd + W n x x) and d/dl (rd + W n x x) = W n x rd is CONSTANT along the ray.
// ---------------------------------------------------------------------------------------------------------------

// Local step from the axis distance rho of the current sample. sinAxis = |n x rd| (<= 1; pass 1.0 if unknown):
// the axis distance grows by at most sinAxis per unit l, so |dq/dl| <= 1 + |W| (rho + sinAxis s).
float stRotatingLocalStep( in float d, in float rho, in float W, in float sinAxis )
{
    return stRetardedStep2(d, 1.0 + abs(W) * rho, abs(W) * sinAxis);
}

// Tighter local step from the exact body-frame ray speed |rd + W n x x| (which can be < 1 when the ray runs
// against the rotation) and its exact rate |W n x rd|.
float stRotatingStep( in float d, in vec3 x, in vec3 rd, in vec3 n, in float W )
{
    return stRetardedStep2(d, length(rd + W * cross(n, x)), abs(W) * length(cross(n, rd)));
}

// mode 0: global bound 1 + |W| rB; mode 1: stRotatingLocalStep; mode 2: stRotatingStep.
// cB = body-frame offset of the box centre from the axis point. All modes march only inside the (rotation
// invariant) bounding sphere |x| <= rB = |cB| + |half extents| + rounding. Returns vec2(l, steps).
vec2 stTraceRotatingBox( in vec3 ro, in float t0, in vec3 rd, in vec3 n, in float W, in vec3 cB, in vec4 shape, in int mode )
{
    float rB = length(cB) + length(shape.xyz) + shape.w;
    float b  = dot(ro, rd);
    float h2 = (rB * rB) - (dot(ro, ro) - b * b);
    if( h2 <= 0.0 ) return vec2(-1.0, 0.0);
    float hh = sqrt(h2);
    float l  = max(-b - hh, 0.0), l1 = -b + hh;
    if( l1 <= 0.0 ) return vec2(-1.0, 0.0);
    float sinAxis = length(cross(n, rd));
    float aW = abs(W);
    for( int i = 0; i < ST_MAX_STEPS; i++ )
    {
        vec3  x = ro + l * rd;
        float h = stSdRoundBox(stRotateAxis(x, n, -W * (t0 - l)) - cB, shape);
        if( h < ST_HIT_EPS ) return vec2(l, float(i));
        if( mode == 0 )      l += stRetardedStep(h, aW * rB);
        else if( mode == 1 ) l += stRotatingLocalStep(h, length(x - n * dot(n, x)), W, sinAxis);
        else                 l += stRetardedStep2(h, length(rd + W * cross(n, x)), aW * sinAxis);
        if( l > l1 ) return vec2(-1.0, float(i + 1));
    }
    return vec2(-1.0, float(ST_MAX_STEPS));
}

// ---------------------------------------------------------------------------------------------------------------
// (d) Rigid oscillation, body centre A sin(w t): d(x, t) = d0(x - A sin(w t))
//     dq/dl = rd + A w cos(w t),  |d^2 q / dl^2| <= |A| w^2.
// ---------------------------------------------------------------------------------------------------------------

// mode 0: global bound 1 + |A| w; mode 1: second-order step from the exact current speed. Returns vec2(l, steps).
vec2 stTraceOscillatingBox( in vec3 ro, in float t0, in vec3 rd, in vec3 A, in float w, in vec4 shape, in float lMax, in int mode )
{
    float aw = length(A) * abs(w), l = 0.0;
    for( int i = 0; i < ST_MAX_STEPS; i++ )
    {
        float ph = w * (t0 - l);
        float h  = stSdRoundBox(ro + l * rd - A * sin(ph), shape);
        if( h < ST_HIT_EPS ) return vec2(l, float(i));
        if( mode == 0 ) l += stRetardedStep(h, aw);
        else            l += stRetardedStep2(h, length(rd + A * (w * cos(ph))), aw * abs(w));
        if( l > lMax ) return vec2(-1.0, float(i + 1));
    }
    return vec2(-1.0, float(ST_MAX_STEPS));
}

#endif // SDF_GEMS_ST_RETARDED_TRACE_GLSL
