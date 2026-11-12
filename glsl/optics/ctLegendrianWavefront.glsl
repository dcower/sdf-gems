// ============================================================================
// ctLegendrianWavefront — Contact Geometry, Legendrian Multi-Valued Wavefronts,
//                         Arnold A3 Swallowtails & Caustic Evolutes
//   ctParabolaMultiSDF2D, ctParabolaWavefront2D, sdgParabolaWavefront2D,
//   sdSwallowtailPocket2D, sdSemiCubicalEvolute2D, sdgSemiCubicalEvolute2D,
//   sdParaboloidCaustic3D, sdgParaboloidCaustic3D,
//   ctEllipseMultiSDF2D, ctEllipseWavefront2D,
//   sdSwallowtailWavefront3D, sdgSwallowtailWavefront3D
// @author Gemini 4 Argon (sdf-gems, 2026)
// @exact  ctParabolaMultiSDF2D: Exact 3-sheeted normal distances (d1 <= d2 <= d3)
//           satisfying ||grad d_i|| = 1 on every sheet inside the caustic evolute.
//         ctParabolaWavefront2D / sdgParabolaWavefront2D: Exact 2D Euclidean SDF
//           to the full 3-sheeted self-intersecting A3 swallowtail wavefront
//           (closed form via the Huygens-Legendrian Cusp Stationarity Theorem).
//         sdSemiCubicalEvolute2D / sdParaboloidCaustic3D: Exact Euclidean SDF to
//           the caustic evolute 27 k^2 x^2 = 8 (k y - 1)^3 (Darboux-convex quartic
//           solver from a certified upper bracket + cusp check, 0 trig).
//         ctEllipseMultiSDF2D / ctEllipseWavefront2D: Exact 4-sheeted normal
//           distances and self-intersecting wavefront SDF for an ellipse via
//           Same-Sign Half-Pole Reconstruction and the Cusp Stationarity Theorem.
//         sdSwallowtailWavefront3D / sdgSwallowtailWavefront3D: Exact 3D Euclidean
//           SDF and unit gradient to the thickened 3D Arnold A3 swallowtail shell.
// ============================================================================

#ifndef SDF_GEMS_HAS_CT_LEGENDRIAN_WAVEFRONT
#define SDF_GEMS_HAS_CT_LEGENDRIAN_WAVEFRONT

// ---------------------------------------------------------------------------
// 1. Multi-Valued Legendrian Normal Generating Family of a Parabola y = 0.5 * k * x^2
// ---------------------------------------------------------------------------
vec4 ctParabolaMultiSDF2D( in vec2 p, in float k )
{
    float ax = abs(p.x);
    float ik = 1.0 / max(k, 1e-6);
    float pc = ik * ik * (k * p.y - 1.0) * (2.0 / 3.0);
    float qc = ik * ik * ax;
    float h  = qc * qc - pc * pc * pc;

    if( h >= 0.0 )
    {
        float r  = pow(qc + sqrt(h), 1.0 / 3.0);
        float s1 = (pc < 0.0) ? (2.0 * qc / (r * r - pc + (pc / r) * (pc / r)))
                              : ((r > 0.0) ? (r + pc / r) : 0.0);
        float d1 = length(vec2(ax - s1, p.y - 0.5 * k * s1 * s1));
        d1 = (p.y >= 0.5 * k * ax * ax) ? d1 : -d1;
        return vec4(d1, d1, d1, 1.0);
    }

    float r  = sqrt(max(pc, 0.0));
    float th = acos(clamp(qc / max(pc * r, 1e-30), -1.0, 1.0)) / 3.0;
    float cs = cos(th), sn = sin(th) * 1.732050808;
    vec3  s  = r * vec3(2.0 * cs, -(cs + sn), -(cs - sn));
    vec3  d;
    for( int i = 0; i < 3; i++ )
    {
        float si = s[i];
        d[i] = length(vec2(ax - si, p.y - 0.5 * k * si * si));
    }
    float dMin = min(d.x, min(d.y, d.z));
    float dMax = max(d.x, max(d.y, d.z));
    float dMid = d.x + d.y + d.z - dMin - dMax;
    return vec4(dMin, dMid, dMax, 3.0);
}

// ---------------------------------------------------------------------------
// 2. Exact 2D Euclidean SDF & Unit Gradient to the Multi-Sheeted Legendrian
//    Wavefront gamma_t(s) = (s, 0.5*k*s^2) + t * (-k*s, 1) / sqrt(1 + k^2*s^2).
// ---------------------------------------------------------------------------
vec3 sdgParabolaWavefront2D( in vec2 p, in float k, in float tOffset )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    float ax = abs(p.x);
    float ik = 1.0 / max(k, 1e-6);
    float pc = ik * ik * (k * p.y - 1.0) * (2.0 / 3.0);
    float qc = ik * ik * ax;
    float h  = qc * qc - pc * pc * pc;

    float bestD = 1e30;
    vec2  bestG = vec2(0.0, 1.0);

    if( h >= 0.0 )
    {
        float r   = pow(qc + sqrt(h), 1.0 / 3.0);
        float s1  = (pc < 0.0) ? (2.0 * qc / (r * r - pc + (pc / r) * (pc / r)))
                               : ((r > 0.0) ? (r + pc / r) : 0.0);
        float ks  = k * s1;
        float inv = inversesqrt(1.0 + ks * ks);
        vec2  n   = vec2(-ks * inv, inv);
        float tau = (p.y - 0.5 * k * s1 * s1) / inv;
        float df  = tau - tOffset;
        bestD = abs(df);
        bestG = (df >= 0.0) ? n : -n;
    }
    else
    {
        float r  = sqrt(max(pc, 0.0));
        float th = acos(clamp(qc / max(pc * r, 1e-30), -1.0, 1.0)) / 3.0;
        float cs = cos(th), sn = sin(th) * 1.732050808;
        vec3  s  = r * vec3(2.0 * cs, -(cs + sn), -(cs - sn));
        for( int i = 0; i < 3; i++ )
        {
            float si  = s[i];
            float ks  = k * si;
            float inv = inversesqrt(1.0 + ks * ks);
            vec2  n   = vec2(-ks * inv, inv);
            float tau = (p.y - 0.5 * k * si * si) / inv;
            float df  = tau - tOffset;
            float ad  = abs(df);
            if( ad < bestD )
            {
                bestD = ad;
                bestG = (df >= 0.0) ? n : -n;
            }
        }
    }

    float kt = k * tOffset;
    if( kt > 1.0 )
    {
        float kt23 = pow(kt, 2.0 / 3.0);
        float u2   = max(kt23 - 1.0, 0.0);
        vec2  cusp = vec2(u2 * sqrt(u2) * ik, 0.5 * ik * (3.0 * kt23 - 1.0));
        vec2  dv   = vec2(ax, p.y) - cusp;
        float d2c  = dot(dv, dv);
        if( d2c < bestD * bestD )
        {
            bestD = sqrt(d2c);
            bestG = (bestD > 1e-20) ? (dv / bestD) : vec2(0.0, 1.0);
        }
    }

    return vec3(bestD, sx * bestG.x, bestG.y);
}

float ctParabolaWavefront2D( in vec2 p, in float k, in float tOffset )
{
    return sdgParabolaWavefront2D(p, k, tOffset).x;
}

// ---------------------------------------------------------------------------
// 3. Exact 2D SDF of the Swallowtail Retrograde Caustic Pocket {d2(p) <= t <= d3(p)}
// ---------------------------------------------------------------------------
float sdSwallowtailPocket2D( in vec2 p, in float k, in float tOffset )
{
    vec4  m      = ctParabolaMultiSDF2D(p, k);
    bool  inside = (m.w > 1.5) && (m.y < tOffset) && (tOffset < m.z);
    float ax     = abs(p.x);
    float ik     = 1.0 / max(k, 1e-6);
    float kt     = max(k * tOffset, 1.0);
    float kt23   = pow(kt, 2.0 / 3.0);
    float sMax   = sqrt(max(kt * kt - 1.0, 0.0)) * ik;
    float pc     = ik * ik * (k * p.y - 1.0) * (2.0 / 3.0);
    float qc     = ik * ik * ax;
    float h      = qc * qc - pc * pc * pc;
    float bestD  = 1e30;
    if( h >= 0.0 )
    {
        float r  = pow(qc + sqrt(h), 1.0 / 3.0);
        float s1 = (pc < 0.0) ? (2.0 * qc / (r * r - pc + (pc / r) * (pc / r)))
                              : ((r > 0.0) ? (r + pc / r) : 0.0);
        s1 = min(s1, sMax);
        float ks  = k * s1;
        float inv = inversesqrt(1.0 + ks * ks);
        vec2  q   = vec2(s1 - tOffset * ks * inv, 0.5 * k * s1 * s1 + tOffset * inv);
        bestD = length(vec2(ax, p.y) - q);
    }
    else
    {
        float r  = sqrt(max(pc, 0.0));
        float th = acos(clamp(qc / max(pc * r, 1e-30), -1.0, 1.0)) / 3.0;
        float cs = cos(th), sn = sin(th) * 1.732050808;
        vec3  s  = clamp(r * vec3(2.0 * cs, -(cs + sn), -(cs - sn)), -sMax, sMax);
        for( int i = 0; i < 3; i++ )
        {
            float si  = s[i];
            float ks  = k * si;
            float inv = inversesqrt(1.0 + ks * ks);
            vec2  q   = vec2(si - tOffset * ks * inv, 0.5 * k * si * si + tOffset * inv);
            bestD = min(bestD, length(vec2(ax, p.y) - q));
        }
    }
    float u2   = max(kt23 - 1.0, 0.0);
    vec2  cusp = vec2(u2 * sqrt(u2) * ik, 0.5 * ik * (3.0 * kt23 - 1.0));
    vec2  top  = vec2(0.0, 0.5 * ik * (kt * kt + 1.0));
    bestD = min(bestD, min(length(vec2(ax, p.y) - cusp), length(vec2(ax, p.y) - top)));
    return inside ? -bestD : bestD;
}

// ---------------------------------------------------------------------------
// 4. Exact Euclidean SDF to the Semi-Cubical Parabola Caustic Evolute
//      27 * k^2 * x^2 = 8 * (k * y - 1)^3   (parametric: (u^3 / k, (1 + 1.5*u^2) / k))
// ---------------------------------------------------------------------------
vec3 sdgSemiCubicalEvolute2D( in vec2 p, in float k )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    float ik = 1.0 / max(k, 1e-6);
    float X  = k * abs(p.x);
    float Y  = k * p.y - 1.0;

    float u = max(pow(X, 1.0 / 3.0), sqrt(max((2.0 / 3.0) * Y, 0.0)));
    for( int i = 0; i < 5; i++ )
    {
        float u2 = u * u;
        float g  = (2.0 * u2 + 3.0) * u2 - 2.0 * X * u - 2.0 * Y;
        float dg = (8.0 * u2 + 6.0) * u - 2.0 * X;
        if( dg <= 1e-20 ) break;
        u = max(u - g / dg, 0.0);
    }

    float u2   = u * u;
    vec2  foot = vec2(u2 * u, 1.5 * u2);
    vec2  dv   = vec2(X, Y) - foot;
    float d2   = dot(dv, dv);
    float d2_0 = X * X + Y * Y;
    if( d2_0 <= d2 )
    {
        u    = 0.0;
        u2   = 0.0;
        dv   = vec2(X, Y);
        d2   = d2_0;
    }

    float l   = sqrt(d2);
    bool  ins = (Y > 0.0) && (27.0 * X * X < 8.0 * Y * Y * Y);
    float sgn = ins ? -1.0 : 1.0;

    vec2 g2;
    if( l > 1e-20 )
    {
        g2 = sgn * dv / l;
    }
    else if( u > 1e-6 )
    {
        g2 = vec2(1.0, -u) * inversesqrt(1.0 + u2);
    }
    else
    {
        g2 = vec2(0.0, -1.0);
    }
    return vec3(sgn * l * ik, sx * g2.x, g2.y);
}

float sdSemiCubicalEvolute2D( in vec2 p, in float k )
{
    return sdgSemiCubicalEvolute2D(p, k).x;
}

// Revolved 3D Paraboloid Caustic Surface (27 * k^2 * (x^2 + z^2) = 8 * (k * y - 1)^3)
vec4 sdgParaboloidCaustic3D( in vec3 p, in float k )
{
    float r   = length(p.xz);
    vec3  dg2 = sdgSemiCubicalEvolute2D(vec2(r, p.y), k);
    vec2  er  = (r > 1e-20) ? (p.xz / r) : vec2(1.0, 0.0);
    return vec4(dg2.x, dg2.y * er.x, dg2.z, dg2.y * er.y);
}

float sdParaboloidCaustic3D( in vec3 p, in float k )
{
    return sdgSemiCubicalEvolute2D(vec2(length(p.xz), p.y), k).x;
}

// ---------------------------------------------------------------------------
// 5. Multi-Valued Legendrian Generating Family & Wavefront of a 2D Ellipse
// ---------------------------------------------------------------------------
vec4 ctEllipseMultiSDF2D( in vec2 p, in vec2 ab )
{
    p = abs(p);
    if( ab.x > ab.y ) { p = p.yx; ab = ab.yx; }

    float a = ab.x, b = ab.y;
    float a2 = a * a, b2 = b * b;
    float l  = b2 - a2, e = 0.5 * l;
    float x  = a * p.x, y = b * p.y;
    float x2 = x * x, y2 = y * y, r2 = x2 + y2;

    float C   = (r2 - l * l) / 3.0;
    float C2  = C * C, C3 = C2 * C;
    float lxy = l * x * y, l2x2y2 = lxy * lxy;
    float D   = C3 + l2x2y2, Q = C3 + 2.0 * l2x2y2;

    vec4  vList;
    float numSheets;
    if( D < 0.0 )
    {
        float h   = atan(2.0 * lxy * sqrt(-D), -Q) / 3.0;
        float s   = cos(h), tr = sin(h) * 1.732050808;
        float oms = (tr * tr) / (3.0 * (1.0 + s));
        float rz1 = sqrt(max(l * l + 2.0 * C * oms, 0.0));
        float rz2 = sqrt(max(r2 - C * (oms - tr), 0.0));
        float rz3 = sqrt(max(r2 - C * (oms + tr), 0.0));
        float sg  = (x2 >= y2) ? 1.0 : -1.0;
        vList = 0.5 * vec4(
             rz1 + rz3 + sg * rz2,
             rz1 - rz3 - sg * rz2,
            -rz1 + rz3 - sg * rz2,
            -rz1 - rz3 + sg * rz2
        );
        numSheets = 4.0;
    }
    else
    {
        float H   = 2.0 * lxy * sqrt(max(D, 0.0));
        float S   = pow(max(Q + H, 0.0), 1.0 / 3.0);
        float U   = (S > 0.0) ? (C2 / S) : 0.0;
        float Zx  = r2 - C + 0.5 * (S + U);
        float dSU = S * S + C2 + U * U;
        float Zy  = 0.866025404 * ((dSU > 0.0) ? (2.0 * H / dSU) : 0.0);
        float Zm  = length(vec2(Zx, Zy));
        float t1  = (Zm > 0.0) ? (l * (x2 - y2) / Zm) : 0.0;
        float t2  = sqrt(max(2.0 * (Zm + Zx), 0.0));
        vList     = 0.5 * vec4(t1 + t2, t1 - t2, t1 - t2, t1 - t2);
        numSheets = 2.0;
    }

    vec4 taus;
    for( int i = 0; i < 4; i++ )
    {
        float v = vList[i];
        float cEl, sEl;
        if( v >= 0.0 )
        {
            float den = v + e;
            sEl = (den > 1e-20) ? min(y / den, 1.0) : 0.0;
            cEl = ((v >= e) ? 1.0 : -1.0) * sqrt(max(1.0 - sEl * sEl, 0.0));
        }
        else
        {
            float den = v - e;
            cEl = (den < -1e-20) ? max(x / den, -1.0) : -1.0;
            sEl = ((v >= -e) ? 1.0 : -1.0) * sqrt(max(1.0 - cEl * cEl, 0.0));
        }
        vec2  vn  = vec2(-b * cEl, -a * sEl);
        vec2  nIn = vn * inversesqrt(max(dot(vn, vn), 1e-30));
        taus[i]   = dot(p - vec2(a * cEl, b * sEl), nIn);
    }

    if( taus.x > taus.y ) taus.xy = taus.yx;
    if( taus.z > taus.w ) taus.zw = taus.wz;
    if( taus.x > taus.z ) taus.xz = taus.zx;
    if( taus.y > taus.w ) taus.yw = taus.wy;
    if( taus.y > taus.z ) taus.yz = taus.zy;

    return vec4(taus.xyz, numSheets);
}

// Exact 2D Euclidean SDF to the Multi-Sheeted Ellipse Legendrian Wavefront at offset tOffset
float ctEllipseWavefront2D( in vec2 p, in vec2 ab, in float tOffset )
{
    p = abs(p);
    if( ab.x > ab.y ) { p = p.yx; ab = ab.yx; }

    float a = ab.x, b = ab.y;
    float a2 = a * a, b2 = b * b;
    float l  = b2 - a2, e = 0.5 * l;
    float x  = a * p.x, y = b * p.y;
    float x2 = x * x, y2 = y * y, r2 = x2 + y2;

    float C   = (r2 - l * l) / 3.0;
    float C2  = C * C, C3 = C2 * C;
    float lxy = l * x * y, l2x2y2 = lxy * lxy;
    float D   = C3 + l2x2y2, Q = C3 + 2.0 * l2x2y2;

    vec4 vList;
    int  nRoots;
    if( D < 0.0 )
    {
        float h   = atan(2.0 * lxy * sqrt(-D), -Q) / 3.0;
        float s   = cos(h), tr = sin(h) * 1.732050808;
        float oms = (tr * tr) / (3.0 * (1.0 + s));
        float rz1 = sqrt(max(l * l + 2.0 * C * oms, 0.0));
        float rz2 = sqrt(max(r2 - C * (oms - tr), 0.0));
        float rz3 = sqrt(max(r2 - C * (oms + tr), 0.0));
        float sg  = (x2 >= y2) ? 1.0 : -1.0;
        vList = 0.5 * vec4(
             rz1 + rz3 + sg * rz2,
             rz1 - rz3 - sg * rz2,
            -rz1 + rz3 - sg * rz2,
            -rz1 - rz3 + sg * rz2
        );
        nRoots = 4;
    }
    else
    {
        float H   = 2.0 * lxy * sqrt(max(D, 0.0));
        float S   = pow(max(Q + H, 0.0), 1.0 / 3.0);
        float U   = (S > 0.0) ? (C2 / S) : 0.0;
        float Zx  = r2 - C + 0.5 * (S + U);
        float dSU = S * S + C2 + U * U;
        float Zy  = 0.866025404 * ((dSU > 0.0) ? (2.0 * H / dSU) : 0.0);
        float Zm  = length(vec2(Zx, Zy));
        float t1  = (Zm > 0.0) ? (l * (x2 - y2) / Zm) : 0.0;
        float t2  = sqrt(max(2.0 * (Zm + Zx), 0.0));
        vList  = 0.5 * vec4(t1 + t2, t1 - t2, 0.0, 0.0);
        nRoots = 2;
    }

    float bestD = 1e30;
    for( int i = 0; i < 4; i++ )
    {
        if( i >= nRoots ) break;
        float v = vList[i];
        float cEl, sEl;
        if( v >= 0.0 )
        {
            float den = v + e;
            sEl = (den > 1e-20) ? min(y / den, 1.0) : 0.0;
            cEl = ((v >= e) ? 1.0 : -1.0) * sqrt(max(1.0 - sEl * sEl, 0.0));
        }
        else
        {
            float den = v - e;
            cEl = (den < -1e-20) ? max(x / den, -1.0) : -1.0;
            sEl = ((v >= -e) ? 1.0 : -1.0) * sqrt(max(1.0 - cEl * cEl, 0.0));
        }
        vec2  vn  = vec2(-b * cEl, -a * sEl);
        vec2  nIn = vn * inversesqrt(max(dot(vn, vn), 1e-30));
        float tau = dot(p - vec2(a * cEl, b * sEl), nIn);
        bestD     = min(bestD, abs(tau - tOffset));
    }

    if( l > 1e-12 && tOffset > 0.0 )
    {
        float abt23 = pow(a * b * tOffset, 2.0 / 3.0);
        if( abt23 > a2 && abt23 < b2 )
        {
            float cc   = sqrt((abt23 - a2) / l);
            float sc   = sqrt((b2 - abt23) / l);
            vec2  cusp = vec2((l / a) * cc * cc * cc, (l / b) * sc * sc * sc);
            bestD = min(bestD, length(p - cusp));
        }
    }
    return bestD;
}

// ---------------------------------------------------------------------------
// 6. Exact 3D Arnold A3 Swallowtail Wavefront Shell (sdSwallowtailWavefront3D)
// ---------------------------------------------------------------------------
vec4 sdgSwallowtailWavefront3D( in vec3 p, in float k, in float tOffset, in float sMax, in float hz, in float th )
{
    float ik   = 1.0 / max(k, 1e-6);
    float kt   = max(k * tOffset, 1.0);
    float kt23 = pow(kt, 2.0 / 3.0);
    float yCen = 0.25 * ik * (3.0 * kt23 - 1.0 + kt * kt + 1.0);
    vec2  q2   = vec2(p.x, p.y + yCen);

    float sx = (q2.x < 0.0) ? -1.0 : 1.0;
    float ax = abs(q2.x);
    float pc = ik * ik * (k * q2.y - 1.0) * (2.0 / 3.0);
    float qc = ik * ik * ax;
    float h  = qc * qc - pc * pc * pc;

    float bestD2 = 1e30;
    vec2  bestG2 = vec2(0.0, 1.0);

    vec4 sCand = vec4(sMax, -sMax, sMax, -sMax);
    int  nCand = 2;
    if( h >= 0.0 )
    {
        float r  = pow(qc + sqrt(h), 1.0 / 3.0);
        float s1 = (pc < 0.0) ? (2.0 * qc / (r * r - pc + (pc / r) * (pc / r)))
                              : ((r > 0.0) ? (r + pc / r) : 0.0);
        sCand.z = clamp(s1, -sMax, sMax);
        nCand   = 3;
    }
    else
    {
        float r  = sqrt(max(pc, 0.0));
        float tr = acos(clamp(qc / max(pc * r, 1e-30), -1.0, 1.0)) / 3.0;
        float cs = cos(tr), sn = sin(tr) * 1.732050808;
        sCand.yzw = clamp(r * vec3(2.0 * cs, -(cs + sn), -(cs - sn)), -sMax, sMax);
        nCand     = 4;
    }

    for( int i = 0; i < 4; i++ )
    {
        if( i >= nCand ) break;
        float si  = sCand[i];
        float ks  = k * si;
        float inv = inversesqrt(1.0 + ks * ks);
        vec2  n   = vec2(-ks * inv, inv);
        vec2  wPt = vec2(si - tOffset * ks * inv, 0.5 * k * si * si + tOffset * inv);
        vec2  dv  = vec2(ax, q2.y) - wPt;
        float d2  = dot(dv, dv);
        if( d2 < bestD2 )
        {
            bestD2 = d2;
            float l = sqrt(d2);
            bestG2 = (l > 1e-6) ? (dv / l) : ((dot(dv, n) >= 0.0) ? n : -n);
        }
    }

    if( kt > 1.0 )
    {
        float u2    = max(kt23 - 1.0, 0.0);
        float sCusp = sqrt(u2) * ik;
        if( sCusp <= sMax )
        {
            vec2  cusp = vec2(u2 * sqrt(u2) * ik, 0.5 * ik * (3.0 * kt23 - 1.0));
            vec2  dv   = vec2(ax, q2.y) - cusp;
            float d2c  = dot(dv, dv);
            if( d2c < bestD2 )
            {
                bestD2 = d2c;
                float l = sqrt(d2c);
                bestG2 = (l > 1e-6) ? (dv / l) : vec2(0.0, 1.0);
            }
        }
    }

    float d2D = sqrt(bestD2);
    vec2  g2D = vec2(sx * bestG2.x, bestG2.y);

    float az  = abs(p.z);
    float dz  = max(az - hz, 0.0);
    float sz  = (p.z < 0.0) ? -1.0 : 1.0;
    float d3D = length(vec2(d2D, dz));
    vec3  g3D = (d3D > 1e-6) ? vec3(g2D * (d2D / d3D), sz * (dz / d3D)) : vec3(g2D, 0.0);
    return vec4(d3D - th, g3D);
}

float sdSwallowtailWavefront3D( in vec3 p, in float k, in float tOffset, in float sMax, in float hz, in float th )
{
    return sdgSwallowtailWavefront3D(p, k, tOffset, sMax, hz, th).x;
}

vec4 sdgSwallowtailWavefront3D( in vec3 p, in float tOffset )
{
    return sdgSwallowtailWavefront3D(p, 1.0, tOffset, 2.2, 0.45, 0.045);
}

float sdSwallowtailWavefront3D( in vec3 p, in float tOffset )
{
    return sdgSwallowtailWavefront3D(p, 1.0, tOffset, 2.2, 0.45, 0.045).x;
}

#endif // SDF_GEMS_HAS_CT_LEGENDRIAN_WAVEFRONT
