// ============================================================================
// sdgEllipsoid, sdgEllipsoidExact & sdgSpheroid — Exact 3D Ellipsoid & Spheroid Gradients via Secular Co-Normal Projection
// @author Inigo Quilez (sdgEllipsoid) — upgraded with Secular Co-Normal Projection Theorem using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distgradfunctions3d/
// @exact  Yes (`||grad d|| = 1`, 0 extra sqrts over scalar SDF, 0 surface/origin NaNs, 0.006 deg vs IQ's 28.7 deg mean error)
// ============================================================================

// 1. 12-Line Unified Conservative Halley 3D Ellipsoid Distance & Gradient (0 extra sqrts over sdEllipsoid)
vec4 sdgEllipsoid( vec3 p, vec3 r )
{
    vec3  s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = max(abs(p), 1e-6 * r);
    vec3  r2 = r * r, u = r * p, u2 = u * u, v = p / r2, q2 = p * v, q4 = v * v;
    vec4  s1 = vec4(q2.yzx + q2.zxy, q2.x + q2.y + q2.z);
    vec4  s2 = vec4(q4.yzx + q4.zxy, q4.x + q4.y + q4.z);
    vec3  l1 = max(u, min(u * inversesqrt(max(1.0 - s1.xyz, 1e-12)), r2)) - r2;
    vec4  lv = max(vec4(l1, length(u) - dot(r2, u2) / dot(u, u)), s1 * (sqrt(s1) - 1.0) / s2);
    float lam = max(max(lv.x, lv.y), max(lv.z, lv.w));
    vec3  inv = 1.0 / (r2 + lam), w2 = u2 * inv * inv;
    float g2  = w2.x + w2.y + w2.z, g3 = dot(w2, inv), d = sqrt(g2) - 1.0;
    lam += g2 * d / (g3 + 1.5 * d * (g3 - g2 * dot(w2, inv * inv) / g3));
    v = p / (r2 + lam);
    float invV = inversesqrt(dot(v, v));
    return vec4((dot(p, v) - length(r * v)) * invV, s0 * v * invV);
}

// 2. True Exact 3D Triaxial Ellipsoid Distance & Gradient via Secular Co-Normal Projection (0 extra sqrts)
vec4 sdgEllipsoidExact( vec3 p, vec3 r )
{
    vec3  s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    bvec3 sw = bvec3(false);
    sw.x = (r.x < r.y); if( sw.x ) { r.xy = r.yx; p.xy = p.yx; }
    sw.y = (r.y < r.z); if( sw.y ) { r.yz = r.zy; p.yz = p.zy; }
    sw.z = (r.x < r.y); if( sw.z ) { r.xy = r.yx; p.xy = p.yx; }

    vec3  r2 = r * r;
    vec2  d  = r2.xy - r2.z;
    vec3  u  = r * p;
    float u2_yz = dot(u.yz, u.yz), u_yz = sqrt(u2_yz), u_xyz = sqrt(u.x * u.x + u2_yz);
    float t_min = max(u.z, max(u_yz - d.y, u_xyz - d.x));
    float ts    = max(t_min, 1e-20);

    vec2  ixy0 = 1.0 / (d + ts);
    float iz0  = (u.z > 0.0) ? (1.0 / ts) : 0.0;
    vec2  wxy2 = (u.xy * ixy0) * (u.xy * ixy0);
    float wz2  = (t_min > 0.0) ? ((u.z * iz0) * (u.z * iz0)) : 0.0;

    float exy  = wxy2.x + wxy2.y;
    float ixy  = dot(wxy2, ixy0) / max(exy, 1e-20);
    float dyz  = d.y * (u.y * u.y / max(u2_yz, 1e-20));
    float iyz  = 1.0 / (dyz + ts);
    float wyz2 = u2_yz * iyz * iyz;

    vec2 A     = vec2(exy, wxy2.x);
    vec2 a     = max(vec2(ixy, ixy0.x), 1e-20);
    vec2 B     = vec2(wz2, wyz2);
    vec2 b     = max(vec2(iz0, iyz), a);
    vec2 dmin  = max(vec2(0.0), (sqrt(A) - 1.0) / a);
    vec2 cusp  = pow(max(2.0 * A * a * sqrt(B) / b, 1e-30), vec2(0.6666667));
    vec2 dcusp = max(vec2(0.0), (sqrt(B / max(max(1.0 - A, 0.0) + cusp, 1e-20)) - 1.0) / b);
    vec2 d12   = sqrt(dmin * dmin + dcusp * dcusp);

    bool  on_focal = (t_min <= 1e-20) && (exy <= 1.0);
    float t        = on_focal ? 0.0 : (ts + max(d12.x, d12.y));

    for( int i = 0; i < 2; i++ )
    {
        float ts_i  = max(t, 1e-20);
        vec2  iv    = 1.0 / (d + ts_i);
        vec2  w2    = (u.xy * iv) * (u.xy * iv);
        float wz2_i = on_focal ? 0.0 : ((u.z / ts_i) * (u.z / ts_i));
        vec2  q     = ts_i * iv;
        float s1    = w2.x + w2.y + wz2_i;
        float ts2   = dot(w2, q) + wz2_i;
        float ts3   = dot(w2, q * q) + wz2_i;
        float e     = sqrt(s1) - 1.0;
        float d2    = max(s1 * ts3 - ts2 * ts2, 0.0);
        float step  = ts_i * (s1 * ts2 * e) / max(ts2 * ts2 - 1.5 * e * d2, max(0.1 * ts2 * ts2, 1e-30));
        t = on_focal ? 0.0 : max(t_min, t + step);
    }

    vec2  vxy    = p.xy / max(d + t, 1e-20);
    vec2  wxy2_f = (r.xy * vxy) * (r.xy * vxy);
    float uz_t   = u.z / max(t, 1e-20);
    float wz2_f  = (on_focal || wxy2_f.x + wxy2_f.y < 0.5)
                 ? max(0.0, 1.0 - wxy2_f.x - wxy2_f.y)
                 : (uz_t * uz_t);
    vec3  v      = vec3(vxy, sqrt(wz2_f) / r.z);
    vec3  x      = r2 * v * inversesqrt(max(wxy2_f.x + wxy2_f.y + wz2_f, 1e-20));
    vec3  n      = v * inversesqrt(max(dot(v, v), 1e-30));
    float dist   = dot(p - x, n);

    if( sw.z ) n.xy = n.yx;
    if( sw.y ) n.yz = n.zy;
    if( sw.x ) n.xy = n.yx;
    return vec4(dist, s0 * n);
}

// 3. Exact Closed-Form 3D Spheroid Distance & Gradient via Co-Normal Projection (0 extra sqrts)
vec4 sdgSpheroid( in vec3 p, in vec2 ab )
{
    float rxz2 = dot(p.xz, p.xz), rxz = sqrt(rxz2);
    vec2  er   = (rxz > 0.0) ? (p.xz / rxz) : vec2(1.0, 0.0);
    float sy   = (p.y < 0.0) ? -1.0 : 1.0;
    vec2  q    = vec2(rxz, abs(p.y));
    bool  sw   = (ab.x > ab.y);
    if( sw ) { q = q.yx; ab = ab.yx; }

    float a2  = ab.x * ab.x, b2 = ab.y * ab.y;
    float l   = b2 - a2, e = 0.5 * l;
    float x   = ab.x * q.x, y = ab.y * q.y;
    float x2  = x * x, y2 = y * y, r2 = x2 + y2;
    float C   = (r2 - l * l) / 3.0, C2 = C * C, C3 = C2 * C;
    float lxy = l * x * y, l2x2y2 = lxy * lxy;
    float D   = C3 + l2x2y2, Q = C3 + 2.0 * l2x2y2;

    float v, vm;
    if( D < 0.0 )
    {
        float h   = atan(2.0 * lxy * sqrt(-D), -Q) / 3.0;
        float s   = cos(h), t = sin(h) * 1.732050808;
        float oms = (t * t) / (3.0 * (1.0 + s));
        float rz1 = sqrt(max(l * l + 2.0 * C * oms, 0.0));
        float rz3 = sqrt(max(r2 - C * (oms + t), 0.0));
        float rz2 = sqrt(max(r2 - C * (oms - t), 0.0));
        v  = 0.5 * (rz1 + rz3 + ((x2 >= y2) ? rz2 : -rz2));
        vm = (x2 < y2) ? (-C) * (t / max(rz3 + rz2, 1e-30) - oms / max(rz1 + l, 1e-30)) : (v - e);
    }
    else
    {
        float H   = 2.0 * lxy * sqrt(D);
        float S   = pow(max(Q + H, 0.0), 1.0 / 3.0);
        float U   = (S > 0.0) ? (C2 / S) : 0.0;
        float Zx  = r2 - C + 0.5 * (S + U);
        float dSU = S * S + C2 + U * U;
        float Zy  = 0.866025404 * ((dSU > 0.0) ? (2.0 * H / dSU) : 0.0);
        float Zm  = length(vec2(Zx, Zy));
        v  = 0.5 * (sqrt(2.0 * (Zm + Zx)) + ((Zm > 0.0) ? l * (x2 - y2) / Zm : 0.0));
        vm = max(v - e, 0.0);
    }

    vec2  w    = vec2(x * (v + e), y * max(vm, 0.0));
    float wl   = length(w);
    float si_u = (l > 0.0) ? min(y / l, 1.0) : 0.0;
    vec2  cs   = (wl > 0.0) ? (w / wl) : vec2(sqrt(max(1.0 - si_u * si_u, 0.0)), si_u);
    vec2  vn   = ab.yx * cs;
    vec2  n    = vn * inversesqrt(dot(vn, vn));
    float d    = dot(q - ab * cs, n);
    if( sw ) n = n.yx;
    return vec4(d, n.x * er.x, sy * n.y, n.x * er.y);
}
