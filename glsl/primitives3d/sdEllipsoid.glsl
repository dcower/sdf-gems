// ============================================================================
// sdEllipsoid — 3D Triaxial Ellipsoid (Exterior Dual-Support Bound)
// @author Inigo Quilez (V2 bound) — upgraded/fixed by Dillon Cower (sdf-gems, 2026)
// @source https://iquilezles.org/articles/ellipsoids/
// @exact  Conservative Exterior Bound (`0%` exterior overstep, clamped `lam >= 0`)
// ============================================================================

float sdEllipsoid( vec3 p, vec3 r )
{
    p = abs(p);
    vec3 r2 = r * r;
    vec3 rp = r * p;
    vec3 q  = p / r;
    vec3 q2 = q * q;
    float k0_2 = q2.x + q2.y + q2.z;
    vec3 c2 = q / r;
    vec3 c4 = c2 * c2;
    float k1_2 = c4.x + c4.y + c4.z;

    float lam = max(sqrt(k0_2) - 1.0, 0.0) * (k0_2 / max(k1_2, 1e-18));
    lam = max(lam, max(max(rp.x - r2.x, rp.y - r2.y), rp.z - r2.z));

    float sxy = q2.x + q2.y, dxy = c4.x + c4.y;
    if (sxy > 1.0) lam = max(lam, (sqrt(sxy) - 1.0) * sxy / max(dxy, 1e-18));
    float syz = q2.y + q2.z, dyz = c4.y + c4.z;
    if (syz > 1.0) lam = max(lam, (sqrt(syz) - 1.0) * syz / max(dyz, 1e-18));
    float sxz = q2.x + q2.z, dxz = c4.x + c4.z;
    if (sxz > 1.0) lam = max(lam, (sqrt(sxz) - 1.0) * sxz / max(dxz, 1e-18));

    vec3 rp2 = rp * rp;
    float rp_norm2 = rp2.x + rp2.y + rp2.z;
    if (rp_norm2 > 1e-18)
        lam = max(lam, sqrt(rp_norm2) - dot(r2, rp2) / rp_norm2);

    vec3 inv = 1.0 / (r2 + lam);
    vec3 u   = rp * inv;
    float g2 = dot(u, u);
    if (g2 > 1.0) {
        float g = sqrt(g2);
        lam += g2 * (g - 1.0) / max(dot(u * u, inv), 1e-18);
    }

    vec3 v = p / (r2 + lam);
    return (dot(p, v) - length(r * v)) / max(length(v), 1e-18);
}
