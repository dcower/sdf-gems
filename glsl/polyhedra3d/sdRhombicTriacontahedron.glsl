// ============================================================================
// sdRhombicTriacontahedron — Exact 1-Sqrt 3D Catalan Rhombic Triacontahedron
// @author Dillon Cower (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, max error `2.22e-15`)
// ============================================================================

float sdRhombicTriacontahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for (int i = 0; i < 3; i++) {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);
    if (p.x + 1.61803399 * p.y <= 0.61803399 * r) return p.z - r;
    float t = clamp((p.y - 1.61803399 * (p.x - 0.61803399 * r)) * 0.27639320, 0.0, 0.38196601 * r);
    return length(vec3(p.x - 0.61803399 * r + 1.61803399 * t, p.y - t, p.z - r));
}
