import numpy as np
from PIL import Image, ImageDraw, ImageFont
import os

def smoothstep(edge0, edge1, x):
    t = np.clip((x - edge0) / (edge1 - edge0 + 1e-12), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)

def render_2d_iq(sdf_fn, bounds=(-1.5, 1.5), res=512, freq=120.0):
    """
    Renders a 2D SDF using Inigo Quilez's signature Shadertoy color scheme,
    plus returns the Eikonal gradient norm ||grad(d)|| map.
    sdf_fn takes (x, y) arrays of shape (H, W) or p of shape (..., 2).
    """
    lin = np.linspace(bounds[0], bounds[1], res)
    x, y = np.meshgrid(lin, -lin)  # y up
    p = np.stack([x, y], axis=-1)
    d = sdf_fn(p)

    # Central difference gradient for Eikonal check
    eps = (bounds[1] - bounds[0]) / res
    dx = (sdf_fn(p + np.array([eps, 0.0])) - sdf_fn(p - np.array([eps, 0.0]))) / (2.0 * eps)
    dy = (sdf_fn(p + np.array([0.0, eps])) - sdf_fn(p - np.array([0.0, eps]))) / (2.0 * eps)
    grad_norm = np.sqrt(dx * dx + dy * dy)

    # IQ coloring: orange outside (0.9, 0.6, 0.3), blue inside (0.65, 0.85, 1.0)
    col = np.where(
        d[..., None] > 0.0,
        np.array([0.9, 0.6, 0.3])[None, None, :],
        np.array([0.65, 0.85, 1.0])[None, None, :]
    )
    col = col * (1.0 - np.exp(-6.0 * np.abs(d))[..., None])
    col = col * (0.8 + 0.2 * np.cos(freq * d)[..., None])

    # White zero-isoline
    px = (bounds[1] - bounds[0]) / res * 1.8
    edge = 1.0 - smoothstep(0.0, px, np.abs(d))[..., None]
    col = col * (1.0 - edge) + edge * 1.0

    img_rgb = np.clip(col * 255.0, 0, 255).astype(np.uint8)
    return img_rgb, d, grad_norm


def render_eikonal_error_map(grad_norm, d):
    """
    Visualizes Eikonal deviation ||grad(d)|| - 1.0:
    - Dark grey/green: ||grad(d)|| == 1.0 (exact Euclidean SDF!)
    - Blue/Cyan: ||grad(d)|| < 1.0 (underestimate / conservative lower bound)
    - Red/Yellow: ||grad(d)|| > 1.0 (overestimate / raymarcher violation!)
    """
    err = grad_norm - 1.0
    H, W = grad_norm.shape
    col = np.zeros((H, W, 3), dtype=np.float64)

    # Base dark background
    col[:] = np.array([0.08, 0.10, 0.12])

    # Underestimation (err < 0): blue-cyan
    under = np.clip(-err / 0.5, 0.0, 1.0)
    col += under[..., None] * np.array([0.0, 0.55, 0.90])

    # Overestimation (err > 0): bright red-yellow warning!
    over = np.clip(err / 0.3, 0.0, 1.0)
    col += over[..., None] * np.array([1.0, 0.15, 0.05])

    # Exact region (|err| < 0.01): subtle emerald tint
    exact = np.exp(- (err / 0.015) ** 2)
    col = col * (1.0 - 0.3 * exact[..., None]) + exact[..., None] * np.array([0.12, 0.28, 0.18])

    # Zero contour overlay
    zero_line = np.exp(- (d / 0.012) ** 2)
    col = col * (1.0 - zero_line[..., None]) + zero_line[..., None] * np.array([1.0, 1.0, 1.0])

    return np.clip(col * 255.0, 0, 255).astype(np.uint8)


def render_3d_sdf(
    sdf_fn,
    cam_pos=(2.2, 1.7, 2.4),
    cam_target=(0.0, 0.0, 0.0),
    res=(400, 400),
    fov=1.8,
    max_steps=100,
    max_dist=12.0,
    show_floor_slice=True,
    floor_y=-0.85,
):
    """
    Vectorized 3D SDF Raymarcher producing both:
    1. Shaded 3D image with IQ's signature floor SDF isoline slice, soft shadows, and AO
    2. Step-count heatmap image
    sdf_fn(p) takes p of shape (N, 3) and returns (N,) float array.
    """
    W, H = res
    u = (np.linspace(0.5, W - 0.5, W) - 0.5 * W) / H
    v = (0.5 * H - np.linspace(0.5, H - 0.5, H)) / H
    uu, vv = np.meshgrid(u, v)

    ro = np.array(cam_pos, dtype=np.float64)
    ta = np.array(cam_target, dtype=np.float64)
    fw = ta - ro
    fw /= np.linalg.norm(fw)
    up_hint = np.array([0.0, 1.0, 0.0])
    rt = np.cross(fw, up_hint)
    rt /= np.linalg.norm(rt)
    up = np.cross(rt, fw)

    rd = uu[..., None] * rt + vv[..., None] * up + fov * fw
    rd /= np.linalg.norm(rd, axis=-1, keepdims=True)
    rd_flat = rd.reshape(-1, 3)
    N = rd_flat.shape[0]

    t = np.zeros(N, dtype=np.float64)
    steps = np.zeros(N, dtype=np.int32)
    hit = np.zeros(N, dtype=bool)
    active = np.ones(N, dtype=bool)

    for i in range(max_steps):
        idx = np.where(active)[0]
        if idx.size == 0:
            break
        pts = ro[None, :] + t[idx, None] * rd_flat[idx]
        d = sdf_fn(pts)
        steps[idx] = i + 1

        newly_hit = np.abs(d) < 0.0008
        hit[idx[newly_hit]] = True
        t[idx] += d
        escaped = t[idx] > max_dist
        active[idx[newly_hit | escaped]] = False

    # Base sky/background
    col = np.zeros((N, 3), dtype=np.float64)
    col[:] = np.array([0.09, 0.11, 0.14]) - 0.04 * rd_flat[:, 1:2]

    # Floor plane intersection (for rays that didn't hit the 3D object closer than floor)
    if show_floor_slice:
        with np.errstate(divide='ignore', invalid='ignore'):
            t_floor = (floor_y - ro[1]) / rd_flat[:, 1]
        hit_floor = (rd_flat[:, 1] < -1e-4) & (t_floor > 0.0) & (t_floor < max_dist) & ((~hit) | (t_floor < t))
        f_idx = np.where(hit_floor)[0]
        if f_idx.size > 0:
            f_pts = ro[None, :] + t_floor[f_idx, None] * rd_flat[f_idx]
            fd = sdf_fn(f_pts)
            fcol = np.where(
                fd[:, None] > 0.0,
                np.array([0.9, 0.6, 0.3])[None, :],
                np.array([0.65, 0.85, 1.0])[None, :]
            )
            fcol = fcol * (1.0 - np.exp(-6.0 * np.abs(fd))[:, None])
            fcol = fcol * (0.8 + 0.2 * np.cos(80.0 * fd)[:, None])
            edge = 1.0 - smoothstep(0.0, 0.015, np.abs(fd))[:, None]
            fcol = fcol * (1.0 - edge) + edge * 1.0
            # Distance fog on floor
            fog = np.exp(-0.04 * t_floor[f_idx])[:, None]
            col[f_idx] = fcol * fog + col[f_idx] * (1.0 - fog)
            # Update hit mask if floor occluded object
            hit[f_idx] = False

    # Shade hit points on the 3D primitive
    h_idx = np.where(hit)[0]
    if h_idx.size > 0:
        pos = ro[None, :] + t[h_idx, None] * rd_flat[h_idx]
        # Tetrahedron normal (4 evaluations!)
        eps = 0.001
        k0 = np.array([1.0, -1.0, -1.0])
        k1 = np.array([-1.0, -1.0, 1.0])
        k2 = np.array([-1.0, 1.0, -1.0])
        k3 = np.array([1.0, 1.0, 1.0])
        nor = (
            k0 * sdf_fn(pos + eps * k0)[:, None] +
            k1 * sdf_fn(pos + eps * k1)[:, None] +
            k2 * sdf_fn(pos + eps * k2)[:, None] +
            k3 * sdf_fn(pos + eps * k3)[:, None]
        )
        nor /= (np.linalg.norm(nor, axis=-1, keepdims=True) + 1e-12)

        # Ambient occlusion (5 steps)
        occ = np.zeros(h_idx.size, dtype=np.float64)
        sca = 1.0
        for aoi in range(1, 6):
            hr = 0.01 + 0.11 * (aoi - 1) / 4.0
            dd = sdf_fn(pos + hr * nor)
            occ += (hr - dd) * sca
            sca *= 0.95
        ao = np.clip(1.0 - 2.5 * occ, 0.0, 1.0)

        # Lighting
        lig = np.array([0.55, 0.72, 0.42])
        lig /= np.linalg.norm(lig)
        hal = lig[None, :] - rd_flat[h_idx]
        hal /= (np.linalg.norm(hal, axis=-1, keepdims=True) + 1e-12)

        dif = np.clip(np.sum(nor * lig[None, :], axis=-1), 0.0, 1.0)
        # Fast soft shadow (16 steps)
        sh = np.ones(h_idx.size, dtype=np.float64)
        st = np.full(h_idx.size, 0.02)
        for _ in range(16):
            sh_pts = pos + st[:, None] * lig[None, :]
            sd = sdf_fn(sh_pts)
            sh = np.minimum(sh, np.clip(8.0 * sd / st, 0.0, 1.0))
            st += np.clip(sd, 0.02, 0.25)
        dif *= sh

        spe = np.power(np.clip(np.sum(nor * hal, axis=-1), 0.0, 1.0), 16.0) * dif
        amb = np.clip(0.5 + 0.5 * nor[:, 1], 0.0, 1.0) * ao
        bac = np.clip(np.sum(nor * np.array([-lig[0], 0.0, -lig[2]])[None, :], axis=-1), 0.0, 1.0) * np.clip(1.0 - pos[:, 1], 0.0, 1.0) * ao

        # Base material color with normal-tinted highlights
        mat = 0.28 + 0.10 * nor
        lin_light = (
            1.55 * dif[:, None] * np.array([1.00, 0.85, 0.68]) +
            0.40 * amb[:, None] * np.array([0.40, 0.62, 1.00]) +
            0.22 * bac[:, None] * np.array([0.25, 0.25, 0.30])
        )
        scol = mat * lin_light + 0.28 * spe[:, None] * np.array([1.0, 0.95, 0.85])
        # ACES-like tone mapping to prevent blown-out highlights and preserve facet contrast
        scol = (scol * (2.51 * scol + 0.03)) / (scol * (2.43 * scol + 0.59) + 0.14)
        col[h_idx] = scol

    # Gamma correction
    col = np.power(np.clip(col, 0.0, 1.0), 0.4545)
    img_rgb = np.clip(col.reshape(H, W, 3) * 255.0, 0, 255).astype(np.uint8)

    # Step count heatmap (blue -> green -> yellow -> red)
    s_norm = np.clip(steps.reshape(H, W) / float(max_steps), 0.0, 1.0)
    heat = np.zeros((H, W, 3), dtype=np.float64)
    heat[..., 0] = np.clip(2.0 * s_norm - 0.4, 0.0, 1.0)
    heat[..., 1] = np.sin(np.pi * s_norm)
    heat[..., 2] = np.clip(1.0 - 2.0 * s_norm, 0.0, 1.0)
    heat_rgb = np.clip(heat * 255.0, 0, 255).astype(np.uint8)

    return img_rgb, heat_rgb, steps.reshape(H, W)
