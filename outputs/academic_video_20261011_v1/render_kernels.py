"""CPU-only display kernels adapted from the existing 20261007 analytic renderer.
White background, no invented ground surface. Does not run any physical solver.
"""
import math
import numpy as np
from numba import njit, prange

@njit(cache=True)
def shade(nx, ny, nz, vx, vy, vz, r, g, b, metal=False):
    # Fixed studio illumination in world coordinates.
    lx, ly, lz = -.436, -.545, .716
    diffuse = max(0., nx*lx+ny*ly+nz*lz)
    fill = max(0., nx*.68+ny*.35+nz*.64)
    hx, hy, hz = lx+vx, ly+vy, lz+vz
    hlen = math.sqrt(hx*hx+hy*hy+hz*hz)
    spec = max(0., (nx*hx+ny*hy+nz*hz)/hlen)
    facing = max(0., nx*vx+ny*vy+nz*vz)
    rim = (1-facing)**3
    if metal:
        power = .30 + .50*diffuse + .16*fill
        shine = 72*spec**35 + 18*rim
    else:
        power = .30 + .64*diffuse + .14*fill
        shine = 88*spec**65 + 12*rim
    return min(255., r*power+shine), min(255., g*power+shine), min(255., b*power+shine)

@njit(cache=True, parallel=True)
def base_scene(w, h, scale, target, right, up, view, origin, rot, radius, bore, height, cut):
    rgb = np.empty((h, w, 3), np.uint8)
    depth = np.full((h, w), 1.e9, np.float32)
    direction = -view
    dx = direction[0]*rot[0, 0]+direction[1]*rot[1, 0]
    dy = direction[0]*rot[0, 1]+direction[1]*rot[1, 1]
    dz = direction[2]
    for row in prange(h):
        v = (h/2-row-.5)/scale
        for col in range(w):
            u = (col+.5-w/2)/scale
            ox = target[0]+20*view[0]+u*right[0]+v*up[0]
            oy = target[1]+20*view[1]+u*right[1]+v*up[1]
            oz = target[2]+20*view[2]+u*right[2]+v*up[2]
            rr, gg, bb = 255., 255., 255.
            # Transform rays to the actual saved tool frame.
            qx, qy = ox-origin[0], oy-origin[1]
            ax = qx*rot[0, 0]+qy*rot[1, 0]
            ay = qx*rot[0, 1]+qy*rot[1, 1]
            az = oz-origin[2]
            best = 1.e9
            bn0, bn1, bn2 = 0., 0., 0.
            kind = 0
            # Finite cylindrical outside. Tool rotates, round outside does not.
            a = dx*dx+dy*dy
            b = 2*(ax*dx+ay*dy)
            c = ax*ax+ay*ay-radius*radius
            disc = b*b-4*a*c
            if disc >= 0 and a > 1.e-12:
                for sign in (-1., 1.):
                    t = (-b+sign*math.sqrt(disc))/(2*a)
                    z = az+t*dz
                    world_y = oy+t*direction[1]-origin[1]
                    if t > 0 and t < best and 0 <= z <= height and world_y >= cut:
                        x, y = ax+t*dx, ay+t*dy
                        best = t
                        bn0 = (x*rot[0, 0]+y*rot[0, 1])/radius
                        bn1 = (x*rot[1, 0]+y*rot[1, 1])/radius
                        bn2 = 0.
                        kind = 1
            # Top and underside, with the square opening left empty.
            if abs(dz) > 1.e-12:
                for zcap in (0., height):
                    t = (zcap-az)/dz
                    x, y = ax+t*dx, ay+t*dy
                    world_y = oy+t*direction[1]-origin[1]
                    if t > 0 and t < best and x*x+y*y <= radius*radius and (abs(x) >= bore or abs(y) >= bore) and world_y >= cut:
                        best = t
                        bn0, bn1, bn2 = 0., 0., (1. if zcap > 0 else -1.)
                        kind = 2
            # Four planar walls of the genuinely open, rotating square bore.
            for axis in range(2):
                da = dx if axis == 0 else dy
                aa = ax if axis == 0 else ay
                if abs(da) < 1.e-12:
                    continue
                for sign in (-1., 1.):
                    t = (sign*bore-aa)/da
                    other = ay+t*dy if axis == 0 else ax+t*dx
                    z = az+t*dz
                    world_y = oy+t*direction[1]-origin[1]
                    if t > 0 and t < best and abs(other) <= bore and 0 <= z <= height and world_y >= cut:
                        best = t
                        bn0, bn1, bn2 = -sign*rot[0, axis], -sign*rot[1, axis], 0.
                        kind = 3
            # Display-only world-y section plane; no material particles clipped.
            if cut > -radius and abs(direction[1]) > 1.e-12:
                t = (origin[1]+cut-oy)/direction[1]
                x, y, z = ax+t*dx, ay+t*dy, az+t*dz
                if t > 0 and t < best and 0 <= z <= height and x*x+y*y <= radius*radius and (abs(x) >= bore or abs(y) >= bore):
                    best = t
                    bn0, bn1, bn2 = 0., -1., 0.
                    kind = 4
            if best < depth[row, col]:
                if kind == 4:
                    cr, cg, cb = 192., 202., 211.
                elif kind == 3:
                    cr, cg, cb = 185., 197., 207.
                else:
                    cr, cg, cb = 210., 217., 224.
                rr, gg, bb = shade(bn0, bn1, bn2, view[0], view[1], view[2], cr, cg, cb, True)
                depth[row, col] = best
            rgb[row, col, 0] = int(rr)
            rgb[row, col, 1] = int(gg)
            rgb[row, col, 2] = int(bb)
    return rgb, depth

@njit(cache=True)
def spheres(rgb, depth, positions, colors, scale, target, right, up, view, radius):
    h, w = depth.shape
    for i in range(len(positions)):
        p = positions[i]-target
        cx = w/2+(p[0]*right[0]+p[1]*right[1]+p[2]*right[2])*scale
        cy = h/2-(p[0]*up[0]+p[1]*up[1]+p[2]*up[2])*scale
        ct = 20-(p[0]*view[0]+p[1]*view[1]+p[2]*view[2])
        rad = radius*scale
        for y in range(max(0, int(cy-rad-1)), min(h, int(cy+rad+2))):
            b = (cy-y-.5)/scale
            for x in range(max(0, int(cx-rad-1)), min(w, int(cx+rad+2))):
                a = (x+.5-cx)/scale
                disc = radius*radius-a*a-b*b
                if disc <= 0:
                    continue
                root = math.sqrt(disc)
                t = ct-root
                if t >= depth[y, x]:
                    continue
                nx = (a*right[0]+b*up[0]+root*view[0])/radius
                ny = (a*right[1]+b*up[1]+root*view[1])/radius
                nz = (a*right[2]+b*up[2]+root*view[2])/radius
                r, g, bc = shade(nx, ny, nz, view[0], view[1], view[2], colors[i, 0], colors[i, 1], colors[i, 2])
                rgb[y, x, 0], rgb[y, x, 1], rgb[y, x, 2] = int(r), int(g), int(bc)
                depth[y, x] = t

@njit(cache=True)
def depth_lines(rgb, depth, points, color, alpha, scale, target, right, up, view, thickness):
    h, w = depth.shape
    for i in range(1, len(points)):
        pa, pb = points[i-1]-target, points[i]-target
        xa, ya = w/2+np.dot(pa, right)*scale, h/2-np.dot(pa, up)*scale
        xb, yb = w/2+np.dot(pb, right)*scale, h/2-np.dot(pb, up)*scale
        ta, tb = 20-np.dot(pa, view), 20-np.dot(pb, view)
        n = max(1, int(max(abs(xb-xa), abs(yb-ya))*1.5))
        for j in range(n+1):
            f = j/n
            x, y, t = int(xa+(xb-xa)*f), int(ya+(yb-ya)*f), ta+(tb-ta)*f
            for yy in range(max(0, y-thickness), min(h, y+thickness+1)):
                for xx in range(max(0, x-thickness), min(w, x+thickness+1)):
                    if t <= depth[yy, xx]+.014:
                        for channel in range(3):
                            rgb[yy, xx, channel] = int(rgb[yy, xx, channel]*(1-alpha)+color[channel]*alpha)
