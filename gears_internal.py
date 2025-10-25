
from math import sin, cos, tan, atan, asin, acos, pi, sqrt

from make_stl import triangulate_polyhedron, triangulate_prism

# later maybe??
#from mpmath import mp
#mp.dps = 500
#sin = mp.sin
#cos = mp.cos
#tan = mp.tan
#atan = mp.atan
#asin = mp.asin
#acos = mp.acos
#pi = mp.pi
#sqrt = mp.sqrt

# all angles in this file are in radians.

def interpolate_line(a, b, n, endpoints=True):
    ans = []
    if endpoints:
        ans.append(a)
    for k in range(1, n):
        ans.append(tuple((a[i]*(n-k) + b[i]*k)/n for i in range(3)))
    if endpoints:
        ans.append(b)
    return ans

def center(points):
    n = 0
    ax, ay, az = 0, 0, 0
    for (x, y, z) in points:
        n += 1
        ax += (x-ax)/n
        ay += (y-ay)/n
        az += (z-az)/n
    return ax, ay, az

def norm(x, y, z):
    return sqrt(x*x + y*y + z*z)

def project(points, center, dist):
    cx, cy, cz = center
    ans = []
    for (x, y, z) in points:
        v = norm(x-cx, y-cy, z-cz)
        x_ = cx + dist*(x-cx)/v
        y_ = cy + dist*(y-cy)/v
        z_ = cz + dist*(z-cz)/v
        ans.append((x_, y_, z_))
    return ans


#same behavior as openscad rotate(a=...) but rotates a list of points
def rotate(a, l):

    for i in range(len(l)):
        x, y, z = l[i]

        # rotate a[2] radians around z axis
        x, y = x*cos(a[2]) - y*sin(a[2]), x*sin(a[2]) + y*cos(a[2])

        #rotate a[1] radians around y axis
        z, x = z*cos(a[1]) - x*sin(a[1]), z*sin(a[1]) + x*cos(a[1])

        #rotate a[0] radians around x axis
        y, z = y*cos(a[0]) - z*sin(a[0]), y*sin(a[0]) + z*cos(a[0])
    
        l[i] = (x, y, z)

def translate(a, l):
    for i in range(len(l)):
        x, y, z = l[i]
        l[i] = x+a[0], y+a[1], z+a[2]

def rad(t):
    return t/180*pi

def deg(t):
    return t/pi*180


def sphere_ev(t0, t):
    return acos(cos(t)/cos(t0))/sin(t0) - acos(tan(t0)/tan(t));

def sph_to_cart(v):
    r, theta, phi = v
    return (r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta))

clearance = 0.05

def bevel_gear_data(modul, tooth_number, partial_cone_angle, tooth_width, pressure_angle, helix_angle, tooth_step):

    d_outside = modul*tooth_number
    r_outside = d_outside / 2
    rg_outside = r_outside / sin(partial_cone_angle)
    rg_inside = rg_outside - tooth_width
    r_inside = r_outside*rg_inside / rg_outside
    alpha_spur = atan(tan(pressure_angle)/cos(helix_angle))
    da_outside = d_outside + (modul * (2.2 if modul < 1 else 2)) * cos(partial_cone_angle)
    ra_outside = da_outside / 2
    c = modul / 6
    df_outside = d_outside - (modul +c) * 2 * cos(partial_cone_angle)
    rf_outside = df_outside / 2
    #rkf = rg_outside*sin(delta_f)
    delta_f = asin(rf_outside/rg_outside)
    delta_a = asin(ra_outside/rg_outside)
    delta_b = asin(cos(alpha_spur)*sin(partial_cone_angle))

    height_f = rg_outside*cos(delta_f)

    height_k = (rg_outside-tooth_width)/cos(partial_cone_angle)
    rk = (rg_outside-tooth_width)/sin(partial_cone_angle)
    #rfk = rk*height_k*tan(delta_f)/(rk+height_k*tan(delta_f))
    #height_fk = rk*height_k/(height_k*tan(delta_f)+rk)


    phi_r = sphere_ev(delta_b, partial_cone_angle)

    gamma_g = 2*atan(tooth_width*tan(helix_angle)/(2*rg_outside-tooth_width))
    gamma = 2*asin(rg_outside/r_outside*sin(gamma_g/2))

    tau = 2*pi/tooth_number

    mirrpoint = (pi*(1-clearance))/tooth_number+2*phi_r

    tooth_nw = []
    tooth_ne = []
    tooth_sw = []
    tooth_se = []



    start = delta_f
    step = (delta_a - delta_f)/tooth_step # check if this is good
    if delta_b > delta_f:
        flankpoint_under = 1*mirrpoint

        tooth_ne.append(sph_to_cart((rg_outside, delta_f, flankpoint_under)))
        tooth_se.append(sph_to_cart((rg_inside, delta_f, flankpoint_under+gamma)))
        tooth_sw.append(sph_to_cart((rg_inside, delta_f, mirrpoint-flankpoint_under+gamma)))
        tooth_nw.append(sph_to_cart((rg_outside, delta_f, mirrpoint-flankpoint_under)))

        start = delta_b
        step = (delta_a - delta_b)/tooth_step

    #for delta in range(start, delta_a, step):
    for i in range(tooth_step+1):
        delta = (start*(tooth_step - i) + delta_a*i) / tooth_step
        flankpoint_under = sphere_ev(delta_b, delta)

        tooth_nw.append(sph_to_cart((rg_outside, delta, flankpoint_under)))
        tooth_sw.append(sph_to_cart((rg_inside, delta, flankpoint_under+gamma)))
        tooth_se.append(sph_to_cart((rg_inside, delta, mirrpoint-flankpoint_under+gamma)))
        tooth_ne.append(sph_to_cart((rg_outside, delta, mirrpoint-flankpoint_under)))


    for pt_list in (tooth_nw, tooth_ne, tooth_sw, tooth_se):
        rotate([0,pi,0], pt_list)
        translate([0,0,height_f], pt_list)
        rotate([0,0,phi_r+pi/2*(1-clearance)/tooth_number], pt_list)

    return (tooth_nw, tooth_ne, tooth_sw, tooth_se, -tau)

# do translation/rotation next
def bevel_gear_assembly(modul, tooth_number, partial_cone_angle, tooth_width, bore, pressure_angle, helix_angle, tooth_step, flat_step):
    tooth_nw, tooth_ne, tooth_sw, tooth_se, tau = bevel_gear_data(modul, tooth_number, partial_cone_angle, tooth_width, pressure_angle, helix_angle, tooth_step)

    tooth_top = tooth_nw + interpolate_line(tooth_nw[-1], tooth_ne[-1], flat_step, endpoints=False) + tooth_ne[::-1]
    tooth_bottom = tooth_sw + interpolate_line(tooth_sw[-1], tooth_se[-1], flat_step, endpoints=False) + tooth_se[::-1]

    top_face, bottom_face = [], []

    ans = []

    # create top and bottom teeth faces, teeth open prisms

    i = 0
    while True:
        ans += triangulate_polyhedron(tooth_top)
        ans += triangulate_polyhedron(tooth_bottom, reverse=True)
        ans += triangulate_prism(tooth_top, tooth_bottom, closed=False)

        if len(top_face) > 0:
            top_line = interpolate_line(top_face[-1], tooth_top[0], flat_step, endpoints=True)
            bottom_line = interpolate_line(bottom_face[-1], tooth_bottom[0], flat_step, endpoints=True)
            ans += triangulate_prism(top_line, bottom_line, closed=False)
            top_face.pop()
            bottom_face.pop()
            top_face += top_line
            bottom_face += bottom_line
            top_face.append(tooth_top[-1])
            bottom_face.append(tooth_bottom[-1])
        else:
            top_face = [tooth_top[0], tooth_top[-1]]
            bottom_face = [tooth_bottom[0], tooth_bottom[-1]]

        i += 1
        if i == tooth_number:
            break
        rotate((0, 0, tau), tooth_top)
        rotate((0, 0, tau), tooth_bottom)

    top_line = interpolate_line(top_face[-1], top_face[0], flat_step, endpoints=True)
    bottom_line = interpolate_line(bottom_face[-1], bottom_face[0], flat_step, endpoints=True)
    ans += triangulate_prism(top_line, bottom_line, closed=False)

    top_line.pop()
    bottom_line.pop()
    top_face.pop()
    bottom_face.pop()
    top_face += top_line
    bottom_face += bottom_line

    top_center = center(top_face)
    bottom_center = center(bottom_face)

    if bore == 0:
        ans += triangulate_polyhedron(top_face, top_center)
        ans += triangulate_polyhedron(bottom_face, bottom_center, reverse=True)

    else:
        top_bore, bottom_bore = [], []
        for k in range(len(top_face)):
            top_bore.append((-cos(2*pi*k/len(top_face)), sin(2*pi*k/len(top_face)), top_face[k][2]))
            bottom_bore.append((-cos(2*pi*k/len(top_face)), sin(2*pi*k/len(top_face)), bottom_face[k][2]))


        ans += triangulate_prism(top_bore, top_face, closed=True)
        ans += triangulate_prism(bottom_bore, top_bore, closed=True)
        ans += triangulate_prism(bottom_face, bottom_bore, closed=True)

    return ans

def bevel_herringbone_gear_data(modul, tooth_number, partial_cone_angle, tooth_width, pressure_angle, helix_angle, tooth_step):

    tooth_width = tooth_width / 2
    d_outside = modul * tooth_number
    r_outside = d_outside / 2
    rg_outside = r_outside/sin(partial_cone_angle)
    c = modul / 6
    df_outside = d_outside - (modul +c) * 2 * cos(partial_cone_angle)
    rf_outside = df_outside / 2
    delta_f = asin(rf_outside/rg_outside)
    height_f = rg_outside*cos(delta_f)

    gamma_g = 2*atan(tooth_width*tan(helix_angle)/(2*rg_outside-tooth_width))
    gamma = 2*asin(rg_outside/r_outside*sin(gamma_g/2))

    height_k = (rg_outside-tooth_width)/cos(partial_cone_angle)
    rk = (rg_outside-tooth_width)/sin(partial_cone_angle)
    rfk = rk*height_k*tan(delta_f)/(rk+height_k*tan(delta_f))
    height_fk = rk*height_k/(height_k*tan(delta_f)+rk)

    modul_inside = modul*(1-tooth_width/rg_outside)
    
    # I think the -1 degree correction added here is a band-aid to the real issue
    # of the gears not meshing.
    # For now I address this by averaging the two different middle gears instead.
    # Later I would like to try to address this differently:
    # I think it's a math error, not due to floating points -- the mean square
    # difference between the two gears didn't change at all when I tried using
    # more precise floating points with mpmath.

    lower_cone_angle = partial_cone_angle #- pi/180

    tooth_aw, tooth_ae, tooth_bw_upper, tooth_be_upper, tau = bevel_gear_data( \
            modul, tooth_number, lower_cone_angle, tooth_width, pressure_angle, helix_angle, tooth_step)

    tooth_bw_lower, tooth_be_lower, tooth_cw, tooth_ce, tau = bevel_gear_data( \
            modul_inside, tooth_number, partial_cone_angle, tooth_width, pressure_angle, -helix_angle, tooth_step)

    for pt_list in (tooth_bw_lower, tooth_be_lower, tooth_ce, tooth_cw):
        rotate([0, 0, -gamma], pt_list)
        translate([0, 0, height_f - height_fk], pt_list)

    tooth_bw, tooth_be = [], []
    for i in range(len(tooth_be_upper)):
        tooth_bw.append(center((tooth_bw_upper[i], tooth_bw_lower[i])))
        tooth_be.append(center((tooth_be_upper[i], tooth_be_lower[i])))

    return tooth_aw, tooth_ae, tooth_bw, tooth_be, tooth_cw, tooth_ce, tau


def bevel_herringbone_gear_assembly(modul, tooth_number, partial_cone_angle, tooth_width, bore, pressure_angle, helix_angle, tooth_step, flat_step):
    tooth_aw, tooth_ae, tooth_bw, tooth_be, tooth_cw, tooth_ce, tau = bevel_herringbone_gear_data(modul, tooth_number, partial_cone_angle, tooth_width, pressure_angle, helix_angle, tooth_step)

    tooth_a = tooth_aw + interpolate_line(tooth_aw[-1], tooth_ae[-1], flat_step, endpoints=False) + tooth_ae[::-1]
    tooth_b = tooth_bw + interpolate_line(tooth_bw[-1], tooth_be[-1], flat_step, endpoints=False) + tooth_be[::-1]
    tooth_c = tooth_cw + interpolate_line(tooth_cw[-1], tooth_ce[-1], flat_step, endpoints=False) + tooth_ce[::-1]

    a_face, b_face, c_face = [], [], []

    ans = []

    # create top and bottom teeth faces, teeth open prisms

    i = 0
    while True:
        ans += triangulate_polyhedron(tooth_a)
        ans += triangulate_polyhedron(tooth_c, reverse=True)
        ans += triangulate_prism(tooth_a, tooth_b, closed=False)
        ans += triangulate_prism(tooth_b, tooth_c, closed=False)

        if len(a_face) > 0:
            a_line = interpolate_line(a_face[-1], tooth_a[0], flat_step, endpoints=True)
            b_line = interpolate_line(b_face[-1], tooth_b[0], flat_step, endpoints=True)
            c_line = interpolate_line(c_face[-1], tooth_c[0], flat_step, endpoints=True)
            ans += triangulate_prism(a_line, b_line, closed=False)
            ans += triangulate_prism(b_line, c_line, closed=False)
            a_face.pop()
            b_face.pop()
            c_face.pop()
            a_face += a_line
            b_face += b_line
            c_face += c_line
            a_face.append(tooth_a[-1])
            b_face.append(tooth_b[-1])
            c_face.append(tooth_c[-1])
        else:
            a_face = [tooth_a[0], tooth_a[-1]]
            b_face = [tooth_b[0], tooth_b[-1]]
            c_face = [tooth_c[0], tooth_c[-1]]

        i += 1
        if i == tooth_number:
            break
        rotate((0, 0, tau), tooth_a)
        rotate((0, 0, tau), tooth_b)
        rotate((0, 0, tau), tooth_c)

    a_line = interpolate_line(a_face[-1], a_face[0], flat_step, endpoints=True)
    b_line = interpolate_line(b_face[-1], b_face[0], flat_step, endpoints=True)
    c_line = interpolate_line(c_face[-1], c_face[0], flat_step, endpoints=True)
    ans += triangulate_prism(a_line, b_line, closed=False)
    ans += triangulate_prism(b_line, c_line, closed=False)

    a_line.pop()
    c_line.pop()
    a_face.pop()
    c_face.pop()
    a_face += a_line
    c_face += c_line

    a_center = center(a_face)
    c_center = center(c_face)

    if bore == 0:
        ans += triangulate_polyhedron(a_face, a_center)
        ans += triangulate_polyhedron(c_face, c_center, reverse=True)

    else:
        top_bore, bottom_bore = [], []
        for k in range(len(a_face)):
            top_bore.append((-cos(2*pi*k/len(a_face)), sin(2*pi*k/len(a_face)), a_face[k][2]))
            bottom_bore.append((-cos(2*pi*k/len(a_face)), sin(2*pi*k/len(a_face)), c_face[k][2]))


        ans += triangulate_prism(top_bore, a_face, closed=True)
        ans += triangulate_prism(bottom_bore, top_bore, closed=True)
        ans += triangulate_prism(c_face, bottom_bore, closed=True)

    return ans


def bevel_gear_pair_assembly(modul, gear_teeth, pinion_teeth, axis_angle, tooth_width, gear_bore, pinion_bore, pressure_angle, helix_angle, together_built, tooth_step, flat_step):


    r_gear = modul*gear_teeth/2
    delta_gear = atan(sin(axis_angle)/(pinion_teeth/gear_teeth+cos(axis_angle)))
    delta_pinion = atan(sin(axis_angle)/(gear_teeth/pinion_teeth+cos(axis_angle)))
    rg = r_gear/sin(delta_gear)
    c = modul / 6
    df_pinion = 2*rg*delta_pinion - 2 * (modul + c)
    rf_pinion = df_pinion / 2
    delta_f_pinion = rf_pinion/(pi*rg) * pi
    rkf_pinion = rg*sin(delta_f_pinion)
    height_f_pinion = rg*cos(delta_f_pinion)

    df_gear = 2*rg*delta_gear - 2 * (modul + c)
    rf_gear = df_gear / 2
    delta_f_gear = rf_gear/rg
    rkf_gear = rg*sin(delta_f_gear)
    height_f_gear = rg*cos(delta_f_gear)

    gear_1 = bevel_gear_assembly(modul, gear_teeth, delta_gear, tooth_width, gear_bore, pressure_angle, helix_angle, tooth_step, flat_step)

    if pinion_teeth % 2 == 0:
        for tri in gear_1:
            rotate([0, 0, pi*(1-clearance)/gear_teeth], tri)

    gear_2 = bevel_gear_assembly(modul, pinion_teeth, delta_pinion, tooth_width, pinion_bore, pressure_angle, -helix_angle, tooth_step, flat_step)

    if together_built:
        for tri in gear_2:
            rotate([0, axis_angle, 0], tri)
            dx = -height_f_pinion*cos(pi/2-axis_angle)
            dz = height_f_gear-height_f_pinion*sin(pi/2-axis_angle)
            translate([dx, 0, dz], tri)
    else:
        for tri in gear_2:
            translate([rkf_pinion*2 + modul + rkf_gear, 0, 0], tri)

    # you can have rotate and translate take list slices, not lists,
    # and add the option to pass ans into bevel_gear
    # so you don't have to do this copy.
    return gear_1 + gear_2

def bevel_herringbone_gear_pair_assembly(modul, gear_teeth, pinion_teeth, axis_angle, tooth_width, gear_bore, pinion_bore, pressure_angle, helix_angle, together_built, tooth_step, flat_step):


    r_gear = modul*gear_teeth/2
    delta_gear = atan(sin(axis_angle)/(pinion_teeth/gear_teeth+cos(axis_angle)))
    delta_pinion = atan(sin(axis_angle)/(gear_teeth/pinion_teeth+cos(axis_angle)))
    rg = r_gear/sin(delta_gear)
    c = modul / 6
    df_pinion = rg*delta_pinion*2 - 2 * (modul + c)
    rf_pinion = df_pinion / 2
    delta_f_pinion = rf_pinion/rg
    rkf_pinion = rg*sin(delta_f_pinion)
    height_f_pinion = rg*cos(delta_f_pinion)

    df_gear = 2*rg*delta_gear - 2 * (modul + c)
    rf_gear = df_gear / 2

    delta_f_gear = rf_gear/rg
    rkf_gear = rg*sin(delta_f_gear)
    height_f_gear = rg*cos(delta_f_gear)

    gear_1 = bevel_herringbone_gear_assembly(modul, gear_teeth, delta_gear, tooth_width, gear_bore, pressure_angle, helix_angle, tooth_step, flat_step)
    
    gear_2 = bevel_herringbone_gear_assembly(modul, pinion_teeth, delta_pinion, tooth_width, pinion_bore, pressure_angle, -helix_angle, tooth_step, flat_step)


    if pinion_teeth % 2 == 0:
        for tri in gear_1:
            rotate([0,0,pi*(1-clearance)/gear_teeth], tri)

    if together_built:
        for tri in gear_2:
            rotate([0, axis_angle, 0], tri)
            dx = -height_f_pinion*cos(pi/2-axis_angle)
            dz = height_f_gear-height_f_pinion*sin(pi/2-axis_angle)
            translate([dx, 0, dz], tri)
    else:
        for tri in gear_2:
            translate([rkf_pinion*2 + modul + rkf_gear, 0, 0], tri)

    return gear_1 + gear_2

