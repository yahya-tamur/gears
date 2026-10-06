from math import sin, cos, tan, atan, asin, acos, pi, floor
from .math import (
    center,
    rotate,
    translate,
    sphere_ev,
    sph_to_cart,
    polar_ev,
    pol_to_cart,
)

# later maybe??
# from mpmath import mp
# mp.dps = 500
# sin = mp.sin
# cos = mp.cos
# tan = mp.tan
# atan = mp.atan
# asin = mp.asin
# acos = mp.acos
# pi = mp.pi
# sqrt = mp.sqrt

# all angles in this file are in radians.

clearance = 0.05


def bevel_gear_data(
    modul,
    tooth_number,
    partial_cone_angle,
    tooth_width,
    pressure_angle,
    helix_angle,
    tooth_step,
    helix_step,
):

    d_outside = modul * tooth_number
    r_outside = d_outside / 2
    rg_outside = r_outside / sin(partial_cone_angle)
    rg_inside = rg_outside - tooth_width
    # r_inside = r_outside * rg_inside / rg_outside
    alpha_spur = atan(tan(pressure_angle) / cos(helix_angle))
    da_outside = d_outside + (modul * 2) * cos(partial_cone_angle)
    ra_outside = da_outside / 2
    c = modul / 6
    df_outside = d_outside - (modul + c) * 2 * cos(partial_cone_angle)
    rf_outside = df_outside / 2
    # rkf = rg_outside*sin(delta_f)
    delta_f = asin(rf_outside / rg_outside)
    delta_a = asin(ra_outside / rg_outside)
    delta_b = asin(cos(alpha_spur) * sin(partial_cone_angle))

    #height_f = rg_outside * cos(delta_f)

    # height_k = (rg_outside - tooth_width) / cos(partial_cone_angle)
    # rk = (rg_outside - tooth_width) / sin(partial_cone_angle)
    # rfk = rk*height_k*tan(delta_f)/(rk+height_k*tan(delta_f))
    # height_fk = rk*height_k/(height_k*tan(delta_f)+rk)

    phi_r = sphere_ev(delta_b, partial_cone_angle)

    gamma_g = 2 * atan(tooth_width * tan(helix_angle) / (2 * rg_outside - tooth_width))
    gamma = 2 * asin(rg_outside / r_outside * sin(gamma_g / 2))

    tau = 2 * pi / tooth_number

    mirrpoint = (pi * (1 - clearance)) / tooth_number + 2 * phi_r

    teeth_west = [[] for _ in range(helix_step + 1)]
    teeth_east = [[] for _ in range(helix_step + 1)]
    ## $0 -> north, $-1: south

    start = delta_f
    step = (delta_a - delta_f) / tooth_step  # check if this is good
    if delta_b > delta_f:
        flankpoint_under = 1 * mirrpoint

        for j in range(helix_step + 1):
            rrr = (rg_outside * (helix_step - j) + rg_inside * j) / helix_step
            ggg = j * gamma / helix_step
            teeth_east[j].append(sph_to_cart((rrr, delta_f, flankpoint_under + ggg)))
            teeth_west[j].append(
                sph_to_cart((rrr, delta_f, mirrpoint - flankpoint_under + ggg))
            )

        #        tooth_ne.append(sph_to_cart((rg_outside, delta_f, flankpoint_under)))
        #        tooth_se.append(sph_to_cart((rg_inside, delta_f, flankpoint_under + gamma)))
        #        tooth_sw.append(
        #            sph_to_cart((rg_inside, delta_f, mirrpoint - flankpoint_under + gamma))
        #        )
        #        tooth_nw.append(
        #            sph_to_cart((rg_outside, delta_f, mirrpoint - flankpoint_under))
        #        )

        start = delta_b
        step = (delta_a - delta_b) / tooth_step

    for i in range(tooth_step + 1):
        delta = start + i * step
        flankpoint_under = sphere_ev(delta_b, delta)

        for j in range(helix_step + 1):
            rrr = (rg_outside * (helix_step - j) + rg_inside * j) / helix_step
            ggg = j * gamma / helix_step
            teeth_west[j].append(sph_to_cart((rrr, delta, flankpoint_under + ggg)))
            teeth_east[j].append(
                sph_to_cart((rrr, delta, mirrpoint - flankpoint_under + ggg))
            )

    #        tooth_nw.append(sph_to_cart((rg_outside, delta, flankpoint_under)))
    #        tooth_sw.append(sph_to_cart((rg_inside, delta, flankpoint_under + gamma)))
    #        tooth_se.append(
    #            sph_to_cart((rg_inside, delta, mirrpoint - flankpoint_under + gamma))
    #        )
    #        tooth_ne.append(sph_to_cart((rg_outside, delta, mirrpoint - flankpoint_under)))

    for line in teeth_west:
        #    for pt_list in (tooth_nw, tooth_ne, tooth_sw, tooth_se):
        rotate([0, pi, 0], line)
        rotate([0, 0, phi_r + pi / 2 * (1 - clearance) / tooth_number], line)

    for line in teeth_east:
        #    for pt_list in (tooth_nw, tooth_ne, tooth_sw, tooth_se):
        rotate([0, pi, 0], line)
        rotate([0, 0, phi_r + pi / 2 * (1 - clearance) / tooth_number], line)

    return (teeth_west, teeth_east, -tau)


def bevel_gear_data_old(
    modul,
    tooth_number,
    partial_cone_angle,
    tooth_width,
    pressure_angle,
    helix_angle,
    tooth_step,
):

    d_outside = modul * tooth_number
    r_outside = d_outside / 2
    rg_outside = r_outside / sin(partial_cone_angle)
    rg_inside = rg_outside - tooth_width
    # r_inside = r_outside * rg_inside / rg_outside
    alpha_spur = atan(tan(pressure_angle) / cos(helix_angle))
    da_outside = d_outside + (modul * (2.2 if modul < 1 else 2)) * cos(
        partial_cone_angle
    )
    ra_outside = da_outside / 2
    c = modul / 6
    df_outside = d_outside - (modul + c) * 2 * cos(partial_cone_angle)
    rf_outside = df_outside / 2
    # rkf = rg_outside*sin(delta_f)
    delta_f = asin(rf_outside / rg_outside)
    delta_a = asin(ra_outside / rg_outside)
    delta_b = asin(cos(alpha_spur) * sin(partial_cone_angle))

    height_f = rg_outside * cos(delta_f)

    # height_k = (rg_outside - tooth_width) / cos(partial_cone_angle)
    # rk = (rg_outside - tooth_width) / sin(partial_cone_angle)
    # rfk = rk*height_k*tan(delta_f)/(rk+height_k*tan(delta_f))
    # height_fk = rk*height_k/(height_k*tan(delta_f)+rk)

    phi_r = sphere_ev(delta_b, partial_cone_angle)

    gamma_g = 2 * atan(tooth_width * tan(helix_angle) / (2 * rg_outside - tooth_width))
    gamma = 2 * asin(rg_outside / r_outside * sin(gamma_g / 2))

    tau = 2 * pi / tooth_number

    mirrpoint = (pi * (1 - clearance)) / tooth_number + 2 * phi_r

    tooth_nw = []
    tooth_ne = []
    tooth_sw = []
    tooth_se = []

    start = delta_f
    step = (delta_a - delta_f) / tooth_step  # check if this is good
    if delta_b > delta_f:
        flankpoint_under = 1 * mirrpoint

        tooth_ne.append(sph_to_cart((rg_outside, delta_f, flankpoint_under)))
        tooth_se.append(sph_to_cart((rg_inside, delta_f, flankpoint_under + gamma)))
        tooth_sw.append(
            sph_to_cart((rg_inside, delta_f, mirrpoint - flankpoint_under + gamma))
        )
        tooth_nw.append(
            sph_to_cart((rg_outside, delta_f, mirrpoint - flankpoint_under))
        )

        start = delta_b
        step = (delta_a - delta_b) / tooth_step

    # for delta in range(start, delta_a, step):
    # what was I thinking here: ????
    # for i in range(tooth_step + 1):
    #     delta = (start * (tooth_step - i) + delta_a * i) / tooth_step
    for i in range(tooth_step + 1):
        delta = start + i * step
        flankpoint_under = sphere_ev(delta_b, delta)

        tooth_nw.append(sph_to_cart((rg_outside, delta, flankpoint_under)))
        tooth_sw.append(sph_to_cart((rg_inside, delta, flankpoint_under + gamma)))
        tooth_se.append(
            sph_to_cart((rg_inside, delta, mirrpoint - flankpoint_under + gamma))
        )
        tooth_ne.append(sph_to_cart((rg_outside, delta, mirrpoint - flankpoint_under)))

    for pt_list in (tooth_nw, tooth_ne, tooth_sw, tooth_se):
        rotate([0, pi, 0], pt_list)
        translate([0, 0, height_f], pt_list)
        rotate([0, 0, phi_r + pi / 2 * (1 - clearance) / tooth_number], pt_list)

    return (tooth_nw, tooth_ne, tooth_sw, tooth_se, -tau)


def flat_herringbone_gear_data(
    modul, tooth_number, tooth_width, pressure_angle, helix_angle, tooth_step
):

    # since the original calculations are in spherical coordinates,
    # I couldn't easily accomodate for this case.
    # This value isn't too small for numerical accuracy.
    partial_cone_angle = 0.000001

    tooth_width = tooth_width / 2
    d_outside = modul * tooth_number
    r_outside = d_outside / 2
    rg_outside = r_outside / sin(partial_cone_angle)

    gamma_g = 2 * atan(tooth_width * tan(helix_angle) / (2 * rg_outside - tooth_width))
    gamma = 2 * asin(rg_outside / r_outside * sin(gamma_g / 2))

    tooth_aw, tooth_ae, _, _, tau = bevel_gear_data(
        modul,
        tooth_number,
        partial_cone_angle,
        tooth_width,
        pressure_angle,
        helix_angle,
        tooth_step,
    )

    tooth_bw = [0 for _ in range(len(tooth_aw))]
    tooth_be = [0 for _ in range(len(tooth_ae))]
    tooth_cw = [0 for _ in range(len(tooth_aw))]
    tooth_ce = [0 for _ in range(len(tooth_ae))]

    for i in range(len(tooth_aw)):
        a, b, _c = tooth_aw[i]
        tooth_aw[i] = (a, b, 0)
        tooth_bw[i] = (a, b, tooth_width)
        tooth_cw[i] = (a, b, 2 * tooth_width)

    for i in range(len(tooth_ae)):
        a, b, _c = tooth_ae[i]
        tooth_ae[i] = (a, b, 0)
        tooth_be[i] = (a, b, tooth_width)
        tooth_ce[i] = (a, b, 2 * tooth_width)

    rotate((0, 0, -gamma), tooth_aw)
    rotate((0, 0, -gamma), tooth_ae)
    rotate((0, 0, -gamma), tooth_cw)
    rotate((0, 0, -gamma), tooth_ce)

    return tooth_aw, tooth_ae, tooth_bw, tooth_be, tooth_cw, tooth_ce, tau


def bevel_herringbone_gear_data(
    modul,
    tooth_number,
    partial_cone_angle,
    tooth_width,
    pressure_angle,
    helix_angle,
    tooth_step,
    helix_step,
):

    tooth_width = tooth_width / 2
    d_outside = modul * tooth_number
    r_outside = d_outside / 2
    rg_outside = r_outside / sin(partial_cone_angle)
    c = modul / 6
    df_outside = d_outside - (modul + c) * 2 * cos(partial_cone_angle)
    #rf_outside = df_outside / 2
    #delta_f = asin(rf_outside / rg_outside)
    #height_f = rg_outside * cos(delta_f)

    gamma_g = 2 * atan(tooth_width * tan(helix_angle) / (2 * rg_outside - tooth_width))
    gamma = 2 * asin(rg_outside / r_outside * sin(gamma_g / 2))

    #height_k = (rg_outside - tooth_width) / cos(partial_cone_angle)
    #rk = (rg_outside - tooth_width) / sin(partial_cone_angle)
    # rfk = rk * height_k * tan(delta_f) / (rk + height_k * tan(delta_f))
    #height_fk = rk * height_k / (height_k * tan(delta_f) + rk)

    modul_inside = modul * (1 - tooth_width / rg_outside)

    # I think the -1 degree correction added here is a band-aid to the real issue
    # of the gears not meshing.
    # For now I address this by averaging the two different middle gears instead.
    # Later I would like to try to address this differently:
    # I think it's a math error, not due to floating points -- the mean square
    # difference between the two gears didn't change at all when I tried using
    # more precise floating points with mpmath.

    lower_cone_angle = partial_cone_angle  # - pi/180

    tooth_top_west, tooth_top_east, tau = bevel_gear_data(
        modul,
        tooth_number,
        lower_cone_angle,
        tooth_width,
        pressure_angle,
        helix_angle,
        tooth_step,
        helix_step,
    )

    tooth_bottom_west, tooth_bottom_east, tau = bevel_gear_data(
        modul_inside,
        tooth_number,
        partial_cone_angle,
        tooth_width,
        pressure_angle,
        -helix_angle,
        tooth_step,
        helix_step,
    )

    for pt_list in tooth_bottom_west:
        rotate([0, 0, -gamma], pt_list)

    for pt_list in tooth_bottom_east:
        rotate([0, 0, -gamma], pt_list)
    #        translate([0, 0, height_f - height_fk], pt_list)

    #tooth_bw, tooth_be = [], []
    v_west = tooth_top_west.pop()
    v_east = tooth_top_east.pop()

    tooth_top_west.append([])
    tooth_top_east.append([])
    for i in range(len(v_west)):
        tooth_top_west[-1].append(center((v_west[i], tooth_bottom_west[0][i])))
        tooth_top_east[-1].append(center((v_east[i], tooth_bottom_east[0][i])))

    tooth_top_west += tooth_bottom_west[1:]
    tooth_top_east += tooth_bottom_east[1:]

    return tooth_top_west, tooth_top_east


#    return tooth_aw, tooth_ae, tooth_bw, tooth_be, tooth_cw, tooth_ce, tau


#def bevel_herringbone_gear_assembly(
#    modul,
#    tooth_number,
#    partial_cone_angle,
#    tooth_width,
#    bore,
#    pressure_angle,
#    helix_angle,
#    tooth_step,
#):
#    if partial_cone_angle == 0:
#        tooth_aw, tooth_ae, tooth_bw, tooth_be, tooth_cw, tooth_ce, tau = (
#            flat_herringbone_gear_data(
#                modul,
#                tooth_number,
#                tooth_width,
#                pressure_angle,
#                helix_angle,
#                tooth_step,
#            )
#        )
#    else:
#        tooth_aw, tooth_ae, tooth_bw, tooth_be, tooth_cw, tooth_ce, tau = (
#            bevel_herringbone_gear_data(
#                modul,
#                tooth_number,
#                partial_cone_angle,
#                tooth_width,
#                pressure_angle,
#                helix_angle,
#                tooth_step,
#            )
#        )

#    tooth_a = tooth_aw + tooth_ae[::-1]
#    tooth_b = tooth_bw + tooth_be[::-1]
#    tooth_c = tooth_cw + tooth_ce[::-1]
#
#    a_face, b_face, c_face = [], [], []#
#
#    ans = []
##
    # create top and bottom teeth faces, teeth open prisms

#    i = 0
#    while True:
#        ans += triangulate_polyhedron(tooth_a)
#        ans += triangulate_polyhedron(tooth_c, reverse=True)
#        ans += triangulate_prism(tooth_b, tooth_a, closed=False)
#        ans += triangulate_prism(tooth_c, tooth_b, closed=False)#

#        if len(a_face) > 0:
#            a_line = [a_face[-1], tooth_a[0]]
#            b_line = [b_face[-1], tooth_b[0]]
#            c_line = [c_face[-1], tooth_c[0]]
#            ans += triangulate_prism(b_line, a_line, closed=False)
#            ans += triangulate_prism(c_line, b_line, closed=False)
#            a_face.pop()
#            b_face.pop()
#            c_face.pop()
#            a_face += a_line
#            b_face += b_line
#            c_face += c_line
#            a_face.append(tooth_a[-1])
#            b_face.append(tooth_b[-1])
#            c_face.append(tooth_c[-1])
#        else:
#            a_face = [tooth_a[0], tooth_a[-1]]
#            b_face = [tooth_b[0], tooth_b[-1]]
#            c_face = [tooth_c[0], tooth_c[-1]]

#        i += 1
#        if i == tooth_number:
#            break
#        rotate((0, 0, tau), tooth_a)
#        rotate((0, 0, tau), tooth_b)
#        rotate((0, 0, tau), tooth_c)

#    a_line = [a_face[-1], a_face[0]]
#    b_line = [b_face[-1], b_face[0]]
#    c_line = [c_face[-1], c_face[0]]
#    ans += triangulate_prism(b_line, a_line, closed=False)
#    ans += triangulate_prism(c_line, b_line, closed=False)

#    a_line.pop()
#    c_line.pop()
#    a_face.pop()
#    c_face.pop()
#    a_face += a_line
#    c_face += c_line

#    if bore == 0:
#        ans += triangulate_polyhedron(a_face)
#        ans += triangulate_polyhedron(c_face, reverse=True)

#    else:
#        top_bore, bottom_bore = [], []
#        for k in range(len(a_face)):
#            top_bore.append(
#                (
#                    -bore / 2 * cos(2 * pi * k / len(a_face)),
#                    bore / 2 * sin(2 * pi * k / len(a_face)),
#                    a_face[k][2],
#                )
#            )
#            bottom_bore.append(
#                (
#                    -bore / 2 * cos(2 * pi * k / len(a_face)),
#                    bore / 2 * sin(2 * pi * k / len(a_face)),
#                    c_face[k][2],
#                )
#            )

#        ans += triangulate_prism(a_face, top_bore, closed=True)
#        ans += triangulate_prism(top_bore, bottom_bore, closed=True)
#        ans += triangulate_prism(bottom_bore, c_face, closed=True)

#    return ans


#def bevel_gear_pair_assembly(
#    modul,
#    gear_teeth,
#    pinion_teeth,
#    axis_angle,
#    tooth_width,
#    gear_bore,
#    pinion_bore,
#    pressure_angle,
#    helix_angle,
#    together_built,
#    tooth_step,
#):

#    r_gear = modul * gear_teeth / 2
#    delta_gear = atan(sin(axis_angle) / (pinion_teeth / gear_teeth + cos(axis_angle)))
#    delta_pinion = atan(sin(axis_angle) / (gear_teeth / pinion_teeth + cos(axis_angle)))
#    rg = r_gear / sin(delta_gear)
#    c = modul / 6
#    df_pinion = 2 * rg * delta_pinion - 2 * (modul + c)
#    rf_pinion = df_pinion / 2
#    delta_f_pinion = rf_pinion / (pi * rg) * pi
#    rkf_pinion = rg * sin(delta_f_pinion)
#    height_f_pinion = rg * cos(delta_f_pinion)

#    df_gear = 2 * rg * delta_gear - 2 * (modul + c)
#    rf_gear = df_gear / 2
#    delta_f_gear = rf_gear / rg
#    rkf_gear = rg * sin(delta_f_gear)
#    height_f_gear = rg * cos(delta_f_gear)

#    gear_1 = bevel_gear_assembly(
#        modul,
#        gear_teeth,
#        delta_gear,
#        tooth_width,
#        gear_bore,
#        pressure_angle,
#        helix_angle,
#        tooth_step,
#    )

#    if pinion_teeth % 2 == 0:
#        for tri in gear_1:
#            rotate([0, 0, pi * (1 - clearance) / gear_teeth], tri)

#    gear_2 = bevel_gear_assembly(
#        modul,
#        pinion_teeth,
#        delta_pinion,
#        tooth_width,
#        pinion_bore,
#        pressure_angle,
#        -helix_angle,
#        tooth_step,
#    )

#    if together_built:
#        for tri in gear_2:
#            rotate([0, axis_angle, 0], tri)
#            dx = -height_f_pinion * cos(pi / 2 - axis_angle)
#            dz = height_f_gear - height_f_pinion * sin(pi / 2 - axis_angle)
#            translate([dx, 0, dz], tri)
#    else:
#        for tri in gear_2:
#            translate([rkf_pinion * 2 + modul + rkf_gear, 0, 0], tri)#

    # you can have rotate and translate take list slices, not lists,
    # and add the option to pass ans into bevel_gear
    # so you don't have to do this copy.
#    return gear_1 + gear_2


def bevel_herringbone_gear_pair_assembly(
    modul,
    gear_teeth,
    pinion_teeth,
    axis_angle,
    tooth_width,
    gear_bore,
    pinion_bore,
    pressure_angle,
    helix_angle,
    together_built,
    tooth_step,
):

    if axis_angle == 0:
        gear_1 = bevel_herringbone_gear_assembly(
            modul,
            gear_teeth,
            0,
            tooth_width,
            gear_bore,
            pressure_angle,
            helix_angle,
            tooth_step,
        )
        gear_2 = bevel_herringbone_gear_assembly(
            modul,
            pinion_teeth,
            0,
            tooth_width,
            pinion_bore,
            pressure_angle,
            -helix_angle,
            tooth_step,
        )

        if pinion_teeth % 2 == 0:
            for tri in gear_1:
                rotate([0, 0, pi * (1 - clearance) / gear_teeth], tri)

        if together_built:
            for tri in gear_2:
                dx = -(pinion_teeth + gear_teeth) * modul / 2
                translate([dx, 0, 0], tri)
        else:
            for tri in gear_2:
                translate([modul * (pinion_teeth + gear_teeth / 2 - 2.5), 0, 0], tri)

    else:
        r_gear = modul * gear_teeth / 2
        delta_gear = atan(
            sin(axis_angle) / (pinion_teeth / gear_teeth + cos(axis_angle))
        )
        delta_pinion = atan(
            sin(axis_angle) / (gear_teeth / pinion_teeth + cos(axis_angle))
        )
        rg = r_gear / sin(delta_gear)
        c = modul / 6
        df_pinion = rg * delta_pinion * 2 - 2 * (modul + c)
        rf_pinion = df_pinion / 2
        delta_f_pinion = rf_pinion / rg
        rkf_pinion = rg * sin(delta_f_pinion)
        height_f_pinion = rg * cos(delta_f_pinion)

        df_gear = 2 * rg * delta_gear - 2 * (modul + c)
        rf_gear = df_gear / 2

        delta_f_gear = rf_gear / rg
        rkf_gear = rg * sin(delta_f_gear)
        height_f_gear = rg * cos(delta_f_gear)

        gear_1 = bevel_herringbone_gear_assembly(
            modul,
            gear_teeth,
            delta_gear,
            tooth_width,
            gear_bore,
            pressure_angle,
            helix_angle,
            tooth_step,
        )

        gear_2 = bevel_herringbone_gear_assembly(
            modul,
            pinion_teeth,
            delta_pinion,
            tooth_width,
            pinion_bore,
            pressure_angle,
            -helix_angle,
            tooth_step,
        )

        if pinion_teeth % 2 == 0:
            for tri in gear_1:
                rotate([0, 0, pi * (1 - clearance) / gear_teeth], tri)

        if together_built:
            for tri in gear_2:
                rotate([0, axis_angle, 0], tri)
                dx = -height_f_pinion * cos(pi / 2 - axis_angle)
                dz = height_f_gear - height_f_pinion * sin(pi / 2 - axis_angle)
                translate([dx, 0, dz], tri)
        else:
            for tri in gear_2:
                translate([rkf_pinion * 2 + modul + rkf_gear, 0, 0], tri)

    return gear_1 + gear_2


def herringbone_ring_gear_data(
    modul,
    tooth_number,
    width,
    rim_width,
    pressure_angle,
    helix_angle,
    shortening_factor,
    tooth_step,
    helix_step,
):
    width = width / 2  # !!!
    ha = shortening_factor
    d = modul * tooth_number
    r = d / 2
    alpha_spur = atan(tan(pressure_angle) / cos(helix_angle))
    db = d * cos(alpha_spur)
    rb = db / 2
    c = modul / 6
    da = d + (modul + c) * 2.2 if (modul < 1) else d + (modul + c) * 2
    ra = da / 2
    df = d - 2 * modul * ha
    rf = df / 2
    rho_ra = acos(rb / ra)

    rho_r = acos(rb / r)

    phi_r = tan(rho_r) - rho_r
    gamma = width / (r * tan(pi / 2 - helix_angle))
    step = rho_ra / tooth_step
    tau = 2 * pi / tooth_number

    tooth_width = (pi * (1 + clearance)) / tooth_number + 2 * phi_r

    offset = -phi_r - (pi / 2) * (1 + clearance) / tooth_number

    teeth_west = [[] for _ in range((2 * helix_step + 1))]
    teeth_east = [[] for _ in range((2 * helix_step + 1))]

    for i in range(tooth_step + 1):
        rho = i * step
        a, b = polar_ev(rb, rho)
        if a < rf:
            continue
        for j in range(helix_step):
            teeth_west[j].append(
                pol_to_cart(
                    a,
                    b + tau - gamma * (helix_step - j) / helix_step + offset,
                    z=j * width / helix_step,
                )
            )
            teeth_east[j].append(
                pol_to_cart(
                    a,
                    tooth_width - b - gamma * (helix_step - j) / helix_step + offset,
                    z=j * width / helix_step,
                )
            )

        for j in range(helix_step, 2 * helix_step + 1):
            teeth_west[j].append(
                pol_to_cart(
                    a,
                    b + tau - gamma * (j - helix_step) / helix_step + offset,
                    z=j * width / helix_step,
                )
            )
            teeth_east[j].append(
                pol_to_cart(
                    a,
                    tooth_width - b - gamma * (j - helix_step) / helix_step + offset,
                    z=j * width / helix_step,
                )
            )

    teeth = [teeth_east[i][::-1] + teeth_west[i] for i in range(2 * helix_step + 1)]

    outer_top = [
        pol_to_cart(ra + rim_width, -gamma + offset, z=2 * width),
        pol_to_cart(ra + rim_width, -gamma + tau / 2 + offset, z=2 * width),
    ]
    outer_bottom = [
        pol_to_cart(ra + rim_width, -gamma + offset, z=0),
        pol_to_cart(ra + rim_width, -gamma + tau / 2 + offset, z=0),
    ]

    return (
        teeth,
        outer_top,
        outer_bottom,
        tau,
    )


def herringbone_ring_gear_assembly(
    modul,
    tooth_number,
    width,
    rim_width,
    pressure_angle,
    helix_angle,
    shortening_factor,
    tooth_step,
    helix_step,
):
    (
        teeth,
        outer_top,
        outer_bottom,
        tau,
    ) = herringbone_ring_gear_data(
        modul,
        tooth_number,
        width,
        rim_width,
        pressure_angle,
        helix_angle,
        shortening_factor,
        tooth_step,
        helix_step,
    )

    mesh = (
        triangulate_polyhedron(
            [outer_top[0], outer_top[1], teeth[-1][-1], teeth[-1][0]]
        )
        + triangulate_polyhedron(teeth[-1], reverse=True)
        + triangulate_polyhedron(teeth[0])
        + triangulate_polyhedron(
            [outer_bottom[0], outer_bottom[1], teeth[0][-1], teeth[0][0]],
            reverse=True,
        )
        + triangulate_prism(outer_bottom, outer_top, closed=False)
    )

    for i in range(len(teeth) - 1):
        mesh += triangulate_prism(teeth[i + 1], teeth[i], closed=False)

    first_line = [line[0] for line in teeth] + [
        outer_top[0],
        outer_bottom[0],
    ]
    last_line = [line[-1] for line in teeth] + [
        outer_top[1],
        outer_bottom[1],
    ]

    for _ in range(tooth_number - 1):
        for line in teeth:
            rotate([0, 0, tau], line)

        rotate([0, 0, tau], outer_bottom)
        rotate([0, 0, tau], outer_top)

        first_line_ = [line[0] for line in teeth] + [
            outer_top[0],
            outer_bottom[0],
        ]
        last_line_ = [line[-1] for line in teeth] + [
            outer_top[1],
            outer_bottom[1],
        ]

        mesh += (
            triangulate_polyhedron(
                [outer_top[0], outer_top[1], teeth[-1][-1], teeth[-1][0]]
            )
            + triangulate_polyhedron(teeth[-1], reverse=True)
            + triangulate_polyhedron(teeth[0])
            + triangulate_prism(last_line, first_line_, closed=True)
            + triangulate_polyhedron(
                [outer_bottom[0], outer_bottom[1], teeth[0][-1], teeth[0][0]],
                reverse=True,
            )
            + triangulate_prism(outer_bottom, outer_top, closed=False)
        )
        for i in range(len(teeth) - 1):
            mesh += triangulate_prism(teeth[i + 1], teeth[i], closed=False)

        last_line = last_line_
    mesh += triangulate_prism(last_line, first_line, closed=True)

    return mesh


def planetary_gear_assembly(
    modul,
    sun_teeth,
    planet_teeth,
    number_planets,
    width,
    rim_width,
    sun_bore,
    planet_bore,
    pressure_angle,
    helix_angle,
    together_built,
    tooth_step,
    ring_shortening_factor,
):

    d_planet = modul * planet_teeth
    center_distance = modul * (sun_teeth + planet_teeth) / 2
    ring_teeth = sun_teeth + 2 * planet_teeth

    if number_planets == 0:
        max_planets = floor(
            pi / asin(modul * (planet_teeth) / (modul * (sun_teeth + planet_teeth)))
        )
        number_planets = [
            n
            for n in range(2, max_planets + 1)
            if (((ring_teeth + sun_teeth) % n) == 0)
        ][-1]

    sun_gear = bevel_herringbone_gear_assembly(
        modul, sun_teeth, 0, width, sun_bore, pressure_angle, -helix_angle, tooth_step
    )

    if planet_teeth % 2 == 0:
        for tri in sun_gear:
            rotate([0, 0, pi * (1 - clearance) / sun_teeth], tri)

    mesh = sun_gear

    planet_gear = bevel_herringbone_gear_assembly(
        modul,
        planet_teeth,
        0,
        width,
        planet_bore,
        pressure_angle,
        helix_angle,
        tooth_step,
    )

    for n in range(number_planets):
        new_planet = [tri.copy() for tri in planet_gear]
        for tri in new_planet:
            #    rotate([0, 0, n*2*pi*d_sun/d_planet], tri)
            if together_built:
                translate(
                    pol_to_cart(center_distance, 2 * pi * n / number_planets, z=0),
                    tri,
                )
            else:
                planet_distance = ring_teeth * modul / 2 + rim_width + d_planet
                translate(
                    [planet_distance, d_planet * (-(number_planets - 1) + 2 * n), 0],
                    tri,
                )
        mesh += new_planet

    ring_gear = herringbone_ring_gear_assembly(
        modul,
        ring_teeth,
        width,
        rim_width,
        pressure_angle,
        helix_angle,
        ring_shortening_factor,
        tooth_step,
    )

    mesh += ring_gear

    return mesh
