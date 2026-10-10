from math import sin, cos, tan, atan, asin, acos, pi, sqrt, floor
from .math import (
    center,
    rotate,
    sphere_ev,
    sph_to_cart,
    polar_ev,
    pol_to_cart,
    translate,
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


def spur_gear_data(
    modul,
    tooth_number,
    tooth_width,
    pressure_angle,
    helix_angle,
    tooth_steps,
    helix_steps,
    da_factor,
):
    tooth_width = tooth_width / 2  # !!!!
    d = modul * tooth_number
    r = d / 2
    alpha_spur = atan(tan(pressure_angle) / cos(helix_angle))
    db = d * cos(alpha_spur)
    rb = db / 2
    da = d + 2 * modul * da_factor
    ra = da / 2
    c = modul / 6
    df = d - 2 * (modul + c)
    rf = df / 2

    rho_r = acos(rb / r)

    phi_r = tan(rho_r) - rho_r
    gamma = tooth_width / (r * tan(pi / 2 - helix_angle))

    mirrpoint = pi / tooth_number * (1 - clearance) + 2 * phi_r
    offset = -(phi_r + (pi) / 2 / tooth_number * (1 - clearance)) + gamma

    teeth_west = [[] for _ in range(2 * helix_steps + 1)]
    teeth_east = [[] for _ in range(2 * helix_steps + 1)]

    for j in range(helix_steps):
        r, theta = rf, 0
        ggg = j * gamma / helix_steps
        zzz = j * tooth_width / helix_steps

        teeth_west[j].append(pol_to_cart(r, 0 + theta - ggg + offset, z=zzz))
        teeth_west[2 * helix_steps - j].append(
            pol_to_cart(r, 0 + theta - ggg + offset, z=2 * tooth_width - zzz)
        )
        teeth_east[j].append(pol_to_cart(r, mirrpoint - theta - ggg + offset, z=zzz))
        teeth_east[2 * helix_steps - j].append(
            pol_to_cart(r, mirrpoint - theta - ggg + offset, z=2 * tooth_width - zzz)
        )

    teeth_west[helix_steps].append(pol_to_cart(rf, -gamma + offset, z=tooth_width))
    teeth_east[helix_steps].append(
        pol_to_cart(rf, mirrpoint - gamma + offset, z=tooth_width)
    )

    # this should be better than the older version of sampling points (this way,
    # they're evenly spaced along the curve). But, I haven't applied this to
    # bevel or ring gears yet.
    # 0.5tan^2(t) is the length integral of the curve.
    rho_ra = acos(rb / ra)
    t_max = 0.5 * (tan(rho_ra) ** 2)

    for i in range(tooth_steps + 1):
        t = atan(sqrt(2 * (i / tooth_steps) * t_max))

        r, theta = polar_ev(rb, t)
        if r < rf:
            continue
        for j in range(helix_steps):
            ggg = j * gamma / helix_steps
            zzz = j * tooth_width / helix_steps

            teeth_west[j].append(pol_to_cart(r, 0 + theta - ggg + offset, z=zzz))
            teeth_west[2 * helix_steps - j].append(
                pol_to_cart(r, 0 + theta - ggg + offset, z=2 * tooth_width - zzz)
            )
            teeth_east[j].append(
                pol_to_cart(r, mirrpoint - theta - ggg + offset, z=zzz)
            )
            teeth_east[2 * helix_steps - j].append(
                pol_to_cart(
                    r, mirrpoint - theta - ggg + offset, z=2 * tooth_width - zzz
                )
            )

        teeth_west[helix_steps].append(
            pol_to_cart(r, 0 + theta - gamma + offset, z=tooth_width)
        )
        teeth_east[helix_steps].append(
            pol_to_cart(r, mirrpoint - theta - gamma + offset, z=tooth_width)
        )

    return teeth_east, teeth_west


def bevel_gear_data(
    modul,
    tooth_number,
    partial_cone_angle,
    tooth_width,
    pressure_angle,
    helix_angle,
    tooth_steps,
    helix_steps,
    da_factor,
):

    d_outside = modul * tooth_number
    r_outside = d_outside / 2
    rg_outside = r_outside / sin(partial_cone_angle)
    rg_inside = rg_outside - tooth_width
    alpha_spur = atan(tan(pressure_angle) / cos(helix_angle))
    da_outside = d_outside + da_factor * (modul * 2) * cos(partial_cone_angle)
    ra_outside = da_outside / 2
    c = modul / 6
    df_outside = d_outside - (modul + c) * 2 * cos(partial_cone_angle)
    rf_outside = df_outside / 2
    delta_f = asin(rf_outside / rg_outside)
    delta_a = asin(ra_outside / rg_outside)
    delta_b = asin(cos(alpha_spur) * sin(partial_cone_angle))

    phi_r = sphere_ev(delta_b, partial_cone_angle)

    gamma_g = 2 * atan(tooth_width * tan(helix_angle) / (2 * rg_outside - tooth_width))
    gamma = 2 * asin(rg_outside / r_outside * sin(gamma_g / 2))

    mirrpoint = pi / tooth_number * (1 - clearance) + 2 * phi_r

    teeth_west = [[] for _ in range(helix_steps + 1)]
    teeth_east = [[] for _ in range(helix_steps + 1)]

    start = delta_f
    step = (delta_a - delta_f) / tooth_steps
    if delta_b > delta_f:
        flankpoint_under = 1 * mirrpoint

        for j in range(helix_steps + 1):
            rrr = (rg_outside * (helix_steps - j) + rg_inside * j) / helix_steps
            ggg = j * gamma / helix_steps
            teeth_east[j].append(sph_to_cart((rrr, delta_f, flankpoint_under + ggg)))
            teeth_west[j].append(
                sph_to_cart((rrr, delta_f, mirrpoint - flankpoint_under + ggg))
            )

        start = delta_b
        step = (delta_a - delta_b) / tooth_steps

    for i in range(tooth_steps):
        delta = start + i * step
        flankpoint_under = sphere_ev(delta_b, delta)

        for j in range(helix_steps + 1):
            rrr = (rg_outside * (helix_steps - j) + rg_inside * j) / helix_steps
            ggg = j * gamma / helix_steps
            teeth_west[j].append(sph_to_cart((rrr, delta, flankpoint_under + ggg)))
            teeth_east[j].append(
                sph_to_cart((rrr, delta, mirrpoint - flankpoint_under + ggg))
            )

    for line in teeth_west:
        rotate([0, pi, 0], line)
        rotate([0, 0, phi_r + pi / 2 / tooth_number * (1 - clearance)], line)

    for line in teeth_east:
        rotate([0, pi, 0], line)
        rotate([0, 0, phi_r + pi / 2 / tooth_number * (1 - clearance)], line)

    return (teeth_west, teeth_east)


def bevel_herringbone_gear_data(
    modul,
    tooth_number,
    partial_cone_angle,
    tooth_width,
    pressure_angle,
    helix_angle,
    tooth_steps,
    helix_steps,
    da_factor,
):

    tooth_width = tooth_width / 2
    d_outside = modul * tooth_number
    r_outside = d_outside / 2
    rg_outside = r_outside / sin(partial_cone_angle)

    gamma_g = 2 * atan(tooth_width * tan(helix_angle) / (2 * rg_outside - tooth_width))
    gamma = 2 * asin(rg_outside / r_outside * sin(gamma_g / 2))
    modul_inside = modul * (1 - tooth_width / rg_outside)

    lower_cone_angle = partial_cone_angle  # - pi/180

    tooth_top_west, tooth_top_east = bevel_gear_data(
        modul,
        tooth_number,
        lower_cone_angle,
        tooth_width,
        pressure_angle,
        helix_angle,
        tooth_steps,
        helix_steps,
        da_factor,
    )

    tooth_bottom_west, tooth_bottom_east = bevel_gear_data(
        modul_inside,
        tooth_number,
        partial_cone_angle,
        tooth_width,
        pressure_angle,
        -helix_angle,
        tooth_steps,
        helix_steps,
        da_factor,
    )

    for pt_list in tooth_bottom_west:
        rotate([0, 0, -gamma], pt_list)

    for pt_list in tooth_bottom_east:
        rotate([0, 0, -gamma], pt_list)
    v_west = tooth_top_west.pop()
    v_east = tooth_top_east.pop()

    tooth_top_west.append([])
    tooth_top_east.append([])
    for i in range(len(v_west)):
        tooth_top_west[-1].append(center((v_west[i], tooth_bottom_west[0][i])))
        tooth_top_east[-1].append(center((v_east[i], tooth_bottom_east[0][i])))

    tooth_west = tooth_top_west + tooth_bottom_west[1:]
    tooth_east = tooth_top_east + tooth_bottom_east[1:]

    z_average = (tooth_west[0][0][2] + tooth_east[0][0][2]) / 2

    for line in tooth_west:
        translate([0, 0, -z_average], line)
    for line in tooth_east:
        translate([0, 0, -z_average], line)

    return tooth_west, tooth_east


def add_gear(
    mesh,
    modul,
    tooth_number,
    partial_cone_angle,
    tooth_width,
    bore,
    pressure_angle,
    helix_angle,
    tooth_steps,
    flat_steps,
    helix_steps,
    bore_steps,
    da_factor,
):
    if partial_cone_angle == 0:
        tooth_west, tooth_east = spur_gear_data(
            modul=modul,
            tooth_number=tooth_number,
            tooth_width=tooth_width,
            pressure_angle=pressure_angle,
            helix_angle=helix_angle,
            tooth_steps=tooth_steps,
            helix_steps=helix_steps,
            da_factor=da_factor,
        )
    else:
        tooth_west, tooth_east = bevel_herringbone_gear_data(
            modul=modul,
            tooth_number=tooth_number,
            partial_cone_angle=partial_cone_angle,
            tooth_width=tooth_width,
            pressure_angle=pressure_angle,
            helix_angle=helix_angle,
            tooth_steps=tooth_steps,
            helix_steps=helix_steps,
            da_factor=da_factor,
        )

    mesh.add_gear_from_teeth(
        tooth_west,
        tooth_east,
        tooth_number=tooth_number,
        flat_steps=flat_steps,
        bore=bore,
        bore_steps=bore_steps,
        ring_gear=False,
    )


def add_gear_pair(
    mesh,
    modul,
    gear_teeth,
    pinion_teeth,
    axis_angle,
    tooth_width,
    bore,
    pressure_angle,
    helix_angle,
    together_built,
    tooth_steps,
    flat_steps,
    helix_steps,
    bore_steps,
    da_factor,
):
    start = mesh.current_length()
    if axis_angle == 0:
        delta_pinion = 0
        delta_gear = 0
    else:
        delta_gear = atan(
            sin(axis_angle) / (pinion_teeth / gear_teeth + cos(axis_angle))
        )
        delta_pinion = atan(
            sin(axis_angle) / (gear_teeth / pinion_teeth + cos(axis_angle))
        )

    add_gear(
        mesh=mesh,
        modul=modul,
        tooth_number=gear_teeth,
        partial_cone_angle=delta_gear,
        tooth_width=tooth_width,
        bore=bore,
        pressure_angle=pressure_angle,
        helix_angle=helix_angle,
        tooth_steps=tooth_steps,
        flat_steps=flat_steps,
        helix_steps=helix_steps,
        bore_steps=bore_steps,
        da_factor=da_factor,
    )
    middle = mesh.current_length()
    add_gear(
        mesh=mesh,
        modul=modul,
        tooth_number=pinion_teeth,
        partial_cone_angle=delta_pinion,
        tooth_width=tooth_width,
        bore=bore,
        pressure_angle=pressure_angle,
        helix_angle=-helix_angle,
        tooth_steps=tooth_steps,
        flat_steps=flat_steps,
        helix_steps=helix_steps,
        bore_steps=bore_steps,
        da_factor=da_factor,
    )

    if pinion_teeth % 2 == 0:
        mesh.rotate([0, 0, pi / gear_teeth * (1 - clearance)], start=start, end=middle)

    if axis_angle == 0:
        if together_built:
            dx = (pinion_teeth + gear_teeth) * modul / 2
        else:
            dx = -modul * (pinion_teeth + gear_teeth / 2 - 2.5)

        mesh.translate([dx, 0, 0], start=middle)

    else:
        r_gear = modul * gear_teeth / 2
        rg = r_gear / sin(delta_gear)
        c = modul / 6
        df_pinion = 2 * rg * delta_pinion - 2 * (modul + c)
        rf_pinion = df_pinion / 2
        delta_f_pinion = rf_pinion / (pi * rg) * pi
        rkf_pinion = rg * sin(delta_f_pinion)
        height_f_pinion = rg * cos(delta_f_pinion)

        df_gear = 2 * rg * delta_gear - 2 * (modul + c)
        rf_gear = df_gear / 2
        delta_f_gear = rf_gear / rg
        rkf_gear = rg * sin(delta_f_gear)
        height_f_gear = rg * cos(delta_f_gear)

        if together_built:
            mesh.rotate([0, axis_angle, 0], start=middle)
            dx = -height_f_pinion * cos(pi / 2 - axis_angle)
            dz = height_f_gear - height_f_pinion * sin(pi / 2 - axis_angle)
            mesh.translate([dx, 0, dz], start=middle)
        else:
            mesh.translate([rkf_pinion * 2 + modul + rkf_gear, 0, 0], start=middle)


def herringbone_ring_gear_data(
    modul,
    tooth_number,
    width,
    rim_width,
    pressure_angle,
    helix_angle,
    shortening_factor,
    tooth_steps,
    helix_steps,
    da_factor,
):
    width = width / 2  # !!!
    ha = shortening_factor
    d = modul * tooth_number
    r = d / 2
    alpha_spur = atan(tan(pressure_angle) / cos(helix_angle))
    db = d * cos(alpha_spur)
    rb = db / 2
    c = modul / 6
    da = d + (modul + c) * 2 * da_factor
    ra = da / 2

    # calculated differently from original!!
    # it made more sense to me to have higher shortening factor = more shortening
    rf = rb + (ra - rb) * ha

    rho_r = acos(rb / r)

    phi_r = tan(rho_r) - rho_r
    gamma = width / (r * tan(pi / 2 - helix_angle))

    tau = 2 * pi / tooth_number

    mirrpoint = pi / tooth_number * (1 + clearance) + 2 * phi_r

    offset = -phi_r - pi / 2 / tooth_number * (1 + clearance) + gamma

    teeth_west = [[] for _ in range((2 * helix_steps + 1))]
    teeth_east = [[] for _ in range((2 * helix_steps + 1))]

    rho_ra = acos(rb / ra)
    t_max = 0.5 * (tan(rho_ra) ** 2)

    for i in range(tooth_steps + 1):
        t = atan(sqrt(2 * (i / tooth_steps) * t_max))

        r, theta = polar_ev(rb, t)

        if r < rf:
            continue
        for j in range(helix_steps):
            ggg = j * gamma / helix_steps
            zzz = j * width / helix_steps
            teeth_west[j].append(pol_to_cart(r, tau + theta - ggg + offset, z=zzz))
            teeth_west[2 * helix_steps - j].append(
                pol_to_cart(r, tau + theta - ggg + offset, z=2 * width - zzz)
            )
            teeth_east[j].append(
                pol_to_cart(r, mirrpoint - theta - ggg + offset, z=zzz)
            )
            teeth_east[2 * helix_steps - j].append(
                pol_to_cart(r, mirrpoint - theta - ggg + offset, z=2 * width - zzz)
            )
        teeth_west[helix_steps].append(
            pol_to_cart(r, tau + theta - gamma + offset, z=width)
        )
        teeth_east[helix_steps].append(
            pol_to_cart(r, mirrpoint - theta - gamma + offset, z=width)
        )

    return teeth_west, teeth_east, 2 * (ra + rim_width)


def add_ring_gear(
    mesh,
    modul,
    tooth_number,
    width,
    rim_width,
    pressure_angle,
    helix_angle,
    shortening_factor,
    tooth_steps,
    flat_steps,
    helix_steps,
    bore_steps,
    da_factor,
):
    tooth_west, tooth_east, radius = herringbone_ring_gear_data(
        modul=modul,
        tooth_number=tooth_number,
        width=width,
        rim_width=rim_width,
        pressure_angle=pressure_angle,
        helix_angle=helix_angle,
        shortening_factor=shortening_factor,
        helix_steps=helix_steps,
        tooth_steps=tooth_steps,
        da_factor=da_factor,
    )

    mesh.add_gear_from_teeth(
        tooth_west,
        tooth_east,
        tooth_number=tooth_number,
        flat_steps=flat_steps,
        bore=radius,
        bore_steps=bore_steps,
        ring_gear=True,
    )


def add_planetary_gear(
    mesh,
    modul,
    sun_teeth,
    planet_teeth,
    number_planets,
    width,
    rim_width,
    bore,
    pressure_angle,
    helix_angle,
    ring_shortening_factor,
    together_built,
    tooth_steps,
    flat_steps,
    helix_steps,
    bore_steps,
    da_factor,
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

    start = mesh.current_length()
    add_gear(
        mesh=mesh,
        modul=modul,
        tooth_number=sun_teeth,
        partial_cone_angle=0,
        tooth_width=width,
        bore=bore,
        pressure_angle=pressure_angle,
        helix_angle=-helix_angle,
        tooth_steps=tooth_steps,
        flat_steps=flat_steps,
        helix_steps=helix_steps,
        bore_steps=bore_steps,
        da_factor=da_factor,
    )
    sun_end = mesh.current_length()

    if planet_teeth % 2 == 0:
        mesh.rotate([0, 0, pi / sun_teeth], start=start)

    add_gear(
        mesh=mesh,
        modul=modul,
        tooth_number=planet_teeth,
        partial_cone_angle=0,
        tooth_width=width,
        bore=bore,
        pressure_angle=pressure_angle,
        helix_angle=helix_angle,
        tooth_steps=tooth_steps,
        flat_steps=flat_steps,
        helix_steps=helix_steps,
        bore_steps=bore_steps,
        da_factor=da_factor,
    )
    planet_length = mesh.current_length() - sun_end
    for n in range(number_planets - 1):
        mesh.clone(start=sun_end, end=sun_end + planet_length)

    for n in range(number_planets):
        start = sun_end + n * planet_length
        end = sun_end + (n + 1) * planet_length
        # for tri in new_planet:
        #    rotate([0, 0, n*2*pi*d_sun/d_planet], tri)
        if together_built:
            mesh.translate(
                pol_to_cart(center_distance, 2 * pi * n / number_planets, z=0),
                start=start,
                end=end,
            )
        else:
            planet_distance = ring_teeth * modul / 2 + rim_width + d_planet
            mesh.translate(
                [planet_distance, d_planet * (-(number_planets - 1) + 2 * n), 0],
                start=start,
                end=end,
            )

    add_ring_gear(
        mesh=mesh,
        modul=modul,
        tooth_number=ring_teeth,
        width=width,
        rim_width=rim_width,
        pressure_angle=pressure_angle,
        helix_angle=helix_angle,
        shortening_factor=ring_shortening_factor,
        tooth_steps=tooth_steps,
        flat_steps=flat_steps,
        helix_steps=helix_steps,
        bore_steps=bore_steps,
        da_factor=da_factor,
    )
