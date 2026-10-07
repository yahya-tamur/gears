from .mesh import Mesh
from .gears import bevel_herringbone_gear_data, herringbone_ring_gear_data
from .math import rad


def bevel_gear(
    modul,
    tooth_number,
    partial_cone_angle,
    tooth_width,
    bore=0,
    pressure_angle=20,
    helix_angle=10,
    tooth_steps=16,
    flat_steps=5,
    helix_steps=10,
    bore_steps=20,
    filename="a.stl",
):

    tooth_west, tooth_east = bevel_herringbone_gear_data(
        modul=modul,
        tooth_number=tooth_number,
        partial_cone_angle=rad(partial_cone_angle),
        tooth_width=tooth_width,
        pressure_angle=rad(pressure_angle),
        helix_angle=rad(helix_angle),
        tooth_steps=tooth_steps,
        helix_steps=helix_steps,
    )

    mesh = Mesh()
    z = mesh.add_gear(
        tooth_west,
        tooth_east,
        tooth_number=tooth_number,
        flat_steps=flat_steps,
        bore=bore,
        bore_steps=bore_steps,
        ring_gear=False,
    )

    mesh.translate((0, 0, -z))

    mesh.save_stl(filename)


def ring_gear(
    modul,
    tooth_number,
    width,
    rim_width=5,
    pressure_angle=20,
    helix_angle=10,
    shortening_factor=1,
    tooth_steps=16,
    flat_steps=5,
    helix_steps=10,
    cylinder_steps=20,
    filename="a.stl",
):
    tooth_west, tooth_east, radius = herringbone_ring_gear_data(
        modul=modul,
        tooth_number=tooth_number,
        width=width,
        rim_width=rim_width,
        pressure_angle=rad(pressure_angle),
        helix_angle=rad(helix_angle),
        shortening_factor=shortening_factor,
        helix_steps=helix_steps,
        tooth_steps=tooth_steps,
    )

    mesh = Mesh()
    z = mesh.add_gear(
        tooth_west,
        tooth_east,
        tooth_number=tooth_number,
        flat_steps=flat_steps,
        bore=radius,
        bore_steps=cylinder_steps,
        return_z_offset=True,
        ring_gear=True,
    )

    print(z)

    mesh.save_stl(filename)


# def bevel_gear(
#    modul,
#    tooth_number,
#    partial_cone_angle,
#    tooth_width,
#    bore=0,
#    pressure_angle=20,
#    helix_angle=0,
#    tooth_step=16,
# ):
#    return bevel_gear_assembly(
#        modul,
#        tooth_number,
#        rad(partial_cone_angle),
#        tooth_width,
#        bore,
#        rad(pressure_angle),
#        rad(helix_angle),
#        tooth_step,
#    )

# def bevel_gear_pair(
#    modul,
#    gear_teeth,
#    pinion_teeth,
#    tooth_width,
#    axis_angle=90,
#    gear_bore=0,
#    pinion_bore=0,
#    pressure_angle=20,
#    helix_angle=0,
#    together_built=True,
#    tooth_step=16,
# ):
#    return bevel_gear_pair_assembly(
#        modul,
#        gear_teeth,
#        pinion_teeth,
#        rad(axis_angle),
#        tooth_width,
#        gear_bore,
#        pinion_bore,
#        rad(pressure_angle),
#        rad(helix_angle),
#        together_built,
#        tooth_step,
#    )


# def bevel_herringbone_gear_pair(
#    modul,
#    gear_teeth,
#    pinion_teeth,
#    tooth_width,
#    axis_angle=90,
#    gear_bore=0,
#    pinion_bore=0,
#    pressure_angle=20,
#    helix_angle=10,
#    together_built=True,
#    tooth_step=16,
# ):
#    return bevel_herringbone_gear_pair_assembly(
#        modul,
#        gear_teeth,
#        pinion_teeth,
#        rad(axis_angle),
#        tooth_width,
#        gear_bore,
#        pinion_bore,
#        rad(pressure_angle),
#        rad(helix_angle),
#        together_built,
#        tooth_step,
#    )


# def planetary_gear(
#    modul,
#    sun_teeth,
#    planet_teeth,
#    width,
#    number_planets=0,
#    rim_width=5,
#    sun_bore=2,
#    planet_bore=1,
#    pressure_angle=20,
#    helix_angle=10,
#    together_built=True,
#    tooth_step=16,
#    ring_shortening_factor=1,
# ):
#    return planetary_gear_assembly(
#        modul,
#        sun_teeth,
#        planet_teeth,
#        number_planets,
#        width,
#        rim_width,
#        sun_bore,
#        planet_bore,
#        rad(pressure_angle),
#        rad(helix_angle),
#        together_built,
#        tooth_step,
#        ring_shortening_factor,
#    )


# might be useful from cli


# def disp(mesh, viewer="fstl", filename="a.stl"):
#    make_stl(mesh, filename)
#    os.system(f"{viewer} {filename}")
