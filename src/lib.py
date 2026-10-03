# There's probably a more canonical way to do this in python, but this is just
# meant to provide a cleaner file you can import.

# the main functionality is providing default values and converting degrees to radians.

import os

from .gears import (
    bevel_gear_assembly,
    bevel_herringbone_gear_assembly,
    bevel_gear_pair_assembly,
    bevel_herringbone_gear_pair_assembly,
    herringbone_ring_gear_assembly,
    planetary_gear_assembly,
)

from .math import rad

from .make_stl import make_stl


def bevel_gear(
    modul,
    tooth_number,
    partial_cone_angle,
    tooth_width,
    bore=0,
    pressure_angle=20,
    helix_angle=0,
    tooth_step=16,
):
    return bevel_gear_assembly(
        modul,
        tooth_number,
        rad(partial_cone_angle),
        tooth_width,
        bore,
        rad(pressure_angle),
        rad(helix_angle),
        tooth_step,
    )


def bevel_herringbone_gear(
    modul,
    tooth_number,
    partial_cone_angle,
    tooth_width,
    bore=0,
    pressure_angle=20,
    helix_angle=10,
    tooth_step=16,
):
    return bevel_herringbone_gear_assembly(
        modul,
        tooth_number,
        rad(partial_cone_angle),
        tooth_width,
        bore,
        rad(pressure_angle),
        rad(helix_angle),
        tooth_step,
    )


def bevel_gear_pair(
    modul,
    gear_teeth,
    pinion_teeth,
    tooth_width,
    axis_angle=90,
    gear_bore=0,
    pinion_bore=0,
    pressure_angle=20,
    helix_angle=0,
    together_built=True,
    tooth_step=16,
):
    return bevel_gear_pair_assembly(
        modul,
        gear_teeth,
        pinion_teeth,
        rad(axis_angle),
        tooth_width,
        gear_bore,
        pinion_bore,
        rad(pressure_angle),
        rad(helix_angle),
        together_built,
        tooth_step,
    )


def bevel_herringbone_gear_pair(
    modul,
    gear_teeth,
    pinion_teeth,
    tooth_width,
    axis_angle=90,
    gear_bore=0,
    pinion_bore=0,
    pressure_angle=20,
    helix_angle=10,
    together_built=True,
    tooth_step=16,
):
    return bevel_herringbone_gear_pair_assembly(
        modul,
        gear_teeth,
        pinion_teeth,
        rad(axis_angle),
        tooth_width,
        gear_bore,
        pinion_bore,
        rad(pressure_angle),
        rad(helix_angle),
        together_built,
        tooth_step,
    )


def herringbone_ring_gear(
    modul,
    tooth_number,
    width,
    rim_width=5,
    pressure_angle=20,
    helix_angle=10,
    shortening_factor=1,
    tooth_step=16,
):
    return herringbone_ring_gear_assembly(
        modul,
        tooth_number,
        width,
        rim_width,
        rad(pressure_angle),
        rad(helix_angle),
        shortening_factor,
        tooth_step,
    )


def planetary_gear(
    modul,
    sun_teeth,
    planet_teeth,
    width,
    number_planets=0,
    rim_width=5,
    sun_bore=2,
    planet_bore=1,
    pressure_angle=20,
    helix_angle=10,
    together_built=True,
    tooth_step=16,
    ring_shortening_factor=1,
):
    return planetary_gear_assembly(
        modul,
        sun_teeth,
        planet_teeth,
        number_planets,
        width,
        rim_width,
        sun_bore,
        planet_bore,
        rad(pressure_angle),
        rad(helix_angle),
        together_built,
        tooth_step,
        ring_shortening_factor,
    )


# might be useful from cli


def disp(mesh, viewer="fstl", filename="a.stl"):
    make_stl(mesh, filename)
    os.system(f"{viewer} {filename}")
