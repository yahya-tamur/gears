# for easier imports: converts degrees to radians and provides some default parameters

from .mesh import Mesh
from .gears import (
    add_gear,
    add_gear_pair,
    add_ring_gear,
    add_planetary_gear,
)
from .math import rad


def make_gear(
    modul,
    tooth_number,
    partial_cone_angle,
    tooth_width,
    bore=0,
    pressure_angle=20,
    helix_angle=10,
    tooth_steps=16,
    flat_steps=3,
    helix_steps=5,
    bore_steps=10,
    da_factor=1,
    filename="a.stl",
):
    mesh = Mesh()
    add_gear(
        mesh=mesh,
        modul=modul,
        tooth_number=tooth_number,
        partial_cone_angle=rad(partial_cone_angle),
        tooth_width=tooth_width,
        bore=bore,
        pressure_angle=rad(pressure_angle),
        helix_angle=rad(helix_angle),
        tooth_steps=tooth_steps,
        flat_steps=flat_steps,
        helix_steps=helix_steps,
        bore_steps=bore_steps,
        da_factor=da_factor,
    )
    mesh.save_stl(filename)


def make_gear_pair(
    modul,
    gear_teeth,
    pinion_teeth,
    axis_angle,
    tooth_width,
    bore=0,
    pressure_angle=20,
    helix_angle=10,
    together_built=True,
    tooth_steps=16,
    flat_steps=3,
    helix_steps=5,
    bore_steps=10,
    da_factor=1,
    filename="a.stl",
):
    mesh = Mesh()
    add_gear_pair(
        mesh=mesh,
        modul=modul,
        gear_teeth=gear_teeth,
        pinion_teeth=pinion_teeth,
        axis_angle=rad(axis_angle),
        tooth_width=tooth_width,
        bore=bore,
        pressure_angle=rad(pressure_angle),
        helix_angle=rad(helix_angle),
        together_built=together_built,
        tooth_steps=tooth_steps,
        flat_steps=flat_steps,
        helix_steps=helix_steps,
        bore_steps=bore_steps,
        da_factor=da_factor,
    )
    mesh.save_stl(filename)


def make_ring_gear(
    modul,
    tooth_number,
    width,
    rim_width=5,
    pressure_angle=20,
    helix_angle=10,
    shortening_factor=0.4,
    tooth_steps=16,
    flat_steps=3,
    helix_steps=5,
    bore_steps=20,
    da_factor=1,
    filename="a.stl",
):
    mesh = Mesh()
    add_ring_gear(
        mesh=mesh,
        modul=modul,
        tooth_number=tooth_number,
        width=width,
        rim_width=rim_width,
        pressure_angle=rad(pressure_angle),
        helix_angle=rad(helix_angle),
        shortening_factor=shortening_factor,
        tooth_steps=tooth_steps,
        flat_steps=flat_steps,
        helix_steps=helix_steps,
        bore_steps=20,
        da_factor=da_factor,
    )
    mesh.save_stl(filename)


def make_planetary_gear(
    modul,
    sun_teeth,
    planet_teeth,
    number_planets,
    width,
    rim_width=5,
    bore=0,
    pressure_angle=20,
    helix_angle=15,
    ring_shortening_factor=0.4,
    together_built=True,
    tooth_steps=16,
    flat_steps=3,
    helix_steps=5,
    bore_steps=10,
    da_factor=1,
    filename="a.stl",
):
    mesh = Mesh()
    add_planetary_gear(
        mesh=mesh,
        modul=modul,
        sun_teeth=sun_teeth,
        planet_teeth=planet_teeth,
        number_planets=number_planets,
        width=width,
        rim_width=rim_width,
        bore=bore,
        pressure_angle=rad(pressure_angle),
        helix_angle=rad(helix_angle),
        ring_shortening_factor=ring_shortening_factor,
        together_built=together_built,
        tooth_steps=tooth_steps,
        flat_steps=flat_steps,
        helix_steps=helix_steps,
        bore_steps=bore_steps,
        da_factor=da_factor,
    )
    mesh.save_stl(filename)
