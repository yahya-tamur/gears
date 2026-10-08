from os import system

from .gears import add_planetary_gear
from .mesh import Mesh
from .math import rad

# save_gear_pair(
#     modul=2,
#     gear_teeth=20,
#     pinion_teeth=21,
#     axis_angle=0,
#     tooth_width=10,
#     helix_angle=30
# )
# save_gear(
#     modul=2,
#     tooth_number=26,
#     tooth_width=10,
#     partial_cone_angle=0,
#     helix_angle=15,
#     filename="a1.stl"
# )
# save_ring_gear(
#     modul=2,
#     tooth_number=26,
#     width=10,
#     rim_width=5,
#     pressure_angle=rad(20),
#     helix_angle=rad(15),
#     shortening_factor=0,
#     tooth_steps=16,
#     flat_steps=10,
#     helix_steps=10,
#     bore_steps=30,
#     da_factor=1,
# )
# save_ring_gear(
#     modul=2,
#     tooth_number=20,
#     width=10,
#     rim_width=5,
#     pressure_angle=20,
#     helix_angle=10,
#     shortening_factor=0.6,
#     tooth_steps=16,
#     flat_steps=10,
#     helix_steps=10,
#     bore_steps=30,
#     da_factor=1,
#     filename="a1.stl"
# )
mesh = Mesh()
add_planetary_gear(
    mesh=mesh,
    modul=2,
    sun_teeth=10,
    planet_teeth=10,
    number_planets=4,
    width=10,
    rim_width=4,
    bore=0,
    pressure_angle=rad(20),
    helix_angle=rad(20),
    ring_shortening_factor=0.6,
    together_built=True,
    tooth_steps=16,
    flat_steps=5,
    helix_steps=10,
    bore_steps=20,
    da_factor=1,
)
mesh.save_stl("a.stl")
system("fstl a.stl")
