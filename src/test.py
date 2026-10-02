
from .gears import ring_gear_assembly
from .lib import disp

mesh = ring_gear_assembly(modul=2, tooth_number=20, width=10, rim_width=10, pressure_angle=0.349066, helix_angle=0.349066, shortening_factor=0.6, tooth_step=16)
disp(mesh)
