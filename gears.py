# There's probably a more canonical way to do this in python, but this is just
# meant to provide a cleaner file you can import.

# the main functionality is providing default values and converting degrees to radians.

from gears_internal import (
    rad,
    bevel_gear_assembly, 
    bevel_herringbone_gear_assembly, 
    bevel_gear_pair_assembly, 
    bevel_herringbone_gear_pair_assembly
)

def bevel_gear(modul, \
               tooth_number, \
               partial_cone_angle, \
               tooth_width, \
               bore = 0, \
               pressure_angle = 20, \
               helix_angle = 0, \
               tooth_step = 16, \
               flat_step = 10):
    return bevel_gear_assembly(modul, \
            tooth_number, \
            rad(partial_cone_angle), \
            tooth_width, \
            bore, \
            rad(pressure_angle), \
            rad(helix_angle), \
            tooth_step, \
            flat_step)

def bevel_herringbone_gear(modul, \
        tooth_number, \
        partial_cone_angle, \
        tooth_width, \
        bore = 0, \
        pressure_angle = 20, \
        helix_angle = 10, \
        tooth_step = 16, \
        flat_step = 10):
    return bevel_herringbone_gear_assembly(modul, \
            tooth_number, \
            rad(partial_cone_angle), \
            tooth_width, \
            bore, \
            rad(pressure_angle), \
            rad(helix_angle), \
            tooth_step, \
            flat_step)

def bevel_gear_pair(modul, \
        gear_teeth, \
        pinion_teeth, \
        tooth_width, \
        axis_angle = 90, \
        gear_bore = 0, \
        pinion_bore = 0, \
        pressure_angle = 20, \
        helix_angle = 0, \
        together_built = True, \
        tooth_step = 16, \
        flat_step = 10):
    return bevel_gear_pair_assembly(modul, \
            gear_teeth, \
            pinion_teeth, \
            rad(axis_angle), \
            tooth_width, \
            gear_bore, \
            pinion_bore, \
            rad(pressure_angle), \
            rad(helix_angle), \
            together_built, \
            tooth_step, \
            flat_step)

def bevel_herringbone_gear_pair(modul, \
        gear_teeth, \
        pinion_teeth, \
        tooth_width, \
        axis_angle = 90, \
        gear_bore = 0, \
        pinion_bore = 0, \
        pressure_angle = 20, \
        helix_angle = 10, \
        together_built = True, \
        tooth_step = 16, \
        flat_step = 10):
    return bevel_herringbone_gear_pair_assembly(modul, \
            gear_teeth, \
            pinion_teeth, \
            rad(axis_angle), \
            tooth_width, \
            gear_bore, \
            pinion_bore, \
            rad(pressure_angle), \
            rad(helix_angle), \
            together_built, \
            tooth_step, \
            flat_step)


