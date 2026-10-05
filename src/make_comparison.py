from .lib import (
    bevel_gear,
    bevel_herringbone_gear,
    bevel_gear_pair,
    bevel_herringbone_gear_pair,
    herringbone_ring_gear,
    planetary_gear,
)
import subprocess

from .make_stl import make_stl

from time import perf_counter
from sys import argv


# must be run from project directory
def run_openscad(command, output):
    command = f'OPENSCADPATH=.  openscad -o {output} <(echo "include <original/gears.scad>; {command};")'
    subprocess.run(command, shell=True, executable="/bin/bash")


def make_comparison(args):
    time_a, time_b = 0, 0

    if "bevel" in args or "all" in args:
        time_a -= perf_counter()
        make_stl(
            bevel_gear(
                modul=1,
                bore=1.5,
                tooth_number=20,
                partial_cone_angle=30,
                tooth_width=2,
                helix_angle=40,
            ),
            "./comparison/bevel_gear.stl",
        )
        time_a += perf_counter()

        time_b -= perf_counter()
        run_openscad(
            "bevel_gear("
            + "modul=1,"
            + "bore=1.5,"
            + "tooth_number=20,"
            + "partial_cone_angle=30,"
            + "tooth_width=2,"
            + "helix_angle=40"
            + ")",
            "./comparison/openscad_bevel_gear.stl",
        )
        time_b += perf_counter()

    if "bevel_herringbone" in args or "all" in args:
        time_a -= perf_counter()
        make_stl(
            bevel_herringbone_gear(
                modul=2,
                tooth_number=21,
                partial_cone_angle=20,
                tooth_width=10,
                bore=0,
                pressure_angle=20,
                helix_angle=-20,
                tooth_step=3,
            ),
            "./comparison/bevel_herringbone_gear.stl",
        )
        time_a += perf_counter()

        time_b -= perf_counter()
        run_openscad(
            "bevel_herringbone_gear("
            + "modul=2,"
            + "bore=0,"
            + "tooth_number=21,"
            + "partial_cone_angle=20,"
            + "tooth_width=10,"
            + "helix_angle=-20"
            + ")",
            "./comparison/openscad_bevel_herringbone_gear.stl",
        )
        time_b += perf_counter()

    if "bevel_pair" in args or "all" in args:
        time_a -= perf_counter()
        make_stl(
            bevel_gear_pair(
                modul=1,
                gear_teeth=50,
                pinion_teeth=18,
                tooth_width=6,
                axis_angle=90,
                gear_bore=0,
                pinion_bore=1.1,
                pressure_angle=15,
                helix_angle=40,
                together_built=True,
                tooth_step=16,
            ),
            "./comparison/bevel_gear_pair.stl",
        )
        time_a += perf_counter()

        time_b -= perf_counter()
        run_openscad(
            "bevel_gear_pair("
            + "modul=1,"
            + "gear_teeth=50,"
            + "pinion_teeth=18,"
            + "tooth_width=6,"
            + "axis_angle=90,"
            + "gear_bore=0,"
            + "pinion_bore=1.1,"
            + "pressure_angle=15,"
            + "helix_angle=40,"
            + "together_built=true"
            + ")",
            "./comparison/openscad_bevel_gear_pair.stl",
        )
        time_b += perf_counter()

    if "bevel_herringbone_pair" in args or "all" in args:
        time_a -= perf_counter()
        make_stl(
            bevel_herringbone_gear_pair(
                modul=0.7,
                gear_teeth=70,
                pinion_teeth=47,
                tooth_width=2,
                axis_angle=30,
                gear_bore=0.5,
                pinion_bore=0,
                pressure_angle=25,
                helix_angle=-10,
                together_built=True,
                tooth_step=16,
            ),
            "./comparison/bevel_herringbone_gear_pair.stl",
        )
        time_a += perf_counter()

        time_b -= perf_counter()
        run_openscad(
            "bevel_herringbone_gear_pair("
            + "modul=0.7,"
            + "gear_teeth=70,"
            + "pinion_teeth=47,"
            + "tooth_width=2,"
            + "axis_angle=30,"
            + "gear_bore=0.5,"
            + "pinion_bore=0,"
            + "pressure_angle=25,"
            + "helix_angle=-10,"
            + "together_built=true"
            + ")",
            "./comparison/openscad_bevel_herringbone_gear_pair.stl",
        )
        time_b += perf_counter()

    if "herringbone_ring" in args or "all" in args:
        time_a -= perf_counter()
        make_stl(
            herringbone_ring_gear(
                modul=1.7,
                tooth_number=27,
                width=47,
                rim_width=15,
                pressure_angle=20,
                helix_angle=30,
                shortening_factor=0.6,
                tooth_step=26,
            ),
            "./comparison/herringbone_ring_gear.stl",
        )
        time_a += perf_counter()

        time_b -= perf_counter()
        run_openscad(
            "herringbone_ring_gear("
            + "modul=1.7,"
            + "tooth_number=27,"
            + "width=47,"
            + "rim_width=15,"
            + "pressure_angle=20,"
            + "helix_angle=30"
            + ")",
            "./comparison/openscad_herringbone_ring_gear.stl",
        )
        time_b += perf_counter()

    if "planetary" in args or "all" in args:
        time_a -= perf_counter()
        make_stl(
            planetary_gear(
                modul=2,
                sun_teeth=64,
                planet_teeth=16,
                width=30,
                number_planets=8,
                rim_width=5,
                sun_bore=10,
                planet_bore=10,
                pressure_angle=20,
                helix_angle=30,
                together_built=True,
                tooth_step=16,
                ring_shortening_factor=1,
            ),
            "./comparison/planetary_gear.stl",
        )
        time_a += perf_counter()

        time_b -= perf_counter()
        run_openscad(
            "planetary_gear("
            + "modul=2,"
            + "sun_teeth=64,"
            + "planet_teeth=16,"
            + "number_planets=8,"
            + "width=30,"
            + "rim_width=5,"
            + "bore=10,"
            + "pressure_angle=20,"
            + "helix_angle=30,"
            + "together_built=true,"
            + "optimized=false"
            + ")",
            "./comparison/openscad_planetary_gear.stl",
        )
        time_b += perf_counter()

    print(f"Total time (python): {time_a:.6f} seconds")
    print(f"Total time (openscad): {time_b:.6f} seconds")


if __name__ == "__main__":
    make_comparison(set(argv[1:]))
