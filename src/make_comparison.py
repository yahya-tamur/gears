from .lib import (
    make_gear,
    make_gear_pair,
    make_ring_gear,
    make_planetary_gear
)
import subprocess

from time import perf_counter
from sys import argv


# must be run from project directory
def run_openscad(command, output):
    command = f'OPENSCADPATH=.  openscad -o {output} <(echo "include <original/gears.scad>;\n {command};")'
    subprocess.run(command, shell=True, executable="/bin/bash")


def make_comparison(args):
    time_a, time_b = 0, 0

    if "flat_gear" in args or "all" in args:
        time_a -= perf_counter()
        # In this version, make_gear selects spur gear or
        # bevel gear based on partial cone angle.
        # In the openscad version, helix_angle for one
        # is equivalent to -helix angle in the other.
        # So, for flat gears, the helix angles in the two
        # versions are opposite. They are the same for bevel
        # gears.
        # The final model is also rotated a little differently.
        make_gear(
            modul=1,
            tooth_number=30,
            partial_cone_angle=0,
            tooth_width=7,
            bore=1.5,
            helix_angle=-40,
            filename="./comparison/flat_gear.stl",
        )
        time_a += perf_counter()

        time_b -= perf_counter()
        run_openscad(
            "herringbone_gear("
            + "modul=1,"
            + "bore=1.5,"
            + "tooth_number=30,"
            + "width=7,"
            + "helix_angle=40,"
            + "optimized=false"
            + ")",
            "./comparison/openscad_flat_gear.stl",
        )
        time_b += perf_counter()

    if "gear_pair" in args or "all" in args:
        time_a -= perf_counter()

        make_gear_pair(
            modul=0.7,
            gear_teeth=20,
            pinion_teeth=17,
            axis_angle=30,
            tooth_width=8,
            bore=0.5,
            helix_angle=-20,
            bore_steps=30,
            together_built=True,
            filename="./comparison/bevel_gear_pair.stl",
        )
        time_a += perf_counter()

        time_b -= perf_counter()
        run_openscad(
            "bevel_herringbone_gear_pair("
            + "modul=0.7,"
            + "gear_teeth=20,"
            + "pinion_teeth=17,"
            + "tooth_width=8,"
            + "axis_angle=30,"
            + "gear_bore=0.5,"
            + "pinion_bore=0.5,"
            + "helix_angle=-20,"
            + "together_built=true"
            + ")",
            "./comparison/openscad_bevel_gear_pair.stl",
        )
        time_b += perf_counter()

    if "ring_gear" in args or "all" in args:
        time_a -= perf_counter()
        make_ring_gear(
            modul=1.7,
            tooth_number=27,
            width=15,
            rim_width=15,
            pressure_angle=20,
            shortening_factor=0.2,
            helix_angle=30,
            filename="./comparison/ring_gear.stl",
        )
        time_a += perf_counter()

        time_b -= perf_counter()
        run_openscad(
            "herringbone_ring_gear("
            + "modul=1.7,"
            + "tooth_number=27,"
            + "width=15,"
            + "rim_width=15,"
            + "pressure_angle=20,"
            + "helix_angle=-30"
            + ")",
            "./comparison/openscad_ring_gear.stl",
        )
        time_b += perf_counter()

    if "planetary" in args or "all" in args:
        time_a -= perf_counter()
        make_planetary_gear(
            modul=2,
            sun_teeth=64,
            planet_teeth=16,
            width=30,
            number_planets=8,
            rim_width=5,
            bore=10,
            helix_angle=10,
            tooth_steps=8,
            flat_steps=2,
            helix_steps=20,
            bore_steps=50,
            ring_shortening_factor=0.6,
            together_built=True,
            filename="./comparison/planetary_gear.stl",
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
            + "helix_angle=-30,"
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
