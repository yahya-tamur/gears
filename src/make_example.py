from .lib import (
    bevel_gear,
    bevel_herringbone_gear,
    bevel_gear_pair,
    bevel_herringbone_gear_pair,
)
from .make_stl import make_stl


def make_comparison():
    make_stl(
        bevel_gear(
            modul=1,
            tooth_number=20,
            partial_cone_angle=20,
            tooth_width=2,
            helix_angle=20,
        ),
        "./comparison/bevel_gear.stl",
    )

    make_stl(
        bevel_herringbone_gear(
            modul=1,
            tooth_number=20,
            partial_cone_angle=20,
            tooth_width=10,
            helix_angle=20,
            bore=2,
        ),
        "./comparison/bevel_herringbone_gear.stl",
    )

    make_stl(
        bevel_gear_pair(
            modul=1, gear_teeth=12, pinion_teeth=7, axis_angle=100, tooth_width=3
        ),
        "./comparison/bevel_gear_pair.stl",
    )


make_stl(
    bevel_herringbone_gear_pair(
        modul=2,
        gear_bore=2,
        pinion_bore=10,
        gear_teeth=40,
        pinion_teeth=22,
        tooth_width=25,
        axis_angle=70,
        helix_angle=40,
    ),
    "./comparison/bevel_herringbone_gear_pair.stl",
)


if __name__ == "__main__":
    make_comparison()
