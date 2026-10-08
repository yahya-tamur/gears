include <gears.scad>

//bevel_gear(modul=1, tooth_number=20, partial_cone_angle=20, tooth_width=2, helix_angle=20);

//

//herringbone_gear(modul=2, tooth_number=20, width=10, helix_angle=30);
//bevel_herringbone_gear(modul=2, partial_cone_angle=0.1, tooth_number=20, width=10, helix_angle=30);

//bevel_herringbone_gear(modul=2, tooth_number=20, width=1, partial_cone_angle=30);
//bevel_herringbone_gear_pair(
//    modul=2,
//    gear_teeth=20,
//    pinion_teeth=21,
//    axis_angle=0,
//    tooth_width=10,
//)//bevel_gear(modul=1, tooth_number=20, partial_cone_angle=20, tooth_width=2, helix_angle=20);
//bevel_gear_pair(modul=1, gear_teeth=12, pinion_teeth=7, axis_angle=100, tooth_width=3);

//bevel_herringbone_gear_pair(modul=2, gear_teeth=40, pinion_teeth=22, tooth_width=25, axis_angle=70, helix_angle=40);

//herringbone_ring_gear(
        modul=2,
        tooth_number=26,
        pressure_angle=20, 
        helix_angle=15,width=10,rim_width=5);
        
//herringbone_gear(modul=2, tooth_number=20, width=20, bore=0);

planetary_gear(modul=2, sun_teeth=16, planet_teeth=16, number_planets=4, width=10, rim_width=4, bore=0, helix_angle=20);

//herringbone_ring_gear(modul=2, tooth_number=20, width=10, rim_width=5, pressure_angle=20, helix_angle=10);