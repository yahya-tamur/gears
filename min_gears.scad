$fn = 50;
// General Variables
pi = 3.14159;
rad = 57.29578;
clearance = 0.05;   // clearance between teeth


/*  Sphere-Involutes-Function
    Returns the Azimuth Angle of an Involute Sphere
    theta0 = Polar Angle of the Cone, where the Cutting Edge of the Large Sphere unrolls the Involute
    theta = Polar Angle for which the Azimuth Angle of the Involute is to be calculated */
function sphere_ev(theta0,theta) = 1/sin(theta0)*acos(cos(theta)/cos(theta0))-acos(tan(theta0)/tan(theta));

/*  Converts Spherical Coordinates to Cartesian
    Format: radius, theta, phi; theta = Angle to z-Axis, phi = Angle to x-Axis on xy-Plane */
function sphere_to_cartesian(vect) = [
    vect[0]*sin(vect[1])*cos(vect[2]),
    vect[0]*sin(vect[1])*sin(vect[2]),
    vect[0]*cos(vect[1])
];

/*  Bevel Gear
    modul = Height of the Tooth Tip over the Partial Cone; Specification for the Outside of the Cone
    tooth_number = Number of Gear Teeth
    partial_cone_angle = (Half)angle of the Cone on which the other Ring Gear rolls
    tooth_width = Width of the Teeth from the Outside toward the apex of the Cone
    bore = Diameter of the Center Hole
    pressure_angle = Pressure Angle, Standard = 20° according to DIN 867. Should not exceed 45°.
    helix_angle = Helix Angle, Standard = 0° */
module bevel_gear(modul, tooth_number, partial_cone_angle, tooth_width, bore, pressure_angle = 20, helix_angle=0) {

    // Dimension Calculations
    d_outside = modul * tooth_number;                                    // Part Cone Diameter at the Cone Base,
                                                                    // corresponds to the Chord in a Spherical Section
    r_outside = d_outside / 2;                                        // Part Cone Radius at the Cone Base
    rg_outside = r_outside/sin(partial_cone_angle);                      // Large-Cone Radius for Outside-Tooth, corresponds to the Length of the Cone-Flank;
    rg_inside = rg_outside - tooth_width;                              // Large-Cone Radius for Inside-Tooth
    r_inside = r_outside*rg_inside/rg_outside;
    alpha_spur = atan(tan(pressure_angle)/cos(helix_angle));// Helix Angle in Transverse Section
    delta_b = asin(cos(alpha_spur)*sin(partial_cone_angle));          // Base Cone Angle
    da_outside = (modul <1)? d_outside + (modul * 2.2) * cos(partial_cone_angle): d_outside + modul * 2 * cos(partial_cone_angle);
    ra_outside = da_outside / 2;
    delta_a = asin(ra_outside/rg_outside);
    c = modul / 6;                                                  // Tip Clearance
    df_outside = d_outside - (modul +c) * 2 * cos(partial_cone_angle);
    rf_outside = df_outside / 2;
    delta_f = asin(rf_outside/rg_outside);
    rkf = rg_outside*sin(delta_f);                                   // Radius of the Cone Foot
    height_f = rg_outside*cos(delta_f);                               // Height of the Cone from the Root Cone

    echo("Part Cone Diameter at the Cone Base = ", d_outside);

    // Sizes for Complementary Truncated Cone
    height_k = (rg_outside-tooth_width)/cos(partial_cone_angle);          // Height of the Complementary Cone for corrected Tooth Length
    rk = (rg_outside-tooth_width)/sin(partial_cone_angle);               // Foot Radius of the Complementary Cone
    rfk = rk*height_k*tan(delta_f)/(rk+height_k*tan(delta_f));        // Tip Radius of the Cylinders for
                                                                    // Complementary Truncated Cone
    height_fk = rk*height_k/(height_k*tan(delta_f)+rk);                // height of the Complementary Truncated Cones

    echo("Bevel Gear Height = ", height_f-height_fk);

    phi_r = sphere_ev(delta_b, partial_cone_angle);                      // Angle to Point of Involute on Partial Cone

    // Torsion Angle gamma from Helix Angle
    gamma_g = 2*atan(tooth_width*tan(helix_angle)/(2*rg_outside-tooth_width));
    gamma = 2*asin(rg_outside/r_outside*sin(gamma_g/2));

    step = (delta_a - delta_b)/16;
    tau = 360/tooth_number;                                             // Pitch Angle
    start = (delta_b > delta_f) ? delta_b : delta_f;
    mirrpoint = (180*(1-clearance))/tooth_number+2*phi_r;

    if (delta_b > delta_f){
        echo("AAAA");
    else {
        echo("BBBB");
    }
    // Drawing
   // rotate([0,0,phi_r+90*(1-clearance)/tooth_number]){                      // Center Tooth on X-Axis;
                                                                    // Makes Alignment with other Gears easier
       // translate([0,0,height_f]) rotate(a=[0,180,0]){
            union(){
            //    translate([0,0,height_f]) rotate(a=[0,180,0]){                               // Truncated Cone
              //      difference(){
                //        linear_extrude(height=height_f-height_fk, scale=rfk/rkf) circle(rkf*1.001); // 1 permille Overlap with Tooth Root
                  //      translate([0,0,-1]){
                    //        cylinder(h = height_f-height_fk+2, r = bore/2);                // bore
                    //    }
                   // }
                //}
                for (rot = [0:tau:360]){
                    rotate (rot) {                                                          // Copy and Rotate "Number of Teeth"
                        union(){
                            if (delta_b > delta_f){
                                // Tooth Root
                                flankpoint_under = 1*mirrpoint;
                                flankpoint_over = sphere_ev(delta_f, start);
                                polyhedron(
                                    points = [
                                        sphere_to_cartesian([rg_outside, start*1.001, flankpoint_under]),    // 1 permille Overlap with Tooth
                                        sphere_to_cartesian([rg_inside, start*1.001, flankpoint_under+gamma]),
                                        sphere_to_cartesian([rg_inside, start*1.001, mirrpoint-flankpoint_under+gamma]),
                                        sphere_to_cartesian([rg_outside, start*1.001, mirrpoint-flankpoint_under]),
                                        sphere_to_cartesian([rg_outside, delta_f, flankpoint_under]),
                                        sphere_to_cartesian([rg_inside, delta_f, flankpoint_under+gamma]),
                                        sphere_to_cartesian([rg_inside, delta_f, mirrpoint-flankpoint_under+gamma]),
                                        sphere_to_cartesian([rg_outside, delta_f, mirrpoint-flankpoint_under])
                                    ],
                                    faces = [[0,1,2],[0,2,3],[0,4,1],[1,4,5],[1,5,2],[2,5,6],[2,6,3],[3,6,7],[0,3,7],[0,7,4],[4,6,5],[4,7,6]],
                                    convexity =1
                                );
                            }
                            // Tooth
                            for (delta = [start:step:delta_a-step]){
                                flankpoint_under = sphere_ev(delta_b, delta);
                                flankpoint_over = sphere_ev(delta_b, delta+step);
                                polyhedron(
                                    points = [
                                        sphere_to_cartesian([rg_outside, delta, flankpoint_under]),
                                        sphere_to_cartesian([rg_inside, delta, flankpoint_under+gamma]),
                                        sphere_to_cartesian([rg_inside, delta, mirrpoint-flankpoint_under+gamma]),
                                        sphere_to_cartesian([rg_outside, delta, mirrpoint-flankpoint_under]),
                                        sphere_to_cartesian([rg_outside, delta+step, flankpoint_over]),
                                        sphere_to_cartesian([rg_inside, delta+step, flankpoint_over+gamma]),
                                        sphere_to_cartesian([rg_inside, delta+step, mirrpoint-flankpoint_over+gamma]),
                                        sphere_to_cartesian([rg_outside, delta+step, mirrpoint-flankpoint_over])
                                    ],
                                    faces = [[0,1,2],[0,2,3],[0,4,1],[1,4,5],[1,5,2],[2,5,6],[2,6,3],[3,6,7],[0,3,7],[0,7,4],[4,6,5],[4,7,6]],
                                    convexity =1
                                );
                            }
                        }
                    }
                }
          //  }
       // }
    }
}

bevel_gear(modul=1, tooth_number=20, partial_cone_angle=20, tooth_width=10, bore=0, pressure_angle = 20, helix_angle=20);
