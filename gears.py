
# get two arrays of points (same length >= 3)
# no top or bottom.
def triangulate_prism(top, bottom):
    ans = []
    

    ans.append([top[-1],bottom[-1], top[0]]) 
    ans.append([top[0],bottom[0], bottom[-1]]) 

    for i in range(len(top)-1):
        ans.append([top[i], bottom[i], top[i+1]])
        ans.append([top[i+1], bottom[i], bottom[i+1]])
    print(len(top), len(bottom), len(ans))
    return ans

def triangulate_polyhedron(p, center=None):
    start = 0
    if center is None:
        center = p[0]
        start = 1

    ans = []#[[center, p[-1], p[0]]]

    for i in range(start, len(p)-1):
        ans.append([center, p[i], p[i+1]])

    return ans


from math import sin, cos, tan, atan, asin, acos, pi

# all angles are in RADIANS here!!!!

def rad(t):
    return t/180*pi

def deg(t):
    return t/pi*180


def sphere_ev(t0, t):
    return acos(cos(t)/cos(t0))/sin(t0) - acos(tan(t0)/tan(t));

def sph_to_cart(v):
    r, theta, phi = v
    return [r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta)]

def bevel_gear(modul, tooth_number, partial_cone_angle, tooth_width, pressure_angle = rad(20), helix_angle=0):
    
    ans = []

    clearance = 0.05
    d_outside = modul*tooth_number
    r_outside = d_outside / 2
    rg_outside = r_outside / sin(partial_cone_angle)
    rg_inside = rg_outside - tooth_width
    r_inside = r_outside*rg_inside / rg_outside
    alpha_spur = atan(tan(pressure_angle)/cos(helix_angle))
    delta_b = asin(cos(alpha_spur)*sin(partial_cone_angle))
    da_outside = d_outside + (modul * (2.2 if modul < 1 else 2)) * cos(partial_cone_angle)
    ra_outside = da_outside / 2
    delta_a = asin(ra_outside/rg_outside)
    c = modul / 6
    df_outside = d_outside - (modul +c) * 2 * cos(partial_cone_angle)
    rf_outside = df_outside / 2
    delta_f = asin(rf_outside/rg_outside)
    rkf = rg_outside*sin(delta_f)
    height_f = rg_outside*cos(delta_f)
    height_k = (rg_outside-tooth_width)/cos(partial_cone_angle)
    rk = (rg_outside-tooth_width)/sin(partial_cone_angle)
    rfk = rk*height_k*tan(delta_f)/(rk+height_k*tan(delta_f))
    height_fk = rk*height_k/(height_k*tan(delta_f)+rk)
    phi_r = sphere_ev(delta_b, partial_cone_angle)

    gamma_g = 2*atan(tooth_width*tan(helix_angle)/(2*rg_outside-tooth_width))
    gamma = 2*asin(rg_outside/r_outside*sin(gamma_g/2))

    step = (delta_a - delta_b)/16
    tau = 2*pi/tooth_number

    start = delta_b if (delta_b > delta_f) else delta_f
    mirrpoint = (pi*(1-clearance))/tooth_number+2*phi_r

    #for rot in range(0, 2*pi, tau):
    # one tooth

    tooth_top_from_left = []
    tooth_top_from_right = []
    tooth_bottom_from_left = []
    tooth_bottom_from_right = []

    if delta_b > delta_f:
        flankpoint_under = 1*mirrpoint

        tooth_top_from_right.append(sph_to_cart([rg_outside, delta_f, flankpoint_under]))
        tooth_bottom_from_right.append(sph_to_cart([rg_inside, delta_f, flankpoint_under+gamma]))
        tooth_bottom_from_left.append(sph_to_cart([rg_inside, delta_f, mirrpoint-flankpoint_under+gamma]))
        tooth_top_from_left.append(sph_to_cart([rg_outside, delta_f, mirrpoint-flankpoint_under]))

    #for delta in range(start, delta_a, step):
    delta = start
    #print(deg(start), deg(delta_a), deg(step))
    #print(delta_a, ra_outside, rg_outside)
    while delta < delta_a+(step/2):
        #print(delta_b, delta)
        flankpoint_under = sphere_ev(delta_b, delta)
        #print(deg(flankpoint_under), deg(gamma))

        tooth_top_from_left.append(sph_to_cart([rg_outside, delta, flankpoint_under]))
        tooth_bottom_from_left.append(sph_to_cart([rg_inside, delta, flankpoint_under+gamma]))
        tooth_bottom_from_right.append(sph_to_cart([rg_inside, delta, mirrpoint-flankpoint_under+gamma]))
        tooth_top_from_right.append(sph_to_cart([rg_outside, delta, mirrpoint-flankpoint_under]))

        delta += step
    tooth_top = tooth_top_from_left + tooth_top_from_right[::-1]
    tooth_bottom = tooth_bottom_from_left + tooth_bottom_from_right[::-1]

    for t in tooth_top:
        print(t)
    print('='*50)
    for t in tooth_bottom:
        print(t)
    #print(tooth_bottom)
    print(len(tooth_top))
    print(len(tooth_bottom))
    ans += triangulate_polyhedron(tooth_top)
    ans += triangulate_polyhedron(tooth_bottom)
    ans += triangulate_prism(tooth_top, tooth_bottom)
    #for t in ans:
     #   print(t)

    return ans

