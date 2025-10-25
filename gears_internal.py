
# try mpmath?
#from math import sin, cos, tan, atan, asin, acos, pi, sqrt

from mpmath import mp
mp.dps = 500
sin = mp.sin
cos = mp.cos
tan = mp.tan
atan = mp.atan
asin = mp.asin
acos = mp.acos
pi = mp.pi
sqrt = mp.sqrt

# get two arrays of points (same length >= 3)
# no top or bottom.
def triangulate_prism(top, bottom, open=False):
    ans = []
    

    if not open:
        ans.append([top[-1],bottom[-1], top[0]]) 
        ans.append([top[0],bottom[0], bottom[-1]]) 

    for i in range(len(top)-1):
        ans.append([top[i], bottom[i], top[i+1]])
        ans.append([top[i+1], bottom[i], bottom[i+1]])
    return ans

def triangulate_polyhedron(p, center=None, reverse=False):
    start = 0
    if center is None:
        center = p[0]
        start = 1

    ans = []#[[center, p[-1], p[0]]]

    for i in range(start, len(p)-1):
        if reverse:
            ans.append([center, p[i+1], p[i]])
        else:
            ans.append([center, p[i], p[i+1]])

    return ans

def interpolate_line(a, b, n, endpoints=True):
    ans = []
    if endpoints:
        ans.append(a)
    for k in range(1, n):
        ans.append(tuple((a[i]*(n-k) + b[i]*k)/n for i in range(3)))
    if endpoints:
        ans.append(b)
    return ans

def center(points):
    n = 0
    ax, ay, az = 0, 0, 0
    for (x, y, z) in points:
        n += 1
        ax += (x-ax)/n
        ay += (y-ay)/n
        az += (z-az)/n
    return ax, ay, az

def norm(x, y, z):
    return sqrt(x*x + y*y + z*z)

def project(points, center, dist):
    cx, cy, cz = center
    ans = []
    for (x, y, z) in points:
        v = norm(x-cx, y-cy, z-cz)
        x_ = cx + dist*(x-cx)/v
        y_ = cy + dist*(y-cy)/v
        z_ = cz + dist*(z-cz)/v
        ans.append((x_, y_, z_))
    return ans


#same behavior as openscad rotate(a=...) but rotates a list of points
def rotate(a, l):

    for i in range(len(l)):
        x, y, z = l[i]

        # rotate a[2] radians around z axis
        x, y = x*cos(a[2]) - y*sin(a[2]), x*sin(a[2]) + y*cos(a[2])

        #rotate a[1] radians around y axis
        x, z = x*cos(a[1]) - z*sin(a[1]), x*sin(a[1]) + z*cos(a[1])

        #rotate a[0] radians around x axis
        y, z = y*cos(a[0]) - z*sin(a[0]), y*sin(a[0]) + z*cos(a[0])
    
        l[i] = (x, y, z)

def translate(a, l):
    for i in range(len(l)):
        x, y, z = l[i]
        l[i] = x+a[0], y+a[1], z+a[2]

# all angles are in RADIANS here!!!!

def rad(t):
    return t/180*pi

def deg(t):
    return t/pi*180


def sphere_ev(t0, t):
    return acos(cos(t)/cos(t0))/sin(t0) - acos(tan(t0)/tan(t));

def sph_to_cart(v):
    r, theta, phi = v
    return (r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta))


def printline(prefix, line):
    print(prefix, end='')

    for (x, y, z) in line:
        print('(%.2f %.2f %.2f)' % (x, y, z), end='')
    print()


def bevel_gear_data(modul, tooth_number, partial_cone_angle, tooth_width, pressure_angle, helix_angle, tooth_step):

    clearance = 0.05
    d_outside = modul*tooth_number
    r_outside = d_outside / 2
    rg_outside = r_outside / sin(partial_cone_angle)
    rg_inside = rg_outside - tooth_width
    r_inside = r_outside*rg_inside / rg_outside
    alpha_spur = atan(tan(pressure_angle)/cos(helix_angle))
    da_outside = d_outside + (modul * (2.2 if modul < 1 else 2)) * cos(partial_cone_angle)
    ra_outside = da_outside / 2
    c = modul / 6
    df_outside = d_outside - (modul +c) * 2 * cos(partial_cone_angle)
    rf_outside = df_outside / 2
    #rkf = rg_outside*sin(delta_f)
    delta_f = asin(rf_outside/rg_outside)
    delta_a = asin(ra_outside/rg_outside)
    delta_b = asin(cos(alpha_spur)*sin(partial_cone_angle))

    height_f = rg_outside*cos(delta_f)

    height_k = (rg_outside-tooth_width)/cos(partial_cone_angle)
    rk = (rg_outside-tooth_width)/sin(partial_cone_angle)
    #rfk = rk*height_k*tan(delta_f)/(rk+height_k*tan(delta_f))
    #height_fk = rk*height_k/(height_k*tan(delta_f)+rk)


    phi_r = sphere_ev(delta_b, partial_cone_angle)

    gamma_g = 2*atan(tooth_width*tan(helix_angle)/(2*rg_outside-tooth_width))
    gamma = 2*asin(rg_outside/r_outside*sin(gamma_g/2))

    tau = 2*pi/tooth_number

    mirrpoint = (pi*(1-clearance))/tooth_number+2*phi_r

    #for rot in range(0, 2*pi, tau):
    # one tooth

    tooth_nw = []
    tooth_ne = []
    tooth_sw = []
    tooth_se = []



    start = delta_f
    step = (delta_a - delta_f)/tooth_step # check if this is good
    if delta_b > delta_f:
        flankpoint_under = 1*mirrpoint

        tooth_ne.append(sph_to_cart((rg_outside, delta_f, flankpoint_under)))
        tooth_se.append(sph_to_cart((rg_inside, delta_f, flankpoint_under+gamma)))
        tooth_sw.append(sph_to_cart((rg_inside, delta_f, mirrpoint-flankpoint_under+gamma)))
        tooth_nw.append(sph_to_cart((rg_outside, delta_f, mirrpoint-flankpoint_under)))

        start = delta_b
        step = (delta_a - delta_b)/tooth_step

    #for delta in range(start, delta_a, step):
    delta = start
    while delta < delta_a+(step/2):
        flankpoint_under = sphere_ev(delta_b, delta)

        tooth_nw.append(sph_to_cart((rg_outside, delta, flankpoint_under)))
        tooth_sw.append(sph_to_cart((rg_inside, delta, flankpoint_under+gamma)))
        tooth_se.append(sph_to_cart((rg_inside, delta, mirrpoint-flankpoint_under+gamma)))
        tooth_ne.append(sph_to_cart((rg_outside, delta, mirrpoint-flankpoint_under)))

        delta += step

    for pt_list in (tooth_nw, tooth_ne, tooth_sw, tooth_se):
        rotate([0,pi,0], pt_list)
        translate([0,0,height_f], pt_list)
        rotate([0,0,phi_r+pi/2*(1-clearance)/tooth_number], pt_list)

    return (tooth_nw, tooth_ne, tooth_sw, tooth_se, -tau)

# do translation/rotation next
def bevel_gear_assembly(modul, tooth_number, partial_cone_angle, tooth_width, bore, pressure_angle, helix_angle, tooth_step, flat_step):
    tooth_nw, tooth_ne, tooth_sw, tooth_se, tau = bevel_gear_data(modul, tooth_number, partial_cone_angle, tooth_width, pressure_angle, helix_angle, tooth_step)

    tooth_top = tooth_nw + interpolate_line(tooth_nw[-1], tooth_ne[-1], flat_step, endpoints=False) + tooth_ne[::-1]
    tooth_bottom = tooth_sw + interpolate_line(tooth_sw[-1], tooth_se[-1], flat_step, endpoints=False) + tooth_se[::-1]

    top_face, bottom_face = [], []

    ans = []

    # create top and bottom teeth faces, teeth open prisms

    i = 0
    while True:
        ans += triangulate_polyhedron(tooth_top)
        ans += triangulate_polyhedron(tooth_bottom, reverse=True)
        ans += triangulate_prism(tooth_top, tooth_bottom, open=True)

        if len(top_face) > 0:
            top_line = interpolate_line(top_face[-1], tooth_top[0], flat_step, endpoints=True)
            bottom_line = interpolate_line(bottom_face[-1], tooth_bottom[0], flat_step, endpoints=True)
            ans += triangulate_prism(top_line, bottom_line, open=True)
            top_face.pop()
            bottom_face.pop()
            top_face += top_line
            bottom_face += bottom_line
            top_face.append(tooth_top[-1])
            bottom_face.append(tooth_bottom[-1])
        else:
            top_face = [tooth_top[0], tooth_top[-1]]
            bottom_face = [tooth_bottom[0], tooth_bottom[-1]]

        i += 1
        if i == tooth_number:
            break
        rotate((0, 0, tau), tooth_top)
        rotate((0, 0, tau), tooth_bottom)

    top_line = interpolate_line(top_face[-1], top_face[0], flat_step, endpoints=True)
    bottom_line = interpolate_line(bottom_face[-1], bottom_face[0], flat_step, endpoints=True)
    ans += triangulate_prism(top_line, bottom_line, open=True)

    top_line.pop()
    bottom_line.pop()
    top_face.pop()
    bottom_face.pop()
    top_face += top_line
    bottom_face += bottom_line

    top_center = center(top_face)
    bottom_center = center(bottom_face)

    if bore == 0:
        ans += triangulate_polyhedron(top_face, top_center)
        ans += triangulate_polyhedron(bottom_face, top_center)

    else:
        top_bore, bottom_bore = [], []
        for k in range(len(top_face)):
            top_bore.append((-cos(2*pi*k/len(top_face)), sin(2*pi*k/len(top_face)), top_face[k][2]))
            bottom_bore.append((-cos(2*pi*k/len(top_face)), sin(2*pi*k/len(top_face)), bottom_face[k][2]))


        ans += triangulate_prism(top_bore, top_face, open=False)
        ans += triangulate_prism(bottom_bore, top_bore, open=False)
        ans += triangulate_prism(bottom_face, bottom_bore, open=False)

    return ans

def bevel_herringbone_gear_data(modul, tooth_number, partial_cone_angle, tooth_width, bore, pressure_angle, helix_angle, tooth_step):

    tooth_width = tooth_width / 2
    d_outside = modul * tooth_number
    r_outside = d_outside / 2
    rg_outside = r_outside/sin(partial_cone_angle)
    c = modul / 6
    df_outside = d_outside - (modul +c) * 2 * cos(partial_cone_angle)
    rf_outside = df_outside / 2
    delta_f = asin(rf_outside/rg_outside)
    height_f = rg_outside*cos(delta_f)

    gamma_g = 2*atan(tooth_width*tan(helix_angle)/(2*rg_outside-tooth_width))
    gamma = 2*asin(rg_outside/r_outside*sin(gamma_g/2))

    height_k = (rg_outside-tooth_width)/cos(partial_cone_angle)
    rk = (rg_outside-tooth_width)/sin(partial_cone_angle)
    rfk = rk*height_k*tan(delta_f)/(rk+height_k*tan(delta_f))
    height_fk = rk*height_k/(height_k*tan(delta_f)+rk)

    modul_inside = modul*(1-tooth_width/rg_outside)
    
    lower_cone_angle = partial_cone_angle# - pi/180

    # I was considering altering the bevel_gear_data program to be able to create
    # the lower mesh only, but since it's O(tooth_step), which is 16 by default,
    # it's really not necessary.
    

    tooth_aw, tooth_ae, tooth_bw_upper, tooth_be_upper, tau = bevel_gear_data( \
            modul, tooth_number, lower_cone_angle, tooth_width, pressure_angle, helix_angle, tooth_step)

    tooth_top = tooth_aw + tooth_ae[::-1]
    print(len(tooth_top))

    tooth_mid_upper = tooth_bw_upper + tooth_be_upper[::-1]


    tooth_bw_lower, tooth_be_lower, tooth_cw, tooth_ce, tau = bevel_gear_data( \
            modul_inside, tooth_number, partial_cone_angle, tooth_width, pressure_angle, -helix_angle, tooth_step)

    tooth_mid_lower = tooth_bw_lower + tooth_be_lower[::-1]
    tooth_bottom = tooth_cw + tooth_ce[::-1]

    for pt_list in (tooth_mid_lower, tooth_bottom):
        rotate([0, 0, -gamma], pt_list)
        translate([0, 0, height_f - height_fk], pt_list)

    print(sum(norm(a[0]-b[0], a[1]-b[1], a[2]-b[2]) for (a,b) in zip(tooth_mid_upper, tooth_mid_lower))/len(tooth_mid_upper))
    ans = []

    ans += triangulate_polyhedron(tooth_top)#tooth_ae[::-1])
    ans += triangulate_polyhedron(tooth_mid_upper)# + tooth_be_upper[::-1])
    ans += triangulate_prism(tooth_top, tooth_mid_upper, open=True)
    ans += triangulate_polyhedron(tooth_mid_lower)# + tooth_be_lower[::-1])
    ans += triangulate_polyhedron(tooth_bottom)# + tooth_ce[::-1])
    ans += triangulate_prism(tooth_mid_lower, tooth_bottom, open=True)

    print(len(tooth_mid_upper), len(tooth_mid_lower))
    return ans
