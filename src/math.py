from math import sin, cos, tan, acos, pi, sqrt


def interpolate_line(a, b, n, endpoints=True):
    ans = []
    if endpoints:
        ans.append(a)
    for k in range(1, n):
        ans.append(tuple((a[i] * (n - k) + b[i] * k) / n for i in range(3)))
    if endpoints:
        ans.append(b)
    return ans


def center(points):
    n = 0
    ax, ay, az = 0, 0, 0
    for x, y, z in points:
        n += 1
        ax += (x - ax) / n
        ay += (y - ay) / n
        az += (z - az) / n
    return ax, ay, az


# same behavior as openscad rotate(a=...) but rotates a list of points
def rotate(a, pointlist):

    for i in range(len(pointlist)):
        x, y, z = pointlist[i]

        # rotate a[2] radians around z axis
        x, y = x * cos(a[2]) - y * sin(a[2]), x * sin(a[2]) + y * cos(a[2])

        # rotate a[1] radians around y axis
        z, x = z * cos(a[1]) - x * sin(a[1]), z * sin(a[1]) + x * cos(a[1])

        # rotate a[0] radians around x axis
        y, z = y * cos(a[0]) - z * sin(a[0]), y * sin(a[0]) + z * cos(a[0])

        pointlist[i] = (x, y, z)


def translate(a, pointlist):
    for i in range(len(pointlist)):
        x, y, z = pointlist[i]
        pointlist[i] = x + a[0], y + a[1], z + a[2]


def rad(t):
    return t / 180 * pi


def deg(t):
    return t / pi * 180


def polar_ev(r, rho):
    return (r / cos(rho), tan(rho) - rho)


def sphere_ev(t0, t):
    return acos(cos(t) / cos(t0)) / sin(t0) - acos(tan(t0) / tan(t))


def sph_to_cart(v):
    r, theta, phi = v
    return (r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta))


def pol_to_cart(r, theta, z=None):
    if z is None:
        return (r * cos(theta), r * sin(theta))
    return (r * cos(theta), r * sin(theta), z)

def diff(v, w):
    return (v[0] - w[0], v[1] - w[1], v[2] - w[2])

def distance(v, w):
    a, b, c = diff(v, w)
    return sqrt(a*a + b*b + c*c)

# line must have at least two points
def diameter(line):
    return max(max(distance(line[i], line[j]) for i in range(j)) for j in range(1, len(line)))

def segment_point_distance(p, q, center):
    px, py, pz = p
    qx, qy, qz = q
    cx, cy, cz = center
    t = ((cx - px)*(qx - px) + (cy - py)*(qy - py) + (cz - pz)*(qz - pz))/(distance(p, q) ** 2)
    if t <= 0:
        return distance(p, center)
    if t >= 1:
        return distance(q, center)
    w = px + t*(qx - px), py + t*(qy - py), pz + t*(qz - pz)
    return distance(w, center)