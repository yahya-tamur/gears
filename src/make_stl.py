import struct

from math import sqrt  # for normal vector
from .math import center

# The two functions below are included for convenience, make_stl can put any
# array of triangles (optionally including a normal vector) into an stl file.


# In case the direction of the normal matters, both functions below create a
# consistent direction, which can reversed:
# triangulate_polyhedron, using the 'reverse' argument, and
# triangulate_prism, by exchanging the top and bottom arguments.

# This can be used to stitch together any two sequences of points, not just the
# sides of a prism. 'closed' means the first and last points will be connected.


def triangulate_prism(top, bottom, closed=True, make_center=True):
    ans = []

    if closed:
        ans += triangulate_polyhedron([top[-1], top[0], bottom[0], bottom[-1]])

    for i in range(len(top) - 1):
        ans += triangulate_polyhedron([top[i], top[i + 1], bottom[i + 1], bottom[i]])

    return ans


# All triangles have a common center. If center=False, this is the first
# element of the list. Else, this is calculated.


def triangulate_polyhedron(p, reverse=False):

    center_point = center(p)

    if reverse:
        ans = [[center_point, p[0], p[-1]]]
    else:
        ans = [[center_point, p[-1], p[0]]]

    for i in range(0, len(p) - 1):
        if reverse:
            ans.append([center_point, p[i + 1], p[i]])
        else:
            ans.append([center_point, p[i], p[i + 1]])

    return ans


def normal(t):
    (x0, y0, z0), (x1, y1, z1), (x2, y2, z2) = t

    d1x, d1y, d1z = x1 - x0, y1 - y0, z1 - z0
    d2x, d2y, d2z = x2 - x0, y2 - y0, z2 - z0

    ans_x, ans_y, ans_z = (
        d1y * d2z - d1z * d2y,
        -d1x * d2z + d1z * d2x,
        d1x * d2y - d1y * d2x,
    )

    n = sqrt(ans_x * ans_x + ans_y * ans_y + ans_z * ans_z)

    return (ans_x / n, ans_y / n, ans_z / n)


def make_stl(triangles, filename):
    with open(filename, "wb") as f:
        f.write(bytearray(80))
        f.write(struct.pack("<i", len(triangles)))

        for tri in triangles:
            if len(tri) == 3:
                for ni in normal(tri):
                    f.write(struct.pack("<f", ni))

            for p in tri:
                for pi in p:
                    f.write(struct.pack("<f", pi))
            f.write(bytearray(2))


cube = [
    [[-1, 0, 0], [0, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, 1, 1], [0, 1, 0], [0, 0, 1]],
    [[1, 0, 0], [1, 0, 0], [1, 1, 0], [1, 0, 1]],
    [[1, 0, 0], [1, 1, 1], [1, 1, 0], [1, 0, 1]],
    [[0, -1, 0], [0, 0, 0], [1, 0, 0], [0, 0, 1]],
    [[0, -1, 0], [1, 0, 1], [1, 0, 0], [0, 0, 1]],
    [[0, 1, 0], [0, 1, 0], [1, 1, 0], [0, 1, 1]],
    [[0, 1, 0], [1, 1, 1], [1, 1, 0], [0, 1, 1]],
    [[0, 0, -1], [0, 0, 0], [1, 0, 0], [0, 1, 0]],
    [[0, 0, -1], [1, 1, 0], [1, 0, 0], [0, 1, 0]],
    [[0, 0, 1], [0, 0, 1], [1, 0, 1], [0, 1, 1]],
    [[0, 0, 1], [1, 1, 1], [1, 0, 1], [0, 1, 1]],
]

make_stl(cube, "cube.stl")
