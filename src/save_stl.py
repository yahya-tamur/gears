import struct

from math import sqrt


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
    if n == 0:
        print("colinear triangle:", t)
        return (1, 0, 0)

    return (ans_x / n, ans_y / n, ans_z / n)


def save_stl(triangles, filename):
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
