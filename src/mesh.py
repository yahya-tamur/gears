from .math import (
    center,
    interpolate_line,
    rotate,
    translate,
)

from .math3d import (
    dot3,
    minus3,
    cross,
    unit,
    distance,
    diameter,
    segment_point_distance,
    gram_schmidt3,
)
from math import pi, sin, cos
import struct

# 'closed' means the first and last points will be connected.


# def mesh_prism(top, bottom, closed=True, make_center=True):
#    ans = []

#    if closed:
#        ans += triangulate_polyhedron([top[-1], top[0], bottom[0], bottom[-1]])

#    for i in range(len(top) - 1):
#        ans += triangulate_polyhedron([top[i], top[i + 1], bottom[i + 1], bottom[i]])

#    return ans
# only to simplify code involving reverse


# factor = 0 -> point,
# factor = 1 -> center
# factor = -1 -> center + 2*(point - center)
def pull(point, center, factor):

    # point + factor*(center - point)
    # (1 - factor)point + factor*center
    a, b, c = point
    aa, bb, cc = center

    return (
        (1 - factor) * a + factor * aa,
        (1 - factor) * b + factor * bb,
        (1 - factor) * c + factor * cc,
    )


class Mesh:
    def __init__(self):
        self.mesh = []

    def save_stl(self, filename):
        def normal(t):
            u, v, w = t

            d1 = minus3(v, u)
            d2 = minus3(w, u)

            ans = cross(d1, d2)

            if dot3(ans, ans) == 0:
                print("colinear triangle:", t)
                return (1, 0, 0)

            return unit(ans)

        with open(filename, "wb") as f:
            f.write(bytearray(80))
            f.write(struct.pack("<i", len(self.mesh)))

            for tri in self.mesh:
                for ni in normal(tri):
                    f.write(struct.pack("<f", ni))

                for p in tri:
                    for pi in p:
                        f.write(struct.pack("<f", pi))

                f.write(bytearray(2))

    def add_triangle(self, a, b, c, reverse=False):
        if reverse:
            self.mesh.append([a, c, b])
        else:
            self.mesh.append([a, b, c])

    def add_polygon(self, p, reverse=False):

        center_point = center(p)

        for i in range(-1, len(p) - 1):
            self.add_triangle(center_point, p[i], p[i + 1], reverse=reverse)

    #  direction = 1 for inside, -1 for outside
    # this will try to move by max_delta, but only move by the
    # max diameter of the smallest segment
    # if full_move is True, if the max diameter of the smallest segment is
    # more than max_delta, instead of moving by max_delta, it will fail.
    # (1 + (direction)max delta) times = most scaling possible
    def scale_reduce(self, line, max_delta, direction, reverse=False):
        print("scale reduce", len(line))
        center_point = center(line)
        segments = []
        line.append(line[0])
        i = 0

        while i < len(line):
            segments.append(line[i : i + 4])
            i += 3

        delta = max_delta
        for segment in segments:
            c = center(segment)
            if len(segment) > 1:  # ?? and distance(c, center_point) > 0:
                d = diameter(segment)
                delta = min(delta, d / distance(c, center_point))

        line_ = []
        for i in range(len(segments)):
            c = pull(center(segments[i]), center_point, delta * direction)
            for j in range(len(segments[i]) - 1):
                self.add_triangle(
                    segments[i][j],
                    segments[i][j + 1],
                    c,
                    reverse=reverse ^ (direction == -1),
                )
            line_.append(c)
        for i in range(len(line_)):
            self.add_triangle(
                line_[i], line_[i - 1], line[i * 3], reverse=reverse ^ (direction == -1)
            )

        line.pop()
        print("scale reduce ans:", len(line_))
        return line_

    # nicely triangulates convex polygon made by line (set of coplanar points)
    # reduces large number of points
    def add_convex_polygon(self, line, reverse=False):
        while len(line) > 6:
            line = self.scale_reduce(line, 0.3, 1, reverse)
        self.add_polygon(line, reverse)

    def simple_add_ring(self, inner, outer, reverse=False):
        def closest(point, ls):
            _, ix = min((distance(point, ls[i]), i) for i in range(len(ls)))
            return ix

        i, j = 0, 0

        while (j != closest(inner[i], outer)) or (i != closest(outer[j], inner)):
            j = closest(inner[i], outer)
            i = closest(outer[j], inner)

        ii, jj = i, j
        i -= len(inner)
        j -= len(outer)
        while i < ii and j < jj:
            if distance(inner[i + 1], outer[j]) < distance(inner[i], outer[j + 1]):
                self.add_triangle(inner[i], outer[j], inner[i + 1], reverse)
                i += 1
            else:
                self.add_triangle(inner[i], outer[j], outer[j + 1], reverse)
                j += 1

        while i < ii:
            self.add_triangle(inner[i], outer[jj], inner[i + 1], reverse)
            i += 1
        while j < jj:
            self.add_triangle(inner[ii], outer[j], outer[j + 1], reverse)
            j += 1

    def add_ring(self, outer, inner, reverse=False):
        # centers should be equal...
        center_point = center(outer + inner)

        def r_inner():
            return max(distance(p, center_point) for p in inner)

        def r_outer():
            return min(
                segment_point_distance(outer[i], outer[i + 1], center_point)
                for i in range(-1, len(outer) - 1)
            )

        def ring_distance():
            return r_outer() - r_inner()

        def point_distance(ring):
            return max(distance(ring[i - 1], ring[i]) for i in range(-1, len(ring)))

        inner_extent = r_inner() + 0.4 * ring_distance()
        outer_extent = r_outer() - 0.4 * ring_distance()

        print(
            2 * point_distance(inner), inner_extent, point_distance(inner) / r_inner()
        )
        while (
            2 * point_distance(inner) < inner_extent - r_inner()
            and point_distance(inner) / r_inner() < 1
            and len(inner) > 18
        ):
            inner = self.scale_reduce(inner, inner_extent / r_inner() - 1, -1, reverse)

        while 2 * point_distance(outer) < r_outer() - outer_extent and len(outer) > 18:
            outer = self.scale_reduce(outer, 1 - outer_extent / r_outer(), 1, reverse)

        self.simple_add_ring(inner, outer, reverse)

    ## grid = i by j array of points
    def add_grid(self, grid, reverse=False):
        for i in range(len(grid) - 1):
            for j in range(len(grid[0]) - 1):
                self.add_triangle(grid[i][j], grid[i + 1][j], grid[i][j + 1], reverse)
                self.add_triangle(
                    grid[i + 1][j], grid[i + 1][j + 1], grid[i][j + 1], reverse
                )

    def add_stitch(self, first_line, second_line, n, reverse=False):
        grid = [interpolate_line(v, w, n) for (v, w) in zip(first_line, second_line)]
        self.add_grid(grid, reverse)
        return grid[0], grid[-1]

    def add_gear(
        self,
        tooth_west,
        tooth_east,
        tooth_number,
        flat_steps,
        bore,
        bore_steps,
        ring_gear,
    ):
        if ring_gear:
            for line in tooth_east:
                rotate([0, 0, 2 * pi / tooth_number], line)
            tooth_west, tooth_east = tooth_east, tooth_west

        teeth_west = [[line.copy() for line in tooth_west] for _ in range(tooth_number)]
        teeth_east = [[line.copy() for line in tooth_east] for _ in range(tooth_number)]

        for i in range(1, tooth_number):
            for line in teeth_west[i]:
                rotate([0, 0, 2 * pi * i / tooth_number], line)
            for line in teeth_east[i]:
                rotate([0, 0, 2 * pi * i / tooth_number], line)

        for tooth_west in teeth_west:
            self.add_grid(tooth_west, reverse=ring_gear ^ False)

        for tooth_east in teeth_east:
            self.add_grid(tooth_east, reverse=ring_gear ^ True)

        top_edge = []
        bottom_edge = []

        for i in range(-1, tooth_number - 1):
            bottom_between, top_between = self.add_stitch(
                [line[0] for line in teeth_east[i + 1]],
                [line[0] for line in teeth_west[i]],
                flat_steps,
                reverse=ring_gear,
            )
            bottom_tip, top_tip = self.add_stitch(
                [line[-1] for line in teeth_west[i]],
                [line[-1] for line in teeth_east[i]],
                flat_steps,
                reverse=ring_gear,
            )

            if ring_gear:
                bottom_base = interpolate_line(
                    teeth_east[i + 1][0][-1],
                    teeth_west[i][0][-1],
                    flat_steps,
                    endpoints=False,
                )
                top_base = interpolate_line(
                    teeth_east[i + 1][-1][-1],
                    teeth_west[i][-1][-1],
                    flat_steps,
                    endpoints=False,
                )
                # self.add_stitch(teeth_west[i][-1], teeth_east[i+1][-1], flat_steps)
                self.add_convex_polygon(
                    teeth_west[i][0]
                    + bottom_base[::-1]
                    + teeth_east[i + 1][0][::-1]
                    + bottom_between[1:-1],
                    reverse=True,
                )
                self.add_convex_polygon(
                    teeth_west[i][-1]
                    + top_base[::-1]
                    + teeth_east[i + 1][-1][::-1]
                    + top_between[1:-1],
                )
                top_edge += top_tip[::-1] + top_base[::-1]
                bottom_edge += bottom_tip[::-1] + bottom_base[::-1]
                # self.add_convex_polygon(
                #     teeth_west[i][0]
                #     + bottom_tip[1:-1]
                #     + teeth_east[i][0][::-1]
                #     + bottom_base,
                # )
                # top_edge += top_base + top_between[::-1]
                # bottom_edge += bottom_base + bottom_between[::-1]
            else:
                bottom_base = interpolate_line(
                    teeth_east[i][0][0],
                    teeth_west[i][0][0],
                    flat_steps,
                    endpoints=False,
                )
                top_base = interpolate_line(
                    teeth_east[i][-1][0],
                    teeth_west[i][-1][0],
                    flat_steps,
                    endpoints=False,
                )
                self.add_convex_polygon(
                    teeth_west[i][-1]
                    + top_tip[1:-1]
                    + teeth_east[i][-1][::-1]
                    + top_base,
                    reverse=True,
                )
                self.add_convex_polygon(
                    teeth_west[i][0]
                    + bottom_tip[1:-1]
                    + teeth_east[i][0][::-1]
                    + bottom_base,
                )
                top_edge += top_base + top_between[::-1]
                bottom_edge += bottom_base + bottom_between[::-1]

        if bore == 0:
            self.add_convex_polygon(top_edge)
            self.add_convex_polygon(bottom_edge, reverse=True)
        else:
            top_1, top_2, top_3 = center(top_edge)
            bottom_1, bottom_2, bottom_3 = center(bottom_edge)
            (v1, v2, v3), (w1, w2, w3) = gram_schmidt3(
                (top_1 - bottom_1, top_2 - bottom_2, top_3 - bottom_3),
                (
                    top_edge[0][0] - top_1,
                    top_edge[0][1] - top_2,
                    top_edge[0][2] - top_3,
                ),
                (
                    top_edge[1][0] - top_edge[0][0],
                    top_edge[1][1] - top_edge[0][1],
                    top_edge[2][2] - top_edge[0][2],
                ),
            )
            top_bore = []
            bottom_bore = []
            for i in range(bore_steps):
                t = i * 2 * pi / bore_steps
                top_bore.append(
                    (
                        top_1 + bore * cos(t) * v1 + bore * sin(t) * w1,
                        top_2 + bore * cos(t) * v2 + bore * sin(t) * w2,
                        top_3 + bore * cos(t) * v3 + bore * sin(t) * w3,
                    )
                )
                bottom_bore.append(
                    (
                        bottom_1 + bore * cos(t) * v1 + bore * sin(t) * w1,
                        bottom_2 + bore * cos(t) * v2 + bore * sin(t) * w2,
                        bottom_3 + bore * cos(t) * v3 + bore * sin(t) * w3,
                    )
                )

            if ring_gear:
                self.add_ring(top_bore, top_edge, reverse=False)
                self.add_ring(bottom_bore, bottom_edge, reverse=True)
            else:
                self.add_ring(top_edge, top_bore, reverse=False)
                self.add_ring(bottom_edge, bottom_bore, reverse=True)

            top_bore.append(top_bore[0])
            bottom_bore.append(bottom_bore[0])

            self.add_stitch(
                top_bore, bottom_bore, len(teeth_west[0]), reverse=ring_gear ^ False
            )

        n = 0
        z = 0
        for _, _, z_ in bottom_edge:
            z = (n / (n + 1)) * z + (1 / (n + 1)) * z_
            n += 1
        return z

    def translate(self, a):
        for tri in self.mesh:
            translate(a, tri)

    def rotate(self, a):
        for tri in self.mesh:
            rotate(a, tri)
