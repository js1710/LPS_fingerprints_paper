from __future__ import print_function, division
import unittest
import numpy as np
from membrane_thickness import *
from scipy.spatial import Delaunay


class MembraneThicknesTest(unittest.TestCase):
    def test_projected_thickness(self):
        triangle_xy = np.array([[0, 0], [1, 0], [1, 1]])
        triangle_z = np.array([4., 4., 4.])
        position = np.array([0.5, 0.5, 2])
        self.assertAlmostEqual(projected_thickness(position, triangle_xy, triangle_z), 2., delta=2/500.)

    def test_protein_overlap(self):
        p = np.arange(0, 1, 0.1)
        xx, yy = np.meshgrid(p, p)
        zz = np.zeros(xx.flatten().shape)
        zz.fill(10.)
        points = np.array([xx.flatten(), yy.flatten()]).T
        tri = Delaunay(points)
        protein_pos = np.array([[0.5, 0.5, 5], [0.5, 0.5, 15], [0.5, 0.5, 11], [0.5, 0.5, 9]])
        tri_offset = np.array([0, 0, 2.0])
        inds = protein_overlap_vectorised(points, zz, tri, protein_pos, tri_offset)
        self.assertListEqual(inds.tolist(), [2, 3])

    def test_volume_shells(self):
        dr = 5.0
        rmax = 10.
        bins = np.array([ 5.        ,  6.29960525,  7.21124785,  7.93700526,  8.54987973, 9.08560296,  9.56465591, 10.        ])
        np.testing.assert_almost_equal(get_volume_shells(dr, rmax), bins, decimal=4, verbose=True)


