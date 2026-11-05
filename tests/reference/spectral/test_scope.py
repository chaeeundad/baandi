#!/usr/bin/env python3
"""Frozen-cap and degeneracy matching checks for the stored-evidence analyzer."""
import importlib.util
from pathlib import Path
import unittest
spec=importlib.util.spec_from_file_location('scope',Path(__file__).with_name('scope.py'))
s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
class Matching(unittest.TestCase):
    def test_degenerate_boundary_is_not_missing(self):
        r=s.compare([-3.,-1.],[-3.,-1.,-1.,2.],[2.,0.])
        self.assertEqual(r['missing_dense_indices_zero_based'],[])
    def test_degenerate_occupied_member_missing(self):
        r=s.compare([-3.,-1.],[-3.,-3.,-1.,2.],[2.,2.])
        self.assertEqual(r['missing_dense_indices_zero_based'],[1])
        self.assertEqual(r['missing_clusters'],[[0,1]])
        self.assertEqual(r['band_sum_excess_ha'],4.)
    def test_loose_residual_not_a_missing_root_proof(self):
        r=s.compare([-2.99,-1.],[-3.,-1.,2.],[2.,0.])
        self.assertIsNone(r['missing_dense_indices_zero_based'])
        self.assertFalse(r['missing_indices_certified'])
    def test_caps_are_frozen(self):
        self.assertEqual(s.CAP,1e-9)
if __name__=='__main__':unittest.main()
