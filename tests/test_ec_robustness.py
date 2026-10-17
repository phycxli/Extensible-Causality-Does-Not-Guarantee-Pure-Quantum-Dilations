"""Boundary and scope rejection controls for the neighborhood certificate."""

from fractions import Fraction
import unittest

from repro.independent_ec_robustness_audit import (
    certify_ball, ordered_split_coefficients, valid_pauli_word,
)


class ECRobustnessTests(unittest.TestCase):
    def test_closed_stated_ball_has_strict_witness_margin(self):
        result = certify_ball("1/22", "1/264")
        self.assertEqual(result["CS_component_minimum_eigenvalue_bound"], "0")
        self.assertEqual(result["perturbed_witness_upper_bound"], "-13/660")

    def test_larger_radius_is_not_certified_by_the_CS_inner_ball(self):
        with self.assertRaisesRegex(ValueError, "sufficient EC bound"):
            certify_ball("1/22", "1/263")

    def test_zero_witness_margin_is_not_strict_separation(self):
        with self.assertRaisesRegex(ValueError, "strict non-CS"):
            certify_ball("1/12", "1/960")

    def test_original_point_is_not_asserted_to_be_an_interior_point(self):
        with self.assertRaisesRegex(ValueError, "positive ball"):
            certify_ball(0, Fraction(1, 10000))

    def test_noise_CS_boundary_is_excluded(self):
        with self.assertRaisesRegex(ValueError, "positive ball"):
            certify_ball("1/11", "1/10000")

    def test_symmetry_breaking_input_word_is_a_valid_direction(self):
        self.assertTrue(valid_pauli_word("YIYI"))
        self.assertEqual(ordered_split_coefficients("YIYI"), (Fraction(1, 2), Fraction(1, 2)))

    def test_unrestricted_two_output_word_is_not_a_valid_direction(self):
        self.assertFalse(valid_pauli_word("IXIX"))
        self.assertFalse(valid_pauli_word("IXII"))


if __name__ == "__main__":
    unittest.main()
