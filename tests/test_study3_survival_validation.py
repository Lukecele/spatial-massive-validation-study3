#!/usr/bin/env python3
"""
Unit Test Suite for Study 3: Massive Population Survival Validation
Author: Luca Celebrano
Tests the validity of the filtered cohort and the statistical significance
of Kaplan-Meier and Cox hazard ratio estimators across N = 218 EGFR mutants.
"""

import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS_JSON = ROOT / "results" / "study3_survival_validation_results.json"
DATA_CSV = ROOT / "data" / "combined_filtered_massive_cohort.csv"
FIGURE_PNG = ROOT / "results" / "figure_study3_kaplan_meier_validation.png"


class TestStudy3SurvivalValidation(unittest.TestCase):
    def test_results_and_figures_exist(self):
        self.assertTrue(DATA_CSV.exists(), f"Missing {DATA_CSV}")
        self.assertTrue(RESULTS_JSON.exists(), f"Missing {RESULTS_JSON}")
        self.assertTrue(FIGURE_PNG.exists(), f"Missing {FIGURE_PNG}")

    def test_cohort_numbers_and_power(self):
        data = json.loads(RESULTS_JSON.read_text())
        self.assertEqual(data["n_total_patients"], 854)
        self.assertEqual(data["n_egfr_mutant"], 218)
        self.assertEqual(data["n_controls"], 636)
        self.assertEqual(data["statistical_power"], "> 97%")

    def test_statistical_significance_achieved(self):
        data = json.loads(RESULTS_JSON.read_text())
        pfi = data["pfi_analysis"]
        os = data["overall_survival_analysis"]

        # PFS / PFI must be highly significant (p < 0.001)
        self.assertTrue(pfi["statistically_significant"])
        self.assertLess(pfi["log_rank_p_value"], 0.001)
        self.assertGreater(pfi["hazard_ratio"], 2.0)

        # OS must be significant (p < 0.05)
        self.assertTrue(os["statistically_significant"])
        self.assertLess(os["log_rank_p_value"], 0.05)
        self.assertGreater(os["hazard_ratio"], 1.5)


if __name__ == "__main__":
    unittest.main()
