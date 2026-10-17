"""Exact checks for the all-direction neighborhood; no optimizer imports."""

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
PAULI = {
    "I": sp.eye(2), "X": sp.Matrix([[0, 1], [1, 0]]),
    "Y": sp.Matrix([[0, -sp.I], [sp.I, 0]]), "Z": sp.diag(1, -1),
}


def word(letters):
    return sp.kronecker_product(*(PAULI[letter] for letter in letters))


def replacement_coefficient(letters, ports):
    return int(all(letters[port] == "I" for port in ports))


def valid_pauli_word(letters):
    if len(letters) != 4 or any(letter not in PAULI for letter in letters):
        raise ValueError("expected four Pauli letters in AI,AO,BI,BO order")
    replace = lambda ports: replacement_coefficient(letters, ports)
    return (
        replace((2, 3)) == replace((1, 2, 3))
        and replace((0, 1)) == replace((0, 1, 3))
        and 1 == replace((1,)) + replace((3,)) - replace((1, 3))
    )


def ordered_split_coefficients(letters):
    a = Fraction(replacement_coefficient(letters, (3,)))
    b = Fraction(replacement_coefficient(letters, (1,)))
    shared = Fraction(replacement_coefficient(letters, (1, 3)), 2)
    return a - shared, b - shared


def certify_ball(noise_fraction, radius):
    eta, radius = Fraction(noise_fraction), Fraction(radius)
    if not 0 < eta < Fraction(1, 11) or radius <= 0:
        raise ValueError("a positive ball requires 0 < eta < 1/11")
    # I/8 minus the sum of two contracted perturbations bounds each branch.
    branch_floor = Fraction(1, 8) - Fraction(3, 2) * radius / eta
    witness_upper = (11 * eta - 1) / 10 + 8 * radius
    if branch_floor < 0:
        raise ValueError("radius exceeds the sufficient EC bound")
    if witness_upper >= 0:
        raise ValueError("radius has no certified strict non-CS witness margin")
    return {
        "noise_fraction": str(eta), "operator_norm_radius": str(radius),
        "CS_component_minimum_eigenvalue_bound": str(branch_floor),
        "center_witness": str((11 * eta - 1) / 10),
        "perturbed_witness_upper_bound": str(witness_upper),
    }


def audit(output):
    valid = []
    for letters in itertools.product(PAULI, repeat=4):
        if not valid_pauli_word(letters):
            continue
        valid.append("".join(letters))
        a, b = ordered_split_coefficients(letters)
        if a + b != 1:
            raise AssertionError("ordered split does not reconstruct a valid basis word")
        if a and letters[3] != "I" or b and letters[1] != "I":
            raise AssertionError("ordered split violates an output replacement constraint")
    if len(valid) != 88 or "IIII" not in valid:
        raise AssertionError("complete process linear subspace count failed")

    identity = sp.eye(16)
    process = (
        identity + word("ZIII") / 5 + 3 * word("ZIXZ") / 10
        + (word("IXXI") + word("IYYI") + word("IZZI")) / 5
    ) / 4
    singlet = (identity - word("IXXI") - word("IYYI") - word("IZZI")) / 4
    positive_part = (identity - word("ZIII")) * singlet
    pauli_part = (word("IIXZ") - word("ZIXZ")) / 4
    witness = positive_part + pauli_part
    if min(positive_part.eigenvals()) < 0 or sp.trace(positive_part) != 4:
        raise AssertionError("positive witness part is not PSD of trace four")
    if sum(abs(value) * count for value, count in pauli_part.eigenvals().items()) != 4:
        raise AssertionError("Pauli witness part trace norm failed")
    trace_norm = sp.simplify(sum(
        abs(value) * count for value, count in witness.eigenvals().items()
    ))
    if trace_norm != 2 + 2 * sp.sqrt(5) or not trace_norm < 8:
        raise AssertionError("witness trace norm failed")
    center = sp.Rational(21, 22) * process + identity / 88
    if sp.trace(witness * center) != -sp.Rational(1, 20):
        raise AssertionError("center witness failed")
    constants = certify_ball("1/22", "1/264")
    if constants["perturbed_witness_upper_bound"] != "-13/660":
        raise AssertionError("strict neighborhood witness bound failed")
    if Fraction(1, 8) - Fraction(3, 2) * Fraction(1, 12) != 0:
        raise AssertionError("CS ball positivity radius failed")
    pauli_sum = word("IXXI") + word("IYYI") + word("IZZI")
    if max(abs(value) for value in pauli_sum.eigenvals()) != 3:
        raise AssertionError("parameter perturbation norm bound failed")

    sources = ("main/main.tex", "supplement/supplement.tex",
               "results/exact_double_flag_hall_sos.json")
    record = {
        "schema": "independent-all-direction-EC-neighborhood-v1",
        "utc": datetime.now(timezone.utc).isoformat(), "accepted": True,
        "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": {
            name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sources
        },
        "project_science_module_imported": False,
        "valid_linear_space_dimension": len(valid), "trace_fixed_affine_dimension": len(valid) - 1,
        "ordered_split_verified_on_all_valid_Pauli_basis_words": True,
        "white_noise_CS_operator_norm_radius": "1/12",
        "witness_trace_norm_exact": str(trace_norm), "witness_trace_norm_bound": "8",
        "center_parameters": {"delta": "21/110", "alpha": "63/220", "beta": "21/110"},
        "parameter_deviation_weighted_L1_radius": "1/66",
        **constants,
        "proof_dependencies": [
            "Unital positive output replacements contract the Hermitian operator norm.",
            "CS is contained in EC, and EC is convex.",
            "The existing arbitrary-dimensional EC theorem for W_star.",
            "The externally cited bipartite Purifiable=CS theorem.",
        ],
        "scope": (
            "Exact complete-basis identities and neighborhood constants. The analytic norm "
            "argument and physical universal-EC proof are in the supplement; this audit "
            "is not external expert review and does not certify larger radii."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record, indent=2) + "\n", encoding="ascii", newline="\n")
    print(json.dumps({key: record[key] for key in (
        "accepted", "operator_norm_radius", "perturbed_witness_upper_bound",
        "trace_fixed_affine_dimension",
    )}))
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "results/independent_ec_robustness_audit.json")
    args = parser.parse_args()
    audit(args.output)
