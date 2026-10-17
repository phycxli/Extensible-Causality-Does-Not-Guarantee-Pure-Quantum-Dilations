"""Independent exact physicality and witness audit of the EC counterexample."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
I2 = sp.eye(2)
X = sp.Matrix([[0, 1], [1, 0]])
Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
Z = sp.diag(1, -1)


def word(letters):
    lookup = {"I": I2, "X": X, "Y": Y, "Z": Z}
    return sp.kronecker_product(*(lookup[letter] for letter in letters))


def bits(index):
    return [(index >> (3 - position)) & 1 for position in range(4)]


def integer_index(values):
    return sum(value << (3 - position) for position, value in enumerate(values))


def replace_trace(matrix, positions):
    positions = tuple(positions)

    def entry(row, column):
        first, second = bits(row), bits(column)
        if any(first[position] != second[position] for position in positions):
            return sp.Integer(0)
        result = sp.Integer(0)
        for traced in itertools.product(range(2), repeat=len(positions)):
            left, right = list(first), list(second)
            for position, value in zip(positions, traced):
                left[position] = right[position] = value
            result += matrix[integer_index(left), integer_index(right)]
        return result / 2 ** len(positions)

    return sp.ImmutableMatrix(16, 16, entry)


def audit(output):
    identity = sp.eye(16)
    forward = word("IXXI") + word("IYYI") + word("IZZI")
    process = (identity + word("ZIII") / 5 + sp.Rational(3, 10) * word("ZIXZ") + forward / 5) / 4
    if process != process.adjoint() or sp.trace(process) != 4:
        raise AssertionError("process Hermiticity or normalization failed")
    spectrum = process.eigenvals()
    if min(spectrum) != sp.Rational(1, 40) or sum(spectrum.values()) != 16:
        raise AssertionError("exact full-rank positive spectrum failed")
    valid_checks = [
        replace_trace(process, (2, 3)) == replace_trace(process, (1, 2, 3)),
        replace_trace(process, (0, 1)) == replace_trace(process, (0, 1, 3)),
        process == replace_trace(process, (1,)) + replace_trace(process, (3,)) - replace_trace(process, (1, 3)),
    ]
    if not all(valid_checks):
        raise AssertionError("process linear-subspace constraint failed")
    singlet = sp.Matrix([0, 1, -1, 0])
    singlet_projector = singlet * singlet.T / 2
    witness = (sp.kronecker_product(I2 - Z, singlet_projector, I2)
               + (word("IIXZ") - word("ZIXZ")) / 4)
    ordered_replacements = {"A_before_B": replace_trace(witness, (3,)),
                            "B_before_A": replace_trace(witness, (1,))}
    for matrix in ordered_replacements.values():
        if min(matrix.eigenvals()) < 0:
            raise AssertionError("witness is not positive on an ordered-process cone")
    expectation = sp.trace(witness * process)
    noise_expectation = sp.trace(witness * identity / 4)
    if expectation != -sp.Rational(1, 10) or noise_expectation != 1:
        raise AssertionError("witness expectations failed")
    visibility = sp.Symbol("v", real=True)
    noisy = visibility * process + (1 - visibility) * identity / 4
    first = (3 * visibility * identity / 5 + visibility * forward / 5) / 4
    second = ((1 - 3 * visibility / 5) * identity + visibility * word("ZIII") / 5
              + 3 * visibility * word("ZIXZ") / 10) / 4
    if (first + second - noisy).applyfunc(sp.expand) != sp.zeros(16):
        raise AssertionError("white-noise CS decomposition identity failed")
    if first != replace_trace(first, (3,)) or second != replace_trace(second, (1,)):
        raise AssertionError("white-noise decomposition has the wrong causal orders")
    if min(first.subs(visibility, 1).eigenvals()) < 0:
        raise AssertionError("forward branch is not PSD for nonnegative visibility")
    smallest_second = (1 - sp.Rational(11, 10) * visibility) / 4
    for s, t in itertools.product((1, -1), repeat=2):
        sector = (identity + s * word("ZIII")) * (identity + t * word("IIXZ")) / 4
        eigenvalue = (1 - 3 * visibility / 5 + s * visibility / 5 + s * t * 3 * visibility / 10) / 4
        if (second * sector - eigenvalue * sector).applyfunc(sp.expand) != sp.zeros(16):
            raise AssertionError("symbolic reverse-order branch spectrum failed")
    noisy_witness = sp.factor(sp.trace(witness * noisy))
    if sp.expand(noisy_witness - (1 - 11 * visibility / 10)) != 0:
        raise AssertionError("exact noise witness threshold failed")
    payload = {
        "schema": "independent-input-biased-process-and-CS-witness-audit-v1",
        "utc": datetime.now(timezone.utc).isoformat(), "accepted": True,
        "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "project_science_module_imported": False,
        "process": "[I+(1/5)ZIII+(3/10)ZIXZ+(1/5)(IXXI+IYYI+IZZI)]/4",
        "port_order": ["AI", "AO", "BI", "BO"], "trace": "4", "rank": 16,
        "linear_process_constraints_exact": valid_checks, "minimum_eigenvalue": "1/40",
        "spectrum": [{"value": str(value), "multiplicity": int(spectrum[value])} for value in sorted(spectrum)],
        "witness": "(I-Z_AI) tensor P_singlet_AO_BI tensor I_BO + (IIXZ-ZIXZ)/4",
        "output_replaced_witnesses_PSD_exact": True,
        "witness_expectation": str(expectation), "white_noise_witness_expectation": str(noise_expectation),
        "noise_parameter": "W_v=v W_star+(1-v)I/4, 0<=v<=1",
        "noise_witness": str(noisy_witness), "noise_CS_threshold_visibility": "10/11",
        "white_noise_tolerance_nonpurifiable": "0<=noise_fraction<1/11",
        "explicit_noise_CS_decomposition_exact": True,
        "reverse_order_branch_minimum_eigenvalue": str(smallest_second),
        "Purifiable_equals_CS_reference": "Yokojima et al., Quantum 5, 441 (2021), Theorem 4 and Corollary 1, trivial original global past and future; ordered-comb Stinespring construction for the converse.",
        "scope": "Exact physicality, non-CS witness and complete white-noise CS classification. EC uses the separately verified universal Hall SOS and transportation proof; nonpurifiability also uses the cited pure-process structure theorem."
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii", newline="\n")
    print(json.dumps({key: payload[key] for key in ("accepted", "minimum_eigenvalue", "witness_expectation", "noise_CS_threshold_visibility")}, indent=2))
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "results/independent_counterexample_process_audit.json")
    arguments = parser.parse_args()
    audit(arguments.output.resolve())
