"""Independently verify a rational, arbitrary-dimensional Hall SOS.

The verifier imports no solver or project science module. Words commute only
between the two laboratories. Every coefficient and Schur pivot is rational.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ADJOINT = (4, 5, 6, 7, 0, 1, 2, 3, 8, 10, 9, 11, 12, 14, 13, 15)


def ordered(word):
    return tuple(value for value in word if value < 8) + tuple(value for value in word if value >= 8)


def cleaned(poly):
    return {word: value for word, value in poly.items() if value}


def plus(*polynomials):
    answer = defaultdict(Fraction)
    for poly in polynomials:
        for word, value in poly.items():
            answer[ordered(word)] += value
    return cleaned(answer)


def scale(value, poly):
    return cleaned({word: value * coefficient for word, coefficient in poly.items()})


def times(*polynomials):
    answer = {(): Fraction(1)}
    for poly in polynomials:
        product = defaultdict(Fraction)
        for left, a in answer.items():
            for right, b in poly.items():
                product[ordered(left + right)] += a * b
        answer = cleaned(product)
    return answer


def adjoint(poly):
    return {ordered(tuple(ADJOINT[value] for value in reversed(word))): coefficient
            for word, coefficient in poly.items()}


def atom(*word):
    return {tuple(word): Fraction(1)}


def bob(label, row, column):
    if label in (0, 1):
        return atom(8 + 4 * label + 2 * row + column)
    return plus(scale(Fraction(row == column), atom()),
                scale(-1, bob(0, row, column)), scale(-1, bob(1, row, column)))


def alice_slack(flag):
    return plus(atom(), *(scale(-1, times(atom(4 + 2 * flag + output), atom(2 * flag + output)))
                          for output in range(2)))


def independent_minus_hall(q, theta):
    beta, alpha, delta = Fraction(1, 5), Fraction(3, 10), Fraction(1, 5)
    total = [[plus(bob(0, row, column), bob(1, row, column)) for column in range(2)]
             for row in range(2)]
    trace = plus(total[0][0], total[1][1])
    feedback = plus(bob(0, 0, 1), bob(0, 1, 0),
                    scale(-1, bob(1, 0, 1)), scale(-1, bob(1, 1, 0)))
    # Multiplication by Pauli Y on both sides gives these transposed off-diagonals.
    flipped = [[total[1][1], scale(-1, total[1][0])],
               [scale(-1, total[0][1]), total[0][0]]]
    terms = [scale((1 - q) / 2, trace)]
    for flag, sign in enumerate((1, -1)):
        w, a = (q + sign * theta) / 2, (1 + sign * delta + beta) / 4
        for output in range(2):
            norm = times(atom(4 + 2 * flag + output), atom(2 * flag + output))
            terms += [scale(w, norm), scale(-a, times(norm, trace)),
                      scale(-sign * alpha / 4, times(norm, feedback))]
        for row, column in itertools.product(range(2), repeat=2):
            norm = times(atom(4 + 2 * flag + row), atom(2 * flag + column))
            terms.append(scale(beta / 2, times(norm, flipped[row][column])))
    return plus(*terms)


def gram_entries(words, localizer=None):
    localizer = atom() if localizer is None else localizer
    return [[times(adjoint(left), localizer, right) for right in words] for left in words]


def independent_blocks():
    blocks = {}
    complement = [bob(2, row, column) for row, column in itertools.product(range(2), repeat=2)]
    traceless = [plus(complement[0], scale(-1, complement[3])), complement[1], complement[2]]
    blocks["scalar_Bob_traceless_complement_Gram"] = gram_entries(traceless)
    for flag in range(2):
        forward = [times(atom(2 * flag + output), word)
                   for output in range(2) for word in complement]
        blocks[f"flag_{flag}_forward_complement_Gram"] = gram_entries(forward)
        blocks[f"flag_{flag}_Alice_slack_Bob_Gram"] = gram_entries(
            [atom(value) for value in range(8, 16)], alice_slack(flag))
        full = [[times(atom(4 + 2 * flag + a), atom(2 * flag + b), bob(2, i, j))
                 for b, j in itertools.product(range(2), repeat=2)]
                for a, i in itertools.product(range(2), repeat=2)]
        triplet = ((1, 0, 0), (0, 1, 0), (0, 1, 0), (0, 0, 1))
        entries = [[plus(*(scale(triplet[a][i] * triplet[b][j], full[a][b])
                           for a, b in itertools.product(range(4), repeat=2)))
                    for j in range(3)] for i in range(3)]
        blocks[f"flag_{flag}_Bob_complement_localizer"] = entries
        for label in range(2):
            blocks[f"flag_{flag}_Bob_label_{label}_slack_product"] = [
                [times(alice_slack(flag), bob(label, row, column)) for column in range(2)]
                for row in range(2)]
    return blocks


def independent_ideal_terms():
    relations = []
    for left, right, row, column in itertools.product(range(2), repeat=4):
        relations.append(plus(*(times(bob(left, row, inner), bob(right, inner, column))
                                for inner in range(2)),
                              scale(-int(left == right), bob(left, row, column))))
    result = list(relations)
    for flag, row, column in itertools.product(range(2), repeat=3):
        prefix = times(atom(4 + 2 * flag + row), atom(2 * flag + column))
        result += [times(prefix, relation) for relation in relations]
    return result


def exact_ldl(matrix):
    side = len(matrix)
    if any(len(row) != side for row in matrix):
        raise ValueError("non-square SOS matrix")
    if any(matrix[row][column] != matrix[column][row]
           for row, column in itertools.product(range(side), repeat=2)):
        raise ValueError("a rational SOS matrix is not symmetric")
    lower = [[Fraction(int(row == column)) for column in range(side)] for row in range(side)]
    pivots = []
    for column in range(side):
        pivot = matrix[column][column] - sum(lower[column][k] ** 2 * pivots[k] for k in range(column))
        if pivot < 0:
            raise ValueError(f"negative exact Schur pivot at column {column}: {pivot}")
        pivots.append(pivot)
        for row in range(column + 1, side):
            residual = matrix[row][column] - sum(lower[row][k] * lower[column][k] * pivots[k]
                                                 for k in range(column))
            if pivot:
                lower[row][column] = residual / pivot
            elif residual:
                raise ValueError("zero Schur pivot has a nonzero column")
    for row, column in itertools.product(range(side), repeat=2):
        reconstructed = sum(lower[row][k] * pivots[k] * lower[column][k] for k in range(side))
        if reconstructed != matrix[row][column]:
            raise AssertionError("exact LDL reconstruction failed")
    return {"side": side, "positive_semidefinite": True,
            "positive_definite": all(value > 0 for value in pivots),
            "Schur_pivots": [str(value) for value in pivots]}


def audit(certificate_path, output):
    certificate = json.loads(certificate_path.read_text(encoding="utf-8"))
    q = Fraction(certificate.get("q", certificate.get("q_exact")))
    theta = Fraction(certificate.get("theta", certificate.get("theta_exact")))
    if (q, theta) != (Fraction(49, 80), Fraction(3, 25)):
        raise ValueError("this independent certificate audit fixes the reviewed rational margins")
    entries = independent_blocks()
    saved_blocks = certificate.get("blocks", certificate.get("PSD_blocks"))
    if len(saved_blocks) != 11 or set(block["name"] for block in saved_blocks) != set(entries):
        raise ValueError("SOS block inventory mismatch")
    accumulated, psd_reports = {}, []
    for block in saved_blocks:
        polynomial_matrix = entries[block["name"]]
        weight = [[Fraction(value) for value in row] for row in block["weight"]]
        if len(weight) != len(polynomial_matrix):
            raise ValueError(f"SOS block dimension mismatch: {block['name']}")
        psd_report = exact_ldl(weight)
        psd_report["name"] = block["name"]
        psd_reports.append(psd_report)
        accumulated = plus(accumulated, *(scale(weight[row][column], polynomial_matrix[row][column])
                                          for row, column in itertools.product(range(len(weight)), repeat=2)))
    multipliers = certificate["projector_multipliers"]
    relations = independent_ideal_terms()
    if len(multipliers) != len(relations) or len(relations) != 144:
        raise ValueError("projector ideal multiplier inventory mismatch")
    accumulated = plus(accumulated, *(scale(Fraction(value), relation)
                                      for value, relation in zip(multipliers, relations)))
    target = independent_minus_hall(q, theta)
    difference = plus(accumulated, scale(-1, target))
    if difference:
        first = next(iter(difference.items()))
        raise ValueError(f"exact SOS coefficient mismatch: {first}; count={len(difference)}")
    report = {
        "schema": "independent-rational-arbitrary-dimensional-Hall-SOS-audit-v1",
        "utc": datetime.now(timezone.utc).isoformat(), "accepted": True,
        "certificate": certificate_path.resolve().relative_to(ROOT).as_posix(),
        "certificate_sha256": hashlib.sha256(certificate_path.read_bytes()).hexdigest(),
        "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "q": str(q), "theta": str(theta), "process_parameters": {"delta": "1/5", "alpha": "3/10", "beta": "1/5"},
        "solver_imported": False, "project_science_module_imported": False,
        "coefficient_arithmetic": "fractions.Fraction", "independent_target_word_count": len(target),
        "coefficient_identity_exact": True, "coefficient_residual_nonzero_words": 0,
        "projector_ideal_term_count": len(relations), "PSD_blocks": psd_reports,
        "all_PSD_blocks_exactly_verified": True,
        "Hall_universal_bound_exact": "H<=0 after local Bob Naimark dilation",
        "scope": "Exact noncommutative polynomial identity and rational PSD verification. The independently written proof of physical coverage, transportation sufficiency, and mixed-state convexity is required to deduce EC."
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="ascii", newline="\n")
    print(json.dumps({"accepted": True, "exact_coefficient_identity": True,
                      "PSD_block_count": len(psd_reports), "certificate_sha256": report["certificate_sha256"]}, indent=2))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, default=ROOT / "results/reduced_hall_sos_exact.json")
    parser.add_argument("--output", type=Path, default=ROOT / "results/independent_general_ec_certificate_audit.json")
    args = parser.parse_args()
    audit(args.certificate.resolve(), args.output.resolve())
