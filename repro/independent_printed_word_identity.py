"""Expand the printed SOS identity in exact free words and export every coefficient."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ("supplement/supplement.tex", "supplement/hall_sos_matrices.tex",
          "supplement/hall_sos_ideal_multipliers.tex")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unit():
    return {((), ()): F(1)}


def add(*polynomials):
    out = defaultdict(F)
    for polynomial in polynomials:
        for word, coefficient in polynomial.items():
            out[word] += coefficient
    return {word: coefficient for word, coefficient in out.items() if coefficient}


def scale(number, polynomial):
    return {word: F(number) * coefficient for word, coefficient in polynomial.items() if number * coefficient}


def multiply(*polynomials):
    out = unit()
    for polynomial in polynomials:
        expanded = defaultdict(F)
        for (a, b), x in out.items():
            for (c, d), y in polynomial.items():
                expanded[a + c, b + d] += x * y
        out = {word: coefficient for word, coefficient in expanded.items() if coefficient}
    return out


def adjoint(polynomial):
    return {(tuple((i + 4) % 8 for i in reversed(a)), tuple((j, c, r) for j, r, c in reversed(b))): x
            for (a, b), x in polynomial.items()}


def amplitude(sector, output):
    return {((2 * sector + output,), ()): F(1)}


def projector(label, row, column):
    if label == 2:
        return add(scale(int(row == column), unit()), scale(-1, projector(0, row, column)),
                   scale(-1, projector(1, row, column)))
    return {((), ((label, row, column),)): F(1)}


def gram(words, localizer=None):
    localizer = unit() if localizer is None else localizer
    return [[multiply(adjoint(left), localizer, right) for right in words] for left in words]


def contract(matrix, entries):
    assert len(matrix) == len(entries) and all(len(row) == len(matrix) for row in matrix)
    return add(*(scale(matrix[i][j], entries[i][j])
                 for i, j in itertools.product(range(len(matrix)), repeat=2)))


def matrix(body):
    result = [[int(cell.strip()) for cell in row.strip().split("&")]
              for row in body.strip().split(r"\\")]
    assert result and all(len(row) == len(result) for row in result)
    return result


def equation(text, label):
    blocks = re.findall(r"\\begin\{(equation|align)\}(.*?)\\end\{\1\}", text, re.S)
    matches = [body for _, body in blocks if r"\label{" + label + "}" in body]
    assert len(matches) == 1, label
    return matches[0]


def numeric_fractions(text):
    pattern = r"\\(?:tfrac|frac)\s*(?:\{(-?\d+)\}|(-?\d))\s*(?:\{(\d+)\}|(\d))"
    return [F(int(a or b), int(c or d)) for a, b, c, d in re.findall(pattern, text)]


def printed_data(root):
    texts = {name: (root / name).read_text(encoding="utf-8-sig") for name in INPUTS}
    sm = texts[INPUTS[0]]
    physical = numeric_fractions(equation(sm, "eq:Wstar"))
    assert len(physical) == 4 and physical[0] == F(1, 4)
    _, delta, alpha, beta = physical
    a_values = numeric_fractions(equation(sm, "eq:kernel"))
    assert a_values == [(1 + delta + beta) / 4, (1 - delta + beta) / 4]
    margins = numeric_fractions(equation(sm, "eq:margins"))
    assert len(margins) == 5
    q, theta, w_plus, w_minus, c = margins
    assert (w_plus, w_minus, c) == ((q + theta) / 2, (q - theta) / 2, (1 - q) / 2)
    denominator_text = re.search(r"\\kappa=([0-9]+(?:\\,[0-9]+)*)", sm)[1]
    denominator = int(denominator_text.replace(r"\,", ""))
    assert denominator == 4_000_000
    assert "There is no square-root normalization factor." in sm
    hall = equation(sm, "eq:operator-H")
    assert r"a_sT+s\alpha V/4-w_sI" in hall
    assert r"-\frac\beta2\sum_{s,o,p}L_{s,o}^\dagger L_{s,p}F_{op}-cT" in hall
    positive = re.findall(r"N_\{(\d+)\}=\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}", texts[INPUTS[1]], re.S)
    assert [int(index) for index, _ in positive] == list(range(15))
    positive = [matrix(body) for _, body in positive]
    ideal_text = texts[INPUTS[2]]
    bare = tuple(int(re.search(symbol + r"_0=(-?\d+)", ideal_text)[1]) for symbol in ("a", "b", "c"))
    ideal_matrices = {(symbol, 0 if label == "+" else 1): matrix(body)
                      for symbol, label, body in re.findall(r"([UV])_\{([+-])\}=\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}", ideal_text, re.S)}
    assert set(ideal_matrices) == set(itertools.product(("U", "V"), range(2)))
    return positive, bare, ideal_matrices, {
        "kappa": denominator, "delta": delta, "alpha": alpha, "beta": beta,
        "a_plus": a_values[0], "a_minus": a_values[1], "q": q, "theta": theta,
        "w_plus": w_plus, "w_minus": w_minus, "c": c,
    }


def positive_operator_entries():
    complement = [projector(2, r, c) for r, c in itertools.product(range(2), repeat=2)]
    n = [add(complement[0], scale(-1, complement[3])), complement[1], complement[2]]
    e = [projector(j, r, c) for j, r, c in itertools.product(range(2), repeat=3)]
    result = [gram(n)]
    for sector in range(2):
        amplitudes = [amplitude(sector, o) for o in range(2)]
        slack = add(unit(), *(scale(-1, multiply(adjoint(a), a)) for a in amplitudes))
        v = [multiply(a, b) for a in amplitudes for b in complement]
        u = [[add(v[i], scale(sign, v[7 - i])) for i in range(4)] for sign in (1, -1)]
        z = [[add(e[i], scale(sign, e[j])) for i, j in zip((0, 1, 4, 5), (3, 2, 7, 6))]
             for sign in (1, -1)]
        result.extend(gram(words) for words in u)
        result.extend(gram(words, slack) for words in z)
        localizer = [[multiply(adjoint(amplitudes[o]), projector(2, r, c), amplitudes[p])
                      for p, c in itertools.product(range(2), repeat=2)]
                     for o, r in itertools.product(range(2), repeat=2)]
        groups = ((0,), (1, 2), (3,))
        result.append([[add(*(localizer[i][j] for i in left for j in right))
                        for right in groups] for left in groups])
        result.extend([[multiply(slack, projector(label, r, c)) for c in range(2)]
                       for r in range(2)] for label in range(2))
    assert len(result) == 15
    return result


def ideal_operator(bare, matrices):
    a, b, c = bare
    scalar = {(0, 0): ((a, b), (b, a)), (1, 1): ((a, -b), (-b, a)),
              (0, 1): ((c, 0), (0, c)), (1, 0): ((c, 0), (0, c))}
    out, count, terms = {}, 0, []
    for j, ell, r, column in itertools.product(range(2), repeat=4):
        relation = add(*(multiply(projector(j, r, t), projector(ell, t, column)) for t in range(2)),
                       scale(-int(j == ell), projector(j, r, column)))
        multiplier = scale(scalar[j, ell][r][column], unit())
        count += 1
        for sector, o, p in itertools.product(range(2), repeat=3):
            u, v = matrices["U", sector], matrices["V", sector]
            row_index, col_index = 2 * o + r, 2 * p + column
            signs = (1, -1, -1, 1)
            coefficient = (u[row_index][col_index] if (j, ell) == (0, 0) else
                           signs[row_index] * u[row_index][col_index] * signs[col_index] if (j, ell) == (1, 1) else
                           v[row_index][col_index] if (j, ell) == (0, 1) else v[col_index][row_index])
            multiplier = add(multiplier, scale(coefficient, multiply(adjoint(amplitude(sector, o)), amplitude(sector, p))))
            count += 1
        terms.append(multiply(multiplier, relation))
        out = add(out, terms[-1])
    assert count == 144
    return out, terms


def negative_hall(parameters):
    total = [[add(projector(0, r, c), projector(1, r, c)) for c in range(2)] for r in range(2)]
    t = add(total[0][0], total[1][1])
    v = add(projector(0, 0, 1), projector(0, 1, 0), scale(-1, projector(1, 0, 1)), scale(-1, projector(1, 1, 0)))
    f = [[total[1][1], scale(-1, total[1][0])], [scale(-1, total[0][1]), total[0][0]]]
    out = scale(parameters["c"], t)
    for sector, sign in enumerate((1, -1)):
        a = parameters["a_plus"] if sector == 0 else parameters["a_minus"]
        w = parameters["w_plus"] if sector == 0 else parameters["w_minus"]
        for o in range(2):
            norm = multiply(adjoint(amplitude(sector, o)), amplitude(sector, o))
            out = add(out, multiply(norm, add(scale(w, unit()), scale(-a, t), scale(-sign * parameters["alpha"] / 4, v))))
        for o, p in itertools.product(range(2), repeat=2):
            out = add(out, scale(parameters["beta"] / 2,
                      multiply(adjoint(amplitude(sector, o)), amplitude(sector, p), f[o][p])))
    return out


def word_tex(word):
    a, b = word
    factors = [rf"\ell_{{{i % 4}}}" + (r"^\dagger" if i >= 4 else "") for i in a]
    factors.extend(rf"p_{{{'+' if j == 0 else '-'},{r}{c}}}" for j, r, c in b)
    return " ".join(factors) or "1"


def coefficient_row(word, target, q, ideal):
    values = (target.get(word, F()), q.get(word, F()), ideal.get(word, F()))
    residual = values[1] + values[2] - values[0]
    assert all(value.denominator == 1 for value in (*values, residual))
    return {"Alice_word": list(word[0]), "Bob_word": [list(item) for item in word[1]],
            "word_tex": word_tex(word), "target": str(values[0]), "Q": str(values[1]),
            "ideal": str(values[2]), "residual": str(residual), "in_target_support": word in target}



def complete_word_basis():
    alice = [()] + [(4 + 2 * f + o, 2 * f + p)
                    for f, o, p in itertools.product(range(2), repeat=3)]
    letters = list(itertools.product(range(2), repeat=3))
    bob = [()] + [(b,) for b in letters] + [(b, c) for b in letters for c in letters]
    basis = {(a, b) for a in alice for b in bob}
    assert len(alice) == 9 and len(bob) == 73 and len(basis) == 657
    return basis


def reconstruct(root):
    positive, bare, matrices, parameters = printed_data(root)
    entries = positive_operator_entries()
    parts = [contract(weight, entry) for weight, entry in zip(positive, entries)]
    weighted_entries = [scale(weight[r][c], entry[r][c])
                        for weight, entry in zip(positive, entries)
                        for r, c in itertools.product(range(len(weight)), repeat=2) if weight[r][c]]
    q = add(*parts)
    ideal, ideal_parts = ideal_operator(bare, matrices)
    target = scale(parameters["kappa"], negative_hall(parameters))
    generated = set().union(*(set(part) for part in [*weighted_entries, *ideal_parts]))
    aggregate = set(q) | set(ideal) | set(target)
    complete = complete_word_basis()
    assert generated | aggregate <= complete
    residual = add(q, ideal, scale(-1, target))
    assert not residual and len(target) == 48
    internally_cancelled = generated - aggregate
    absent = complete - generated
    ordered = sorted(target) + sorted(aggregate - set(target)) + sorted(internally_cancelled) + sorted(absent)
    assert len(ordered) == len(set(ordered)) == len(complete)
    rows = []
    for word in ordered:
        row = coefficient_row(word, target, q, ideal)
        row["generated_in_expansion"] = word in generated
        row["Q_block_coefficients"] = {f"N_{i}": str(part[word]) for i, part in enumerate(parts) if part.get(word)}
        assert sum(F(value) for value in row["Q_block_coefficients"].values()) == q.get(word, F())
        rows.append(row)
    negative_residuals = {
        "missing_parity_half": add(*(scale(2 if i in {1, 2, 3, 4, 8, 9, 10, 11} else 1, part)
                                       for i, part in enumerate(parts)), ideal, scale(-1, target)),
        "denominator_plus_one": add(q, ideal, scale(-(parameters["kappa"] + 1), negative_hall(parameters))),
        "reverse_ideal_sign": add(q, scale(-1, ideal), scale(-1, target)),
    }
    outside_word = sorted(aggregate - set(target))[0]
    assert ideal.get(outside_word, F()) != 0
    damaged_ideal = {word: value for word, value in ideal.items() if word != outside_word}
    negative_residuals["omit_one_external_ideal_coefficient"] = add(q, damaged_ideal, scale(-1, target))
    cancelled_word = sorted(internally_cancelled)[0]
    cancelled_part = next(part[cancelled_word] for part in parts if part.get(cancelled_word))
    negative_residuals["omit_one_internally_cancelled_Q_coefficient"] = add(
        q, scale(-cancelled_part, {cancelled_word: F(1)}), ideal, scale(-1, target))
    controls = {}
    for name, bad in negative_residuals.items():
        assert bad
        controls[name] = {"rejected": True, "nonzero_residual_words": len(bad),
                          "residuals_on_target_support": len(set(bad) & set(target)),
                          "residuals_outside_target_support": len(set(bad) - set(target))}
    for name in ("omit_one_external_ideal_coefficient", "omit_one_internally_cancelled_Q_coefficient"):
        assert controls[name]["residuals_on_target_support"] == 0
        assert controls[name]["residuals_outside_target_support"] == 1
    outside = aggregate - set(target)
    assert (len(generated), len(aggregate), len(outside), len(internally_cancelled), len(absent)) == (448, 204, 156, 244, 209)
    report = {
        "parameters": {key: str(value) for key, value in parameters.items()},
        "word_convention": {
            "Alice_generators": ["L_plus_0", "L_plus_1", "L_minus_0", "L_minus_1",
                                 "L_plus_0_dagger", "L_plus_1_dagger", "L_minus_0_dagger", "L_minus_1_dagger"],
            "Bob_letter": "[j,r,c] means the (r,c) block of P_plus for j=0 and P_minus for j=1",
            "canonical_order": "All Alice letters precede Bob letters; order within each party is retained",
            "algebra_before_comparison": "Only cross-party commutation and the adjoint rule; no projector relation reduced",
            "complement": "C_rc = delta_rc I - (P_plus)_rc - (P_minus)_rc",
            "denominator": "N_k and every ideal multiplier are unscaled integers; Target is the coefficient of -kappa H",
            "parity_basis": "u_even=v_rep+v_partner, u_odd=v_rep-v_partner; no square-root factor",
            "residual": "Q + Ideal - Target"},
        "positive_block_count": len(parts), "ideal_multiplier_count": 144,
        "complete_basis_words": len(complete), "generated_expansion_words": len(generated),
        "generated_words_outside_target_support": len(generated - set(target)),
        "target_support_words": len(target), "Q_support_words": len(q), "ideal_support_words": len(ideal),
        "aggregate_union_support_words": len(aggregate), "aggregate_words_outside_target_support": len(outside),
        "internally_cancelled_words": len(internally_cancelled), "ungenerated_basis_words": len(absent),
        "all_zero_unprinted_basis_words": len(complete - aggregate),
        "nonzero_residual_words_on_complete_basis": len(residual),
        "all_external_words_cancel": all(q.get(word, F()) + ideal.get(word, F()) == 0 for word in complete - set(target)),
        "all_complete_basis_coefficients": rows, "negative_controls": controls,
    }
    return report, [row for row in rows if row["word_tex"] in {word_tex(word) for word in aggregate}]


def inspect_printed_table(table_path, expected_rows):
    text = table_path.read_text(encoding="utf-8-sig")
    pattern = re.compile(r"^\s*\$(.*?)\$\s*&\s*(-?\d+)\s*&\s*(-?\d+)\s*&\s*(-?\d+)\s*&\s*(-?\d+)\s*\\\\\s*$")
    data_lines = [line for line in text.splitlines() if line.lstrip().startswith("$") and "&" in line]
    parsed = []
    for line in data_lines:
        match = pattern.fullmatch(line)
        assert match is not None, "A printed coefficient row is malformed"
        word, *values = match.groups()
        parsed.append({"word_tex": word, **dict(zip(("target", "Q", "ideal", "residual"), values))})
    assert len(parsed) == len(expected_rows) == 204, "Printed coefficient support is incomplete"
    normalize = lambda value: re.sub(r"\s+", "", value)
    observed = {}
    for row in parsed:
        word = normalize(row["word_tex"])
        assert word not in observed, "A printed free word occurs more than once"
        observed[word] = tuple(F(row[key]) for key in ("target", "Q", "ideal", "residual"))
    expected = {normalize(row["word_tex"]): tuple(F(row[key]) for key in ("target", "Q", "ideal", "residual"))
                for row in expected_rows}
    assert observed == expected, "A printed word or exact coefficient differs from the independent expansion"
    assert r"\mathcal W" in text and r"\mathcal A" in text and r"\mathcal B" in text
    assert all(str(value) in text for value in (657, 448, 204, 48, 156, 244, 209, 453))
    return {"all_rows_match_exactly": True, "printed_rows_checked": len(parsed),
            "word_duplicates": 0, "unlisted_coefficients_zero_on_complete_basis": 453}


def audit(root=ROOT, output=None):
    root = Path(root).resolve()
    source_names = (*INPUTS, "supplement/hall_sos_coefficient_check.tex",
                    "repro/independent_printed_word_identity.py")
    source_hashes = {name: sha(root / name) for name in source_names}
    assert source_hashes["repro/independent_printed_word_identity.py"] == sha(Path(__file__)), "Invoked verifier differs from the supplied public source"
    mathematical, printed_rows = reconstruct(root)
    printed = inspect_printed_table(root / "supplement/hall_sos_coefficient_check.tex", printed_rows)
    assert source_hashes == {name: sha(root / name) for name in source_names}, "Inputs changed during the audit"
    report = {"schema": "paper8-independent-printed-word-identity-v1",
              "utc": datetime.now(timezone.utc).isoformat(), "accepted": True,
              "source_sha256": source_hashes, "arithmetic": "fractions.Fraction",
              "project_science_or_generator_imported": False,
              "source_of_coefficients": "Printed TeX integer matrices and printed rational process parameters",
              "printed_table": printed, **mathematical}
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + "\n", encoding="ascii", newline="\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(root=args.root, output=args.output)
    print("PRINTED_FULL_WORD_IDENTITY_VERIFIED", result["complete_basis_words"], "basis;",
          result["generated_expansion_words"], "generated;", result["aggregate_union_support_words"],
          "printed; residual", result["nonzero_residual_words_on_complete_basis"])


if __name__ == "__main__":
    main()
