"""Check the two laboratory response protocols from exact Choi probabilities."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
I2 = sp.eye(2)
PAULI = {
    "I": I2,
    "X": sp.Matrix([[0, 1], [1, 0]]),
    "Y": sp.Matrix([[0, -sp.I], [sp.I, 0]]),
    "Z": sp.diag(1, -1),
}


def tensor(*matrices):
    return sp.kronecker_product(*matrices)


def word(letters):
    return tensor(*(PAULI[letter] for letter in letters))


def projector(axis, value):
    return (I2 + value * PAULI[axis]) / 2


def choi_measure_prepare(effect, state):
    standard = sp.zeros(4)
    for row in range(2):
        for column in range(2):
            unit = sp.zeros(2)
            unit[row, column] = 1
            standard += tensor(unit, sp.trace(effect * unit) * state)
    return standard.T


def check_instrument(events):
    for event in events:
        if event != event.adjoint() or min(event.eigenvals()) < 0:
            raise ValueError("A local Choi event is not positive")
    total = sum(events, sp.zeros(4))
    reduced = sp.Matrix(2, 2, lambda row, column:
                        sum(total[2 * row + output, 2 * column + output] for output in range(2)))
    if reduced != I2:
        raise ValueError("A local instrument is not trace preserving")


def require_equal(actual, expected, description):
    if sp.simplify(actual - expected) != 0:
        raise ValueError("Operational Born check failed: " + description)


def audit(output):
    delta, alpha, beta = sp.Rational(1, 5), sp.Rational(3, 10), sp.Rational(1, 5)
    process = (word("IIII") + delta * word("ZIII") + alpha * word("ZIXZ")
               + beta * (word("IXXI") + word("IYYI") + word("IZZI"))) / 4
    forward_cases, reverse_cases = [], []
    for axis in ("X", "Y", "Z"):
        for setting in (-1, 1):
            alice = choi_measure_prepare(I2, projector(axis, setting))
            bob_events = [choi_measure_prepare(projector(axis, outcome), I2 / 2) for outcome in (-1, 1)]
            check_instrument([alice])
            check_instrument(bob_events)
            sign = -1 if axis == "Y" else 1
            probabilities = {}
            for outcome, bob in zip((-1, 1), bob_events):
                value = sp.trace(process * tensor(alice, bob))
                require_equal(value, (1 + sign * beta * setting * outcome) / 2,
                              f"forward {axis}, setting {setting}, outcome {outcome}")
                probabilities[outcome] = value
            require_equal(sum(probabilities.values()), 1, "forward normalization")
            mean = sum(value * probabilities[value] for value in (-1, 1))
            require_equal(mean, sign * beta * setting, "forward mean")
            forward_cases.append({"axis": axis, "setting": setting,
                                  "p_plus": str(probabilities[1]), "mean_outcome": str(mean)})

    for setting in (-1, 1):
        branches = {}
        alice_events = [choi_measure_prepare(projector("Z", value), I2 / 2) for value in (-1, 1)]
        bob_events = [choi_measure_prepare(projector("X", value), projector("Z", setting * value))
                      for value in (-1, 1)]
        check_instrument(alice_events)
        check_instrument(bob_events)
        for alice_outcome in (-1, 1):
            alice = choi_measure_prepare(projector("Z", alice_outcome), I2 / 2)
            for bob_outcome in (-1, 1):
                bob = choi_measure_prepare(projector("X", bob_outcome), projector("Z", setting * bob_outcome))
                value = sp.trace(process * tensor(alice, bob))
                require_equal(value, (1 + delta * alice_outcome + alpha * setting * alice_outcome) / 4,
                              "reverse measurement and preparation branch")
                branches[alice_outcome, bob_outcome] = value
        require_equal(sum(branches.values()), 1, "reverse normalization")
        marginal = {value: sum(branches[value, bob] for bob in (-1, 1)) for value in (-1, 1)}
        mean = sum(value * marginal[value] for value in (-1, 1))
        require_equal(mean, delta + alpha * setting, "reverse mean")
        reverse_cases.append({"setting": setting, "p_plus": str(marginal[1]), "mean_outcome": str(mean)})

    wrong_alice = tensor(I2, projector("Y", 1))
    correct_bob = choi_measure_prepare(projector("Y", 1), I2 / 2)
    wrong = sp.trace(process * tensor(wrong_alice, correct_bob))
    correct = (1 - beta) / 2
    if wrong == correct:
        raise ValueError("The complex-transpose negative control is uninformative")

    sources = ("main/main.tex", "supplement/supplement.tex")
    record = {
        "schema": "paper8-exact-operational-responses-v1",
        "utc": datetime.now(timezone.utc).isoformat(), "accepted": True,
        "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sources},
        "process_parameters": {"delta": str(delta), "alpha": str(alpha), "beta": str(beta)},
        "forward_X_mean": "x/5", "forward_Y_mean": "-x/5", "reverse_mean": "1/5+3*y/10",
        "forward_cases": forward_cases, "reverse_cases": reverse_cases,
        "all_probabilities_normalized_and_nonnegative": all(
            0 <= sp.Rational(item["p_plus"]) <= 1 for item in forward_cases + reverse_cases),
        "local_instruments_CP_and_TP": True,
        "negative_control": {"incorrect_single_Y_transpose_p_plus": str(wrong),
                             "correct_p_plus": str(correct), "rejected": wrong != correct},
        "scope": "Two illustrative local protocols on the given process; no hardware realization or causal-inequality violation is inferred.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record, indent=2) + "\n", encoding="ascii", newline="\n")
    print("EXACT_OPERATIONAL_RESPONSES_PASSED", len(forward_cases), "forward cases,", len(reverse_cases), "reverse cases")
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "results/independent_operational_response_audit.json")
    audit(parser.parse_args().output)
