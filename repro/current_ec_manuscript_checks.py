"""Check the current EC manuscript and its exact evidence without a generator."""

from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
FRAGMENTS = (
    "supplement/hall_sos_matrices.tex",
    "supplement/hall_sos_ideal_multipliers.tex",
    "supplement/hall_sos_positive_minors.tex",
)


def require_equal(actual, expected, label):
    if actual != expected:
        raise ValueError(f"current EC audit mismatch: {label}")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read_json(root, relative):
    return json.loads((root / relative).read_text(encoding="utf-8-sig"))


def check_hash(root, relative, recorded):
    require_equal(digest(root / relative), recorded, f"fresh source {relative}")


def matrix(values):
    rows = [[Fraction(value) for value in row] for row in values]
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("current EC audit requires a square matrix")
    return rows


def transpose(values):
    return [list(row) for row in zip(*values)]


def multiply(left, right):
    return [[sum((a * b for a, b in zip(row, column)), Fraction())
             for column in zip(*right)] for row in left]


def determinant(values):
    rows = [row[:] for row in matrix(values)]
    result = Fraction(1)
    for column in range(len(rows)):
        pivot = next((row for row in range(column, len(rows)) if rows[row][column]), None)
        if pivot is None:
            return Fraction()
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            result = -result
        value = rows[column][column]
        result *= value
        for row in range(column + 1, len(rows)):
            scale = rows[row][column] / value
            for index in range(column + 1, len(rows)):
                rows[row][index] -= scale * rows[column][index]
    return result


def tex_matrix(text, symbol):
    pattern = re.escape(symbol) + r"\s*=\s*\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}"
    matches = re.findall(pattern, text, re.S)
    require_equal(len(matches), 1, f"unique displayed {symbol}")
    return [[int(entry.strip()) for entry in row.split("&")]
            for row in matches[0].strip().split(r"\\") if row.strip()]


def check_compact_artifacts(root=ROOT):
    source_name = "results/exact_double_flag_hall_sos_symmetric.json"
    table_name = "results/hall_sos_compact_tables.json"
    source = read_json(root, source_name)
    tables = read_json(root, table_name)
    minors = read_json(root, "results/hall_sos_positive_principal_minors.json")
    ideal = read_json(root, "results/hall_sos_compact_ideal_multipliers.json")
    texts = [(root / relative).read_text(encoding="utf-8-sig") for relative in FRAGMENTS]
    require_equal(tables["schema"], "exact-Hall-SOS-compact-integer-tables-v1", "table schema")
    require_equal(minors["schema"], "exact-Hall-SOS-leading-principal-minors-v1", "minor schema")
    check_hash(root, source_name, tables["source_certificate_sha256"])
    check_hash(root, table_name, minors["source_table_sha256"])
    denominator = tables["common_denominator"]
    require_equal(denominator, 4000000, "common integer denominator")
    require_equal(minors["common_denominator"], denominator, "minor denominator")
    require_equal(ideal["common_denominator"], denominator, "ideal denominator")
    for payload in (source, tables):
        require_equal(Fraction(payload["q"]), Fraction(49, 80), "q")
        require_equal(Fraction(payload["theta"]), Fraction(3, 25), "theta")
    for key, value in (("delta", "1/5"), ("alpha", "3/10"), ("beta", "1/5")):
        require_equal(Fraction(source[key]), Fraction(value), f"SOS process parameter {key}")
    source_blocks = {block["name"]: matrix(block["weight"]) for block in source["blocks"]}
    require_equal(len(source_blocks), len(source["blocks"]), "unique source blocks")
    reconstructed = {name: [[Fraction() for _ in weight] for _ in weight]
                     for name, weight in source_blocks.items()}
    indices = [block["table_index"] for block in tables["blocks"]]
    expected_count = sum(2 if len(weight) == 8 else 1 for weight in source_blocks.values())
    require_equal(indices, list(range(expected_count)), "complete compact matrix indexing")
    expected_source_order = [name for name, weight in source_blocks.items()
                             for _ in range(2 if len(weight) == 8 else 1)]
    require_equal([block["source_block"] for block in tables["blocks"]], expected_source_order,
                  "compact indices match the manuscript Gram/localizer order")
    require_equal(len(minors["blocks"]), expected_count, "complete minor indexing")
    displayed = [int(index) for index in re.findall(r"N_\{(\d+)\}\s*=", texts[0])]
    require_equal(displayed, indices, "complete displayed matrices")
    determinant_rows = re.findall(r"^\s*(\d+)\s*&\s*(.*?)\\\\\s*$", texts[2], re.M)
    require_equal([int(index) for index, _ in determinant_rows], indices, "complete displayed minors")
    for block, minor, (_, minor_body) in zip(tables["blocks"], minors["blocks"], determinant_rows):
        index = block["table_index"]
        numerator = block["integer_numerator"]
        weight = matrix(block["weight"])
        side = len(weight)
        if side > 4:
            raise ValueError("compact SOS matrix exceeds the stated size four")
        require_equal(weight, [[Fraction(value, denominator) for value in row] for row in numerator],
                      f"N_{index}/denominator")
        require_equal(tex_matrix(texts[0], f"N_{{{index}}}"), numerator, f"displayed N_{index}")
        require_equal(weight, transpose(weight), f"symmetric weight {index}")
        determinants = [determinant([row[:size] for row in numerator[:size]])
                        for size in range(1, side + 1)]
        if any(value <= 0 for value in determinants):
            raise ValueError(f"compact SOS weight {index} is not positive definite")
        require_equal(minor["table_index"], index, f"minor index {index}")
        require_equal(list(map(Fraction, minor["integer_numerator_leading_determinants"])),
                      determinants, f"integer determinants {index}")
        previous = [Fraction(1)] + determinants[:-1]
        pivots = [value / before for value, before in zip(determinants, previous)]
        require_equal(list(map(Fraction, minor["integer_numerator_LDL_pivots"])), pivots,
                      f"integer LDL pivots {index}")
        weight_pivots = [value / denominator for value in pivots]
        require_equal(list(map(Fraction, minor["weight_LDL_pivots"])), weight_pivots,
                      f"weight LDL pivots {index}")
        require_equal(list(map(Fraction, block["LDL_positive_diagonal"])), weight_pivots,
                      f"compact LDL pivots {index}")
        lower = matrix(block["LDL_lower"])
        if any(lower[row][column] != int(row == column) for row in range(side)
               for column in range(row, side)):
            raise ValueError(f"non-unit lower LDL matrix {index}")
        scaled = [[value * weight_pivots[column] for column, value in enumerate(row)] for row in lower]
        require_equal(multiply(scaled, transpose(lower)), weight, f"compact LDL reconstruction {index}")
        cells = [cell.strip() for cell in minor_body.split("&")]
        expected_cells = [str(value) for value in determinants] + ["$-$"] * (4 - side)
        require_equal(cells, expected_cells, f"displayed principal minors {index}")
        original = source_blocks[block["source_block"]]
        original_side = len(original)
        if block["basis_transform"] == "identity":
            transform = [[Fraction(int(row == column)) for column in range(side)]
                         for row in range(original_side)]
        else:
            representatives, partners = block["representative_indices"], block["partner_indices"]
            require_equal(sorted(representatives + partners), list(range(original_side)),
                          f"complete symmetry basis {index}")
            if "forward_complement" in block["source_block"]:
                require_equal((representatives, partners), ([0, 1, 2, 3], [7, 6, 5, 4]),
                              f"declared forward symmetry basis {index}")
            elif "Alice_slack_Bob" in block["source_block"]:
                require_equal((representatives, partners), ([0, 1, 4, 5], [3, 2, 7, 6]),
                              f"declared slack symmetry basis {index}")
            else:
                raise ValueError(f"undeclared symmetry basis for compact block {index}")
            if block["basis_transform"] == "u_i=v_rep_i+v_partner_i":
                sign = 1
            elif block["basis_transform"] == "u_i=v_rep_i-v_partner_i":
                sign = -1
            else:
                raise ValueError(f"unknown symmetry transform for compact block {index}")
            transform = [[Fraction() for _ in range(side)] for _ in range(original_side)]
            for column, (representative, partner) in enumerate(zip(representatives, partners)):
                transform[representative][column] = 1
                transform[partner][column] = sign
        contribution = multiply(multiply(transform, weight), transpose(transform))
        accumulated = reconstructed[block["source_block"]]
        for row in range(original_side):
            for column in range(original_side):
                accumulated[row][column] += contribution[row][column]
    require_equal(reconstructed, source_blocks, "all compact weights recover the exact source SOS")

    a, b, c = (ideal[key] for key in ("bare_a", "bare_b", "bare_c"))
    for name, value in (("a", a), ("b", b), ("c", c)):
        found = re.findall(name + r"_0\s*=\s*(-?\d+)", texts[1])
        require_equal(found, [str(value)], f"displayed bare ideal {name}")
    ideal_blocks = {}
    for block in ideal["blocks"]:
        key = (block["flag"], block["kind"])
        if key in ideal_blocks:
            raise ValueError("duplicate compact ideal block")
        ideal_blocks[key] = block["integer_matrix"]
        require_equal(tex_matrix(texts[1], block["symbol"]), block["integer_matrix"],
                      f"displayed {block['symbol']}")
    require_equal(set(ideal_blocks), {(flag, kind) for flag in (0, 1)
                                     for kind in ("same_label", "cross_label")}, "four ideal blocks")
    eta = {(0, 0): [[a, b], [b, a]], (1, 1): [[a, -b], [-b, a]],
           (0, 1): [[c, 0], [0, c]], (1, 0): [[c, 0], [0, c]]}
    multipliers = [Fraction(eta[left, right][row][column], denominator)
                   for left in (0, 1) for right in (0, 1) for row in (0, 1) for column in (0, 1)]
    z = (1, -1, -1, 1)
    for flag in (0, 1):
        same, cross = ideal_blocks[flag, "same_label"], ideal_blocks[flag, "cross_label"]
        for alice_row in (0, 1):
            for alice_column in (0, 1):
                for left in (0, 1):
                    for right in (0, 1):
                        for bob_row in (0, 1):
                            for bob_column in (0, 1):
                                row, column = 2 * alice_row + bob_row, 2 * alice_column + bob_column
                                if left == right == 0:
                                    value = same[row][column]
                                elif left == right == 1:
                                    value = z[row] * same[row][column] * z[column]
                                elif left == 0:
                                    value = cross[row][column]
                                else:
                                    value = cross[column][row]
                                multipliers.append(Fraction(value, denominator))
    require_equal(multipliers, list(map(Fraction, source["projector_multipliers"])),
                  "all compact projection-ideal multipliers recover the exact SOS")
    return {"compact_block_count": expected_count, "common_denominator": denominator,
            "maximum_block_size": max(len(block["weight"]) for block in tables["blocks"]),
            "integer_matrix_transcriptions_exact": True, "independent_Sylvester_checks_exact": True,
            "compact_to_source_reconstruction_exact": True,
            "ideal_multiplier_count": len(multipliers), "ideal_multiplier_transcriptions_exact": True,
            "leading_minor_transcriptions_exact": True}


def check_evidence(root=ROOT):
    receipts = (
        ("results/independent_general_ec_certificate_audit.json", "repro/independent_general_ec_certificate_audit.py"),
        ("results/independent_symmetric_hall_sos_audit.json", "repro/independent_general_ec_certificate_audit.py"),
        ("results/independent_counterexample_process_audit.json", "repro/independent_counterexample_process_audit.py"),
        ("results/independent_ec_robustness_audit.json", "repro/independent_ec_robustness_audit.py"),
        ("results/independent_operational_response_audit.json", "repro/independent_operational_response_audit.py"),
        ("results/input_bias_counterexample_theorem_audit.json", "repro/audit_input_bias_counterexample_theorem.py"),
        ("results/literature_scope/independent_two_certificate_ec_audit.json",
         "results/literature_scope/independent_two_certificate_ec_audit.py"),
        ("results/literature_scope/independent_ec_physical_coverage_audit.json",
         "results/literature_scope/independent_ec_physical_coverage_audit.py"),
    )
    data = {}
    for receipt, verifier in receipts:
        report = read_json(root, receipt)
        require_equal(report["accepted"], True, f"accepted {receipt}")
        check_hash(root, verifier, report["verifier_sha256"])
        if "certificate" in report:
            check_hash(root, report["certificate"].replace("\\", "/"), report["certificate_sha256"])
        if receipt == "results/independent_ec_robustness_audit.json":
            for relative, checksum in report["source_sha256"].items():
                check_hash(root, relative, checksum)
            require_equal(report["operator_norm_radius"], "1/264", "neighborhood radius")
            require_equal(report["noise_fraction"], "1/22", "neighborhood center")
            require_equal(report["perturbed_witness_upper_bound"], "-13/660", "neighborhood witness margin")
            require_equal(report["trace_fixed_affine_dimension"], 87, "complete perturbation dimension")
        if receipt == "results/independent_operational_response_audit.json":
            for relative, checksum in report["source_sha256"].items():
                check_hash(root, relative, checksum)
            require_equal(report["forward_X_mean"], "x/5", "Alice-to-Bob response")
            require_equal(report["reverse_mean"], "1/5+3*y/10", "Bob-to-Alice response")
            require_equal(report["local_instruments_CP_and_TP"], True, "physical response instruments")
            require_equal(report["negative_control"]["rejected"], True, "complex transpose rejection")
        if "q" in report:
            require_equal(Fraction(report["q"]), Fraction(49, 80), f"q in {receipt}")
            require_equal(Fraction(report["theta"]), Fraction(3, 25), f"theta in {receipt}")
        if "process_parameters" in report:
            require_equal({key: Fraction(value) for key, value in report["process_parameters"].items()},
                          {"delta": Fraction(1, 5), "alpha": Fraction(3, 10), "beta": Fraction(1, 5)},
                          f"process parameters in {receipt}")
        data[receipt] = report
    process = data["results/independent_counterexample_process_audit.json"]
    for key, value in (("minimum_eigenvalue", "1/40"), ("witness_expectation", "-1/10"),
                       ("noise_CS_threshold_visibility", "10/11")):
        require_equal(Fraction(process[key]), Fraction(value), key)
    third = data["results/literature_scope/independent_two_certificate_ec_audit.json"]
    for name, certificate in (("canonical", "results/hall_sos_real_face_certificate.json"),
                              ("small_operator", "results/exact_double_flag_hall_sos.json")):
        require_equal(third[name]["accepted"], True, f"third audit {name}")
        check_hash(root, certificate, third[name]["certificate_sha256"])
    check_hash(root, "results/hpc_hall_margin_minimax_20261009/exact_canonical_layout.npz",
               third["canonical"]["layout_sha256"])
    require_equal(third["canonical"]["failed_exact_equations"], 0, "canonical stationarity")
    require_equal(Fraction(third["canonical"]["exact_dual_upper_bound"]), 0, "canonical dual objective")
    require_equal(third["small_operator"]["coefficient_residual_nonzero_words"], 0, "third SOS identity")
    theorem = data["results/input_bias_counterexample_theorem_audit.json"]
    dependencies = theorem["exact_certificate_dependencies"]
    check_hash(root, "results/exact_double_flag_hall_sos.json", dependencies["direct_certificate_sha256"])
    check_hash(root, "results/hall_sos_real_face_certificate.json", dependencies["canonical_certificate_sha256"])
    external = dependencies["external_purification_theorem"]
    check_hash(root, external["source"], external["source_sha256"])
    manual = theorem["manual_logical_audit"]
    require_equal(manual["obligation_count"], len(manual["obligations"]), "manual scope accounting")
    for key in ("arbitrary_finite_settings_and_results", "multi_Kraus_and_mixed_event_outputs",
                "general_Bob_effects_and_Naimark", "all_Hall_subsets", "cross_setting_and_order_weights",
                "mixed_shared_auxiliary_states", "unbounded_finite_auxiliary_dimensions",
                "general_nonpurifiability", "behavioral_causality_only"):
        if not manual["obligations"].get(key):
            raise ValueError(f"missing written logical obligation {key}")
    hall = theorem["symbolic_checks"]["Hall_transportation"]
    for key, value in (("q", "49/80"), ("theta", "3/25"), ("w_plus", "293/800"),
                       ("w_minus", "197/800"), ("c", "31/160")):
        require_equal(Fraction(hall[key]), Fraction(value), f"exact transport {key}")
    physical = data["results/literature_scope/independent_ec_physical_coverage_audit.json"]
    if not physical["cases"]:
        raise ValueError("missing finite Born convention controls")
    for case in physical["cases"]:
        require_equal(case["accepted"], True, "physical control case")
        if any(value > 1e-9 for key, value in case.items() if key.startswith("maximum_")):
            raise ValueError("finite Born convention control residual exceeds tolerance")
        if case["wrong_extra_transpose_negative_control_residual"] <= 1e-4:
            raise ValueError("Born transpose negative control is uninformative")
    lift = read_json(root, "results/independent_alice_event_lift_audit.json")
    check_hash(root, "repro/independent_alice_event_lift_audit.py", lift["source_sha256"])
    if lift["maximum_residual"] > 1e-9:
        raise ValueError("Alice event lift control failed")
    return {"independent_receipts_accepted_and_fresh": True,
            "receipt_paths": [relative for relative, _ in receipts],
            "canonical_stationarity_equation_count": third["canonical"]["equation_count"],
            "manual_logical_obligation_count": manual["obligation_count"],
            "finite_physical_control_cut_count": sum(case["subset_cuts_checked"] for case in physical["cases"]),
            "scope": "Source freshness and stored exact audits; written logical review is separately labelled."}


def compact_tex(text):
    text = re.sub(r"(?<!\\)%[^\n]*", "", text)
    text = re.sub(r"\\(?:notag|quad|qquad)\b|\\[ ,;:!]", "", text)
    text = re.sub(r"\s+|&|\\\\", "", text)
    return re.sub(r"\\frac(?:\{(-?\d+)\}|(-?\d))(?:\{(\d+)\}|(\d))",
                  lambda found: r"\frac{" + (found[1] or found[2]) + "}{" + (found[3] or found[4]) + "}", text)


def equation(text, label):
    blocks = re.findall(r"\\begin\{(equation\*?|align\*?)\}(.*?)\\end\{\1\}", text, re.S)
    found = [body for _, body in blocks if r"\label{" + label + "}" in body]
    require_equal(len(found), 1, f"unique equation {label}")
    return compact_tex(found[0])


def require_tokens(text, tokens, label):
    text = compact_tex(text)
    for token in tokens:
        token = compact_tex(token)
        if token not in text:
            raise ValueError(f"current EC manuscript missing {label}: {token}")


def check_manuscript_texts(main, supplement):
    for name, text, label in (("main", main, "eq:process"),
                              ("supplement", supplement, "eq:Wstar")):
        require_tokens(equation(text, label),
                       (r"\frac14", r"\frac15\sigma_z^{A_I}",
                        r"\frac3{10}\sigma_z^{A_I}\sigma_x^{B_I}\sigma_z^{B_O}",
                        r"\frac15\sum_{k=x,y,z}\sigma_k^{A_O}\sigma_k^{B_I}"),
                       f"{name} process")
        require_tokens(text, (r"\sigma_x,\sigma_y,\sigma_z",), f"{name} Pauli definition")
    require_tokens(equation(supplement, "eq:margins"),
                   (r"q=\frac{49}{80}", r"\theta=\frac3{25}", r"w_+=\frac{293}{800}",
                    r"w_-=\frac{197}{800}"), "fixed margins")
    require_tokens(main, (r"q=\frac{49}{80}",), "main-text order weight")
    for name, text in (("main", main), ("supplement", supplement)):
        require_tokens(equation(text, "eq:directional-response"),
                       (r"\langle b\rangle_x=\frac{x}{5}",
                        r"\langle s\rangle_y=\frac15+\frac{3y}{10}"), f"{name} operational responses")
    require_tokens(equation(main, "eq:causal-model"),
                   (r"p=q\,p^{A\prec B}+(1-q)\,p^{B\prec A}",), "causal probability model")
    require_tokens(equation(main, "eq:witness-value"), (r"\Tr(SW_*)=-\frac1{10}",), "strict witness")
    require_tokens(equation(main, "eq:noise-classification"),
                   (r"\mathrm{Purifiable}", r"\CS", r"\EC", r"\eta\ge\frac1{11}"), "exact noise threshold")
    require_tokens(equation(main, "eq:hall"),
                   (r"p(\mathcalA,\mathcalB|x,y)", r"\le", r"g(\mathcalA|x)+h(\mathcalB|y)"), "all-cut Hall criterion")
    for name, text in (("main", main), ("supplement", supplement)):
        require_tokens(equation(text, "eq:robust-center"),
                       (r"\overline W=\frac{21}{22}W_*+\frac1{22}W_0",), f"{name} neighborhood center")
        require_tokens(equation(text, "eq:robust-ball"),
                       (r"\|\Delta\|_\infty\le\frac1{264}", r"\overline W+\Delta", r"\EC\setminus\CS"),
                       f"{name} exact neighborhood radius")
    require_tokens(equation(supplement, "eq:CS-ball"), (r"\frac1{12}",), "CS inner radius")
    require_tokens(equation(supplement, "eq:robust-witness"), (r"-\frac{13}{660}",), "strict neighborhood margin")
    require_tokens(supplement, (r"\label{tab:proof-map}",), "physical proof map")
    quantifiers = {
        "arbitrary shared ancillary dimension": r"(?:every|all|arbitrary).{0,100}finite.{0,60}(?:shared|ancillar)",
        "arbitrary instruments": r"(?:every|all|arbitrary).{0,90}(?:local instruments|CP instruments)",
        "settings and outcomes": r"(?:settings.{0,60}(?:outcomes|results)|(?:outcomes|results).{0,60}settings)",
        "mixed states": r"mixed.{0,35}states?.{0,350}(?:convex|combine|purif)",
        "multiple Kraus events": r"(?:multiple|multi.?|all).{0,30}Kraus|Kraus.{0,30}environment",
        "Naimark coverage": r"Naimark",
        "all purification dimensions": r"(?:arbitrary|every|all).{0,35}finite.{0,50}past.{0,40}future|past.{0,40}future.{0,70}(?:arbitrary|every|all).{0,35}finite|finite.{0,60}past.{0,50}future.{0,200}unrestricted",
    }
    for name, text in (("main", main), ("supplement", supplement)):
        plain = re.sub(r"\s+", " ", text)
        for label, pattern in quantifiers.items():
            if re.search(pattern, plain, re.I) is None:
                raise ValueError(f"{name} lacks the declared quantifier coverage: {label}")
        normalized = compact_tex(text)
        constants = (r"\frac{49}{80}", r"\frac1{11}")
        if name == "supplement":
            constants += (r"\frac3{25}",)
        require_tokens(normalized, constants,
                       f"{name} exact constants")
    for relative in FRAGMENTS:
        filename = Path(relative).stem
        if re.search(r"\\input\{" + re.escape(filename) + r"(?:\.tex)?\}", supplement) is None:
            raise ValueError(f"supplement does not include {filename}")
    bibliographies = []
    for text in (main, supplement):
        found = re.findall(r"\\begin\{thebibliography\}.*?\\end\{thebibliography\}", text, re.S)
        require_equal(len(found), 1, "unique bibliography")
        bibliographies.append(found[0])
    require_equal(bibliographies[0], bibliographies[1], "main/supplement bibliography agreement")
    defined = re.findall(r"\\bibitem\{([^}]+)\}", bibliographies[0])
    require_equal(len(defined), len(set(defined)), "unique bibliography keys")
    cited = {key.strip() for group in re.findall(r"\\cite\{([^}]+)\}", main + supplement)
             for key in group.split(",")}
    undefined = cited - set(defined)
    if undefined:
        raise ValueError(f"undefined current citation keys: {sorted(undefined)}")
    require_tokens(bibliographies[0], (r"\bibitem{purif}", r"\bibitem{Yokojima}", r"\bibitem{EC}",
                                      r"\bibitem{Wei}", r"\bibitem{Wechs2026}", r"\bibitem{QCQC}",
                                      r"\bibitem{Trusted}", r"\bibitem{Network}"), "scope references")
    return {"counterexample_and_margin_transcriptions_checked": True,
            "exact_noise_threshold_checked": True, "declared_full_quantifier_coverage_checked": True,
            "full_direction_robustness_transcriptions_checked": True,
            "operational_responses_transcriptions_checked": True,
            "all_compact_fragments_included": True, "bibliographies_equal": True,
            "citation_keys_defined": True, "bibliography_entry_count": len(defined),
            "scope": "Textual coverage is checked here; logical coverage follows from the accompanying proof."}


def audit_current(root=ROOT):
    root = Path(root)
    evidence = check_evidence(root)
    compact = check_compact_artifacts(root)
    manuscripts = check_manuscript_texts((root / "main/main.tex").read_text(encoding="utf-8-sig"),
                                         (root / "supplement/supplement.tex").read_text(encoding="utf-8-sig"))
    inputs = ["main/main.tex", "supplement/supplement.tex", *FRAGMENTS,
              "results/exact_double_flag_hall_sos.json", "results/exact_double_flag_hall_sos_symmetric.json",
              "results/hall_sos_real_face_certificate.json", "results/hall_sos_compact_tables.json",
              "results/hall_sos_compact_ideal_multipliers.json", "results/hall_sos_positive_principal_minors.json",
              "repro/manuscript_science_audit.py", "repro/current_ec_manuscript_checks.py",
              *evidence["receipt_paths"]]
    return {"schema": "current-exact-EC-manuscript-science-audit-v1", "accepted": True,
            "utc": datetime.now(timezone.utc).isoformat(), "checks": {**manuscripts, **compact},
            "independent_evidence": evidence,
            "source_sha256": {relative: digest(root / relative) for relative in inputs},
            "proved_theorem": "The explicit full-rank bipartite qubit process is EC and nonpurifiable for all finite ancillary/purification dimensions; its white-noise CS threshold is 1/11, and W_(1/22) has a certified all-valid-direction operator-norm neighborhood of radius 1/264.",
            "not_proved_by_this_audit": ["The externally cited pure-process classification theorem",
                                         "A classification of every EC process", "Novelty completeness",
                                         "External expert approval or journal acceptance"]}
