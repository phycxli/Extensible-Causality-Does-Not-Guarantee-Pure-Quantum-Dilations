"""Single entry point for the portable exact verification distribution."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import unittest


ROOT = Path(__file__).resolve().parent
if not (ROOT / "main/main.tex").is_file():
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked_manifest(root):
    entries = {}
    for line in (root / "SHA256SUMS.txt").read_text(encoding="ascii").splitlines():
        expected, relative = line.split("  ", 1)
        path = root / relative
        if path.resolve().is_relative_to(root.resolve()) is False:
            raise ValueError("manifest path escapes the verification distribution")
        if digest(path) != expected:
            raise ValueError(f"input checksum mismatch: {relative}")
        entries[relative] = expected
    return entries


def verify(output):
    before = checked_manifest(ROOT)
    output.mkdir(parents=True, exist_ok=True)
    record = {
        "schema": "paper8-portable-exact-verification-v1",
        "utc": datetime.now(timezone.utc).isoformat(), "success": False,
        "job_id": os.environ.get("SLURM_JOB_ID", "local"), "hostname": platform.node(),
        "source_root": str(ROOT), "output_root": str(output),
        "manifest_sha256": digest(ROOT / "SHA256SUMS.txt"),
        "input_sha256": before, "audits": {},
    }
    os.chdir(ROOT)
    try:
        from repro.independent_general_ec_certificate_audit import audit as sos_audit
        from repro.independent_counterexample_process_audit import audit as process_audit
        from repro.independent_ec_robustness_audit import audit as robustness_audit
        from repro.independent_operational_response_audit import audit as response_audit
        from repro.independent_printed_word_identity import audit as printed_audit
        from repro.current_ec_manuscript_checks import check_compact_artifacts, check_manuscript_texts

        for name, source in (
            ("original_SOS", "results/exact_double_flag_hall_sos.json"),
            ("symmetric_SOS", "results/exact_double_flag_hall_sos_symmetric.json"),
        ):
            target = output / (name + ".json")
            sos_audit(Path(source), target)
            record["audits"][name] = json.loads(target.read_bytes())["accepted"]
        process_audit(output / "process_witness_noise.json")
        robustness_audit(output / "robustness.json")
        response_audit(output / "operational_responses.json")
        for name in ("process_witness_noise", "robustness", "operational_responses"):
            record["audits"][name] = json.loads((output / (name + ".json")).read_bytes())["accepted"]
        record["transcriptions"] = check_manuscript_texts(
            (ROOT / "main/main.tex").read_text(encoding="utf-8-sig"),
            (ROOT / "supplement/supplement.tex").read_text(encoding="utf-8-sig"),
        )
        record["compact_tables"] = check_compact_artifacts(ROOT)
        record["printed_word_identity"] = printed_audit(root=ROOT, output=output / "printed_word_identity.json")
        suite = unittest.defaultTestLoader.loadTestsFromName("tests.test_ec_robustness")
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        record["rejection_tests"] = {"run": result.testsRun, "passed": result.wasSuccessful()}
        record["input_sha256_after"] = checked_manifest(ROOT)
        record["inputs_unchanged"] = record["input_sha256_after"] == before
        record["success"] = (
            all(record["audits"].values()) and record["printed_word_identity"]["accepted"] and result.wasSuccessful()
            and result.testsRun >= 7 and record["inputs_unchanged"]
        )
    except Exception as error:
        record["error"] = str(error)
        raise
    finally:
        (output / "verification_receipt.json").write_text(
            json.dumps(record, indent=2) + "\n", encoding="ascii", newline="\n"
        )
    if not record["success"]:
        raise SystemExit(1)
    print("ALL_EXACT_CHECKS_PASSED", len(before), "frozen inputs")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "verification_results")
    args = parser.parse_args()
    verify(args.output_dir.resolve())
