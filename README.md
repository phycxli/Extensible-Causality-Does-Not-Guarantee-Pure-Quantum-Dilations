# Extensible Causality Does Not Guarantee Pure Quantum Dilations

Exact rational certificates and independent verification code accompanying the
manuscript by **Chengxi Li, Hao Zhu, Kai Liu, and Wanzi Sun**.
Correspondence: **Wanzi Sun**, wzsun@jnxy.edu.cn.

The work constructs a full-rank process on four qubit ports whose probability
records remain causal under arbitrary finite shared-input extensions, although
the process has no finite pure-process dilation. The original global past and
future are trivial. The complete argument, including the coverage of general
local instruments and arbitrary finite ancillary dimensions, is in the
Supplemental Material.

## Quick Start

Requirements: **Python 3.10 or later** and **SymPy 1.14.0**. No optimizer, GPU,
HPC account, or LaTeX installation is needed to verify the saved certificates.

```bash
git clone https://github.com/phycxli/Extensible-Causality-Does-Not-Guarantee-Pure-Quantum-Dilations.git
cd Extensible-Causality-Does-Not-Guarantee-Pure-Quantum-Dilations
python -m pip install -r environment/verification_requirements.txt
python -B verify.py --output-dir verification_results
```

Use the same Python executable for dependency installation and verification.
The commands also work in Windows PowerShell. A successful run exits with code
zero and ends with:

```text
ALL_EXACT_CHECKS_PASSED 25 frozen inputs
```

The verifier checks the input manifest both before and after execution. Results
are written to `verification_results/`, including `verification_receipt.json`.
An unsuccessful check raises an error or returns a nonzero exit status.
The repository attributes preserve the exact input bytes across operating
systems, including Windows checkouts with automatic line-ending conversion.

## What Is Verified

- The original and symmetry-averaged operator sum-of-squares identities,
  coefficient by coefficient, with exact rational arithmetic.
- Positive semidefiniteness of all rational certificate weights.
- Process normalization, linear constraints, positive spectrum, causal witness,
  and the exact white-noise causal-separability threshold.
- Two local laboratory response protocols and their complex-transpose control.
- The conservative operator-norm neighborhood of radius `1/264` at white-noise
  fraction `1/22`, using the complete valid Pauli basis.
- The printed compact matrices, positive principal minors, ideal multipliers,
  manuscript constants, and matching bibliographies.
- Seven boundary and rejection tests for the neighborhood claim.

The proof remains necessary for universal physical coverage, transportation
sufficiency, convexity, and the externally cited pure-process structure theorem.
The code checks the exact algebra and constants; the universal statements
follow from the accompanying proofs. The radius is a sufficient absolute operator-norm bound;
it is not a percentage of experimental error or an optimal radius.

## Repository Contents

```text
verify.py                      supported verification entry point
environment/                   pinned verification dependency
repro/                         independent exact verifiers
results/                       saved rational certificates and compact tables
tests/                         boundary and rejection controls
main/                          English manuscript, source, and figure PDFs
supplement/                    English proof and exact matrix TeX fragments
SHA256SUMS.txt                  checksum manifest for 25 frozen input files
CITATION.cff                   repository citation metadata
LICENSE                        existing MIT license supplied by the repository owner
```

This repository provides the exact verification distribution, with `verify.py`
as its supported entry point. Optional certificate searches, conic optimizers,
and figure-regeneration workflows are not required to reproduce these checks.

## Citation And License

Use GitHub's **Cite this repository** metadata or `CITATION.cff`. No DOI or formal
journal publication is claimed. When citing a particular code state, also record
the Git commit returned by `git rev-parse HEAD`.

The repository owner's existing MIT license is retained without modification;
see `LICENSE`. `README.txt` is the unchanged README supplied with the journal
code attachment and does not itself assign an additional license.
