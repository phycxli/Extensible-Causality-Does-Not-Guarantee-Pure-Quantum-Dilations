EXACT VERIFICATION: EXTENSIBLE CAUSALITY DOES NOT GUARANTEE PURE QUANTUM DILATIONS

Authors: Chengxi Li, Hao Zhu, and Wanzi Sun.
Corresponding author: Wanzi Sun, wzsun@jnxy.edu.cn.

Requirements: Python 3.10 or later and SymPy 1.14.0.
From the extracted distribution root:

  python -m pip install -r environment/verification_requirements.txt
  python -B verify.py --output-dir verification_results

Success: ALL_EXACT_CHECKS_PASSED and exit code zero.
The input SHA256 manifest is checked before and after verification.
The original and symmetric SOS identities, rational positive weights,
process constraints, spectrum, causal witness, noise threshold, laboratory
responses, and conservative all-direction neighborhood are checked.
No optimizer or HPC account is needed. The proof in the Supplemental
Material is needed for the universal physical coverage and the externally
cited pure-process theorem. Machine verification is not peer review.

This is the minimal exact verification distribution. It does not include
the optional conic-solver, search, or table-regeneration workflows listed
in the full research project. No reuse license is assigned by this file.
