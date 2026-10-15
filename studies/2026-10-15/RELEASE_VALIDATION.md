# Validation performed for this GitHub release

The following checks were run during release preparation on 15 October 2026.

| Check | Result |
|---|---|
| SHA-256 and byte count for every distributed archive member | 1,214 files matched their source archives |
| Python source syntax, without executing module-level analyses | 174 archived Python files parsed |
| Dataset exclusion | No input/data/cache/source/marker directories or compressed/binary dataset formats in the distributed study trees |
| Local links in the new study guides | 15 links resolved |
| Complete Paper 1 archive restoration | SHA-256 verified; 1,520 files restored to the local reproduction workspace |
| Complete Paper 2 archive identity | SHA-256 verified using the restoration helper's check-only mode |
| Existing Paper 1 statistical regression suite on the restored archive | All 11 tests passed |

The statistical suite checks full-design/FWL equivalence, cluster sandwich calculations, whole-animal allocation and region linkage, finite permutation resolution, score/operator identities, spectral gain constraints, invalid/unidentifiable inputs, and frozen-input integrity with matching STAR/RSEM technical-library collapse.

That suite ran in the existing analysis environment: Python 3.12.14, NumPy 2.3.5, pandas 2.2.3 and SciPy 1.16.3. This is the environment used for this release check; it is not a claim that every archived dependency specification or every full analysis stage was rerun in that environment. Original dependency/runtime records remain unchanged.

Both full analysis archives passed ZIP integrity checks before upload preparation. Paper 1's analysis-only Zenodo derivative retained 1,515 existing files byte-for-byte, omitted seven superseded publication documents, retained the original scientific verification report and regenerated its integrity manifests. Paper 2's upload ZIP is unchanged from its existing data/code/results archive.

Full Paper 2 statistical reconstruction, all strengthened Paper 1 model fits, remote data downloads, online literature/enrichment services and publication-document rendering were not rerun for this GitHub release. Their original code, complete archive outputs and dated independent-check records remain available. The new checks verify distribution integrity and the stated statistical primitives; they do not replace the complete study-specific validation histories.

Archived CSV/source files retain their original line endings and whitespace to preserve hashes. `.gitattributes` protects those copies from Git text normalization; the newly written guides are checked as ordinary text.
