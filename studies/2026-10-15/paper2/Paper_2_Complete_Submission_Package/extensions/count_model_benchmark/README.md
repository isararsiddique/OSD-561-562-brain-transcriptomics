# Count-model benchmark

This extension preserves a frozen retrospective plan, all923 clusters/all2769 tests, full-depth DM and multinomial-score exact inference, all 10 original A-half full-matrix sensitivities, independent checks, and an editable four-panel figure. It does not change causal or independent-validation conclusions.

Run with the package Python environment, from any directory:

```sh
python /path/to/package/extensions/count_model_benchmark/scripts/analyse_benchmark.py
python /path/to/package/extensions/count_model_benchmark/scripts/verify_independently.py
python /path/to/package/extensions/count_model_benchmark/scripts/build_report_figure.py
```

See SCIENTIFIC_REPORT.md for complete methods/results and every limitation. No external downloads are needed to reproduce numerical analyses; source metadata/documentation are bundled.

Offline regeneration can change gzip header timestamps in derived inputs. The verifier retains the original compressed SHA256 hashes and accepts a regenerated gzip only when its decompressed bytes exactly match the frozen payload SHA256 in inputs/frozen_payload_sha256.json. Raw inputs and other file types still require their original byte hash.
