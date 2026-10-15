# Confounding/comparator sensitivity extension

Run from the package root using its locked Python environment:

```
python extensions/confounding_sensitivity/scripts/analyse.py
python extensions/confounding_sensitivity/scripts/verify_independently.py
python extensions/confounding_sensitivity/scripts/build_report_figure.py
```

No parent inputs or existing analysis outputs are modified. The dated analysis lock and the separately dated methodological amendment must remain together. Read SCIENTIFIC_REPORT.md for interpretation; results/ contains complete bootstrap, covariate, moment, design and sensitivity tables. figures/ contains the native 183-mm vector PDF/SVG and RGB 600-dpi PNG plus caption and font metadata. Bias benchmarks are scenarios, not date-effect estimates; no causal rescue is claimed.
