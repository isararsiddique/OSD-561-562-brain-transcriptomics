"""Ensure the portable preparation reproduces the supplied audited mapping."""
import os
from pathlib import Path
import pandas as pd
from pipeline.paths import stage_root

expected = pd.read_csv(Path(os.environ["OSDR_INPUT_ROOT"]) / "reference/verified_panel_mapping.csv").fillna("")
actual = pd.read_csv(stage_root("core") / "tables/verified_panel_mapping.csv").fillna("")
pd.testing.assert_frame_equal(expected, actual)
print("Reproduced all 151 audited target mapping records.")
