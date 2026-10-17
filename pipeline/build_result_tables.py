"""Table-derived consensus and summaries used by figures and the report."""
import pandas as pd
from pipeline.paths import stage_root

B, R = stage_root("deep"), stage_root("extended")
old = pd.read_csv(R / "tables/sorted_cell_panel_complete.csv")
new = pd.read_csv(B / "tables/article_5pct_panel_complete.csv")
select = lambda data: data[data.tissue.eq("Myeloid_Forebrain") & data.contrast.eq("flight_balanced_age")]
data = pd.concat([select(old), select(new)], ignore_index=True)
consensus = data.groupby("analyte").agg(n=("q_full_3624_family", "size"),
    max_q=("q_full_3624_family", "max"), min_FC=("log2FC", "min"), max_FC=("log2FC", "max"))
consensus["all_eight_corrected"] = consensus.n.eq(8) & consensus.max_q.lt(.05)
mapping = pd.read_csv(stage_root("core") / "tables/verified_panel_mapping.csv")
consensus = consensus.reset_index().merge(mapping[["analyte", "mouse_symbol_primary"]], on="analyte", validate="one_to_one")
assert len(consensus) == 151 and consensus.n.eq(8).all()
consensus.to_csv(B / "tables/forebrain_all_eight_QC_method_consensus.csv", index=False)
bulk = pd.read_csv(B / "tables/additional_panel_all_regions.csv")
bulk[bulk.analyte.isin(["AIM2", "PYCARD", "AIF1", "SLC1A3", "S100B"])].to_csv(B / "tables/reported_additional_focus_results.csv", index=False)
deletion = pd.read_csv(B / "tables/additional_mouse_deletion_expanded_focus.csv")
summary = deletion[deletion.OSD.eq(deletion.deleted_OSD)].groupby(["OSD", "analyte", "method", "contrast"]).agg(
    n_deletions=("deleted_mouse", "size"), corrected_deletions=("q_expanded_1812_family", lambda q: int(q.lt(.05).sum())),
    maximum_q=("q_expanded_1812_family", "max"), minimum_FC=("log2FC", "min"), maximum_FC=("log2FC", "max"))
summary.to_csv(B / "tables/reported_expanded_deletion_summary.csv")
print("Created complete 151-target consensus and deletion summary tables.")
