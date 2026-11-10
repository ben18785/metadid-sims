# Build the cross-category targets and render the validation report, from a
# store already merged by scripts/merge-targets-store.sh.
#
# Shared by sims.yaml (at the end of a full run) and render-report.yaml (on a
# previous run's uploaded slices) so the target list cannot drift between them.
library(targets)

# Anything outdated beyond this list means the store does not match the code or
# the installed metadid, and tar_make is about to re-run scenarios rather than
# just aggregate. Printed rather than enforced: in sims.yaml a rebuild is
# legitimate recovery for a slice that failed to upload, while in
# render-report.yaml CmdStan is absent, so an attempted refit fails loudly.
outdated <- tar_outdated()
cat("Outdated targets (", length(outdated), "):\n", sep = "")
print(outdated)

tar_make(names = c(
  "scenario_lookup_tbl",
  "diag_plot_A", "diag_plot_B", "diag_plot_C", "diag_plot_D",
  "diag_plot_E", "diag_plot_F", "diag_plot_G", "diag_plot_H",
  "diag_plot_I", "diag_plot_J", "diag_plot_K", "diag_plot_L",
  "diag_plot_M", "diag_plot_N",
  "all_agg", "all_rep",
  "export_agg_csv", "export_rep_rds", "archive_results",
  "scenario_settings_csv", "figure_data_csv",
  "report"
))
