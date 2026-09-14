# Source the pipeline's R/ files once for the whole suite.
#
# metadid-sims is not an installed package (see DESCRIPTION), so there is no
# namespace to load -- the helper sources R/ directly. Paths are resolved
# relative to the repo root so the suite runs from either the root or tests/.
suppressMessages({
  library(tibble); library(dplyr); library(purrr)
})

.sims_root <- local({
  d <- normalizePath(".", mustWork = FALSE)
  while (!file.exists(file.path(d, "DESCRIPTION")) && dirname(d) != d) d <- dirname(d)
  d
})

for (.f in list.files(file.path(.sims_root, "R"), pattern = "[.]R$", full.names = TRUE)) {
  # plots.R pulls in ggplot2 and friends but nothing here tests plotting;
  # skip it so the suite stays fast and dependency-light. R/archive.R is
  # deliberately separate from it for exactly that reason -- its functions
  # guard the pipeline against a malformed archive and must be tested.
  if (basename(.f) == "plots.R") next
  suppressMessages(source(.f))
}

# A minimal dgp list built from default_dgp, for tests that need one.
mk_dgp <- function(...) modifyList(default_dgp, list(...))
