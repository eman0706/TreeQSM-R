# ============================================================
# Adds the GPL notice (docs/INTESTAZIONE_GPL.txt) at the top of every R/*.R file that does not already carry a
# "This file is part of ... TREEQSM" line in its first 8 lines.
# Only comment lines are added before the first line of the file: no function is touched, so the md5 of the function
# bodies in docs/manifesto_validato.csv does not change. Each modified file is re-parsed; if it no longer parses,
# the original bytes are restored.
# Usage (any R session; run from the repository folder, or pass it):
#   source("add_gpl_headers.R")
#   add_gpl_headers("F:/EB/tartarino/TreeQSM-R")                   # dry run: only lists what would change
#   add_gpl_headers("F:/EB/tartarino/TreeQSM-R", dry_run = FALSE)  # writes
# ============================================================
add_gpl_headers <- function(repo = ".", dry_run = TRUE) {
  hfile <- file.path(repo, "docs", "INTESTAZIONE_GPL.txt")
  stopifnot(file.exists(hfile), dir.exists(file.path(repo, "R")))
  hdr <- readLines(hfile, warn = FALSE)
  files <- list.files(file.path(repo, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
  todo <- character(0); done <- character(0); failed <- character(0)
  for (f in files) {
    raw <- readBin(f, "raw", file.size(f))
    first <- readLines(f, n = 8, warn = FALSE)
    if (any(grepl("part of (the R port of )?TREEQSM", first, ignore.case = TRUE))) next
    todo <- c(todo, basename(f))
    if (dry_run) next
    crlf <- length(raw) > 1 && any(raw[-length(raw)] == as.raw(13) & raw[-1] == as.raw(10))
    eol <- if (crlf) "\r\n" else "\n"
    bom <- length(raw) >= 3 && all(raw[1:3] == as.raw(c(0xEF, 0xBB, 0xBF)))
    body <- if (bom) raw[-(1:3)] else raw
    head_bytes <- charToRaw(paste0(paste(hdr, collapse = eol), eol, eol))
    new <- c(if (bom) raw[1:3], head_bytes, body)
    writeBin(new, f)
    ok <- tryCatch({ parse(f, keep.source = FALSE); TRUE }, error = function(e) FALSE)
    if (ok) done <- c(done, basename(f)) else { writeBin(raw, f); failed <- c(failed, basename(f)) }
  }
  cat(sprintf("%d files lack the notice: %s\n", length(todo), paste(todo, collapse = ", ")))
  if (dry_run) cat("Dry run: nothing was written. Rerun with dry_run = FALSE.\n")
  else {
    cat(sprintf("Notice added to %d files.\n", length(done)))
    if (length(failed)) cat("RESTORED (did not parse after the change):", paste(failed, collapse = ", "), "\n")
  }
  invisible(list(todo = todo, done = done, failed = failed))
}
