#!/usr/bin/env Rscript
# Get arguments
args <- commandArgs(trailingOnly = TRUE)

# Help document
print_help <- function() {
  cat("Description:\n")
  cat(" Conversion tip-label of newick tree using look-up table.\n\n")
  cat("Usage:\n")
  cat(" convert_tip_lab.r -i <in_tree> -l <in_lookup> [-o <out_tree>] [-v]\n\n")
  cat("Options:\n")
  cat("  -i, --in_tree <path>   Input tree (required) \n")
  cat("  -l, --in_lookup <path> ID-lookup table (required) \n")  
  cat("  -o, --out_tree <path>  Output tree（default: 'renamed_tree.nwk') \n")
  cat("  -v, --verbose          Enable detailed output(default:FALSE) \n")
  cat("  -h, --help             Print this documents\n\n")
  cat("Example: \n")
  cat(" convert_tip_lab.r -i ./input_tree.nwk -o renamed_tree.nwk -v \n\n")
}

# Default values
in_tree <- NULL
in_lookup <- NULL
out_tree <- "renamed_tree.nwk"
verbose <- FALSE

# Analysis optional arguments
i <- 1
while (i <= length(args)) {
  if (args[i] %in% c("-h", "--help")) {
    print_help()
    quit(status = 0)
  } else if (args[i] %in% c("-i", "--in_tree")) {
    in_tree <- args[i + 1]
    i <- i + 1
  } else if (args[i] %in% c("-t", "--in_lookup")) {
    in_lookup <- args[i + 1]
    i <- i + 1
  } else if (args[i] %in% c("-o", "--out_tree")) {
    out_tree <- args[i + 1]
    i <- i + 1
  } else if (args[i] %in% c("-v", "--verbose")) {
    verbose <- TRUE
  } else {
    cat("[ERROR] Unexpected argument:", args[i], "\n")
    print_help()
    quit(status = 1)
  }
  i <- i + 1
}

# Validate input
if (is.null(in_tree) || is.null(in_lookup)) {
  print_help()
  quit(status = 1)
}

# Validate output tree file
if (file.exists(out_tree)) {
  cat("[ERROR] The output tree already exists.\n")
  quit(status = 1)
}

# Verbose
if (verbose) {
  cat("[INFO] Input tree file :", in_tree, "\n")
  cat("[INFO] Input id-lookup table :", in_lookup, "\n")
  cat("[INFO] Output tree file :", out_tree, "\n")
}

library(ape)
rename_tree <- function(infile, lookup) {
  # Read newick tree
  tree <- ape::read.tree(infile)
  # Index of tip-labels
  idx <- match(tree$tip.label, lookup$old)
  # Validate missing tips
  if (anyNA(idx)) {
    missing <- tree$tip.label[is.na(idx)]
    stop(
      "Missing lookup entries in ", infile, ": ",
      paste(missing, collapse = ", ")
    )
  }
  tree$tip.label <- lookup$new[idx]
  tree

}

# Read lookup table
lookup <- read.delim(
  in_lookup,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
colnames(lookup) <- c("old", "new")

if (anyDuplicated(lookup$old) || anyDuplicated(lookup$new)) {
  stop("Duplicate IDs or replacement labels in lookup table.")
}

renamed_tree <- rename_tree(infile = in_tree, lookup = lookup)

# Write newick tree file
ape::write.tree(renamed_tree, file = out_tree)