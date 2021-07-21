#! /usr/env/Rscript
# This is a script to run the R component of TrackSig pipeline, it improves 
# calling the functions and specifying the inputs over the initial setup
# AUTHOR (of the pipeline and original scripts): Yulia Rubanova
# Author (of the modifications): Kasia Kedzierska

########################################
# OPTIONS AND TESTING                  # 
########################################
# Load necessary packages for script execution
suppressMessages(library("optparse"))

option_list = list(
  make_option(c("-i", "--input_dir"), 
              type = "character", 
              default = NULL, 
              help = paste("Path to the input directory."), 
              metavar = "character"),
  make_option(c("-a", "--annotation_dir"), 
              type = "character", 
              default = NULL, 
              help = paste("Path with signature, active_signature and",
                           "trinucleotide files"), 
              metavar = "character"),
  make_option(c("-c", "--cellularity"), 
              type = "character", 
              default = NULL, 
              help = "Path to cellularity file", 
              metavar = "character"),
  make_option(c("-s", "--subclone"), 
              type = "character", 
              default = NULL, 
              help = "Path to subclone file", 
              metavar = "character"),
  make_option(c("-n", "--sample_name"), 
              type = "character", 
              default = NA, 
              help = paste("Sample name, if not specified, assumes", 
                           "sample-name_rest_of_the_input_snvs_name.vcf.gz",
                           "and exctracts sample-name only."), 
              metavar = "character"),
  make_option(c("-o", "--output_dir"), 
              type = "character", 
              default = ".", 
              help = paste("Output directory in which subdirectory", 
                           "for the sample will be creater. Default: ."), 
              metavar = "character"),
  make_option(c("-v", "--verbose"), 
              action = "store_true", 
              default = FALSE, 
              help = "Print more information")
) 

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


if (opt$verbose) {
  library("logger")
  log_threshold(TRACE)
} else {
  suppressMessages(library("logger"))
  log_threshold(WARN)
}

log_layout(layout_glue_colors)

`%not in%` <- Negate(`%in%`)

# First check if the essential packages are installed
necessary_packages <- c("reshape2", "ggplot2", "NMF")

for (pkg in necessary_packages) {
  if (pkg %not in% rownames(installed.packages())) {
    msg <- paste0(pkg, " not installed! Checked following .libPaths(): ",
                  paste(.libPaths(), collapse = ", "), ".")
    log_fatal(msg)
    stop(msg)
  }
  if (verbose) {
    library(pkg)
    log_debug("Package {pkg} loaded")
  } else {
    suppressMessages(pkg)
  }
}

# Assign variables based on input and test if they are ok

if (opt$simulate_data) {
  simulated_data <- TRUE
  dir_counts <- file.path(opt$input_dir, "counts")
  mutation_order <- NULL
  BOOTSTRAP_COUNTS <- NULL
  purity_file <- NULL
} else {
  dir_counts <- file.path(opt$input_dir, "counts")
  mutation_order <- file.path(opt$input_dir, "mut_order")
  BOOTSTRAP_COUNTS <- file.path(opt$input_dir, "bootstrap")
  purity_file <- file.path(opt$input_dir, "purity.tsv")
}

# if the signatures are specified per cancer type or per sample
cancer_type_signatures = TRUE
# if signatures trajectories need to be computed on bootstrapped signatures as well
# bootstrapping provides the uncertainty estimations on the trajectories
# warning: by default, mutations are bootstrapped 30 times and the script will run 30 time longer
compute_bootstrap = FALSE

sliding_window = FALSE
noise_sig = NULL

# specifies the changepoint detection algorithm.
changepoint_method = "PELT"

# file with cancer types of each sample
tumortype_file <- "data/tumortypes.txt"



# folder to write results to
DIR_RESULTS = opt$output_dir #"results_signature_trajectories/"


# file with signatures definitions
signature_file <- file.path(opt$annotation_dir, "signatures.txt")

# file with trinucleotide context
trinucleotide_file <- file.path(opt$annotation_dir, "trinucleotide.txt")

# specifies active signatures in TCGA cancer types
active_signatures_file <- file.path(opt$annotation_dir,
                                    "active_signatures_transposed.txt")


SAVED_SAMPLES_DIR = "saved_data/"

PLOT_FULL_NAME <- TRUE
mutation_assignments = ""

if (!file.exists(DIR_RESULTS)) {
  dir.create(DIR_RESULTS, recursive = TRUE)
}

src_files <- 
  setdiff(grep(".*R$", list.files(paste0( "src"), full.names = T), value = T), 
                     c(paste0( "src/compute_mutational_signatures.R"),
                       paste0( "src/header.R")))

for (file in src_files) {
  source(file)
}



result_list[alex, 
            tumortypes, 
            active_signatures, 
            active_signatures.our_samples] <- 
  load_annotation(tumortype_file, signature_file, active_signatures_file)

save_data_for_samples()
suppressMessages(compute_signatures_for_all_examples())
if (compute_bootstrap) {
  compute_errorbars_for_all_examples()
}
