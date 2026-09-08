# scripts/gene_modules.R
# Deterministic assignment of gene symbols to functional modules, used to group
# (and order) the genes on the top-25 DEG dot plot in 03_dea_and_plotting.R.
#
# The table below is the single source of truth. Its ORDER matters twice:
#   1. it breaks ties during assignment (first matching module wins), and
#   2. it fixes the left-to-right order of the facet panels on the plot.
#
# Assignment is a two-pass, first-match-wins rule (see assign_gene_modules):
#   pass 1 - literal symbol match against `exact`
#   pass 2 - regex match against `regex`, for whatever pass 1 left over
# Exact-before-regex is what keeps the table maintainable: an exception such as
# S100B (glial, not an S100A-family calcium binder) is fixed by listing the
# symbol explicitly in the right module, without reshuffling the panel order.
#
# Pure string matching - no annotation database, no network, no run-to-run drift.
#
# The module names ARE the facet strip labels. Keep them SHORT: a module that
# ends up holding a single gene gets a one-gene-wide panel, and a long label
# will overhang it (strip.clip = "off" in 03 keeps it readable rather than
# truncated, but it can still collide with its neighbours). Use "\n" to wrap.

GENE_MODULES <- list(
  "Antigen presentation" = list(
    exact = c("CD74", "B2M", "CIITA"),
    regex = c("^HLA-", "^CD1[ABCDE]$")
  ),
  "Glial / neural" = list(
    exact = c("GFAP", "AQP4", "PTPRZ1", "S100B", "MBP", "PLP1", "SOX2", "SOX9",
              "ALDOC", "CLU", "ATP1A2", "SPARCL1", "SLC1A2", "SLC1A3", "FABP7",
              "GJA1", "NTRK2", "MT3", "BCAN", "TTYH1", "PDGFRA"),
    regex = c("^OLIG", "^NEF")
  ),
  "Immune receptors" = list(
    exact = c("HCST", "TYROBP", "TREM1", "TREM2", "MRC1", "CD163", "MSR1",
              "MARCO", "FPR1", "FPR2", "C5AR1", "CD14", "ITGAM", "CSF1R"),
    regex = c("^FC[EGA]R", "^CLEC", "^TLR[0-9]", "^SIGLEC", "^SIRP")
  ),
  "Complement" = list(
    # NB: C1orf162 is an uncharacterised myeloid ORF, NOT a complement gene -
    # do not re-add it here on the strength of the "C1" prefix.
    # ^SERPIN is deliberately absent: it also catches coagulation serpins.
    exact = c("C1QA", "C1QB", "C1QC", "C3", "C3AR1", "C5AR2", "CFD", "FCN1",
              "SERPING1"),
    regex = c("^C1Q")
  ),
  "Coagulation" = list(
    exact = c("F13A1", "PLAUR", "PLAU", "TFPI", "FGL2", "THBD", "SERPINE1"),
    regex = character(0)
  ),
  "Cytokines\n& chemokines" = list(
    exact = c("TNF", "OSM", "SPP1", "IL1B", "IL1RN"),
    regex = c("^CCL[0-9]", "^CXCL[0-9]", "^IL[0-9]+$", "^TGFB[0-9]", "^CSF[0-9]")
  ),
  "Cytokine\nsignalling" = list(
    # PKIB (cAMP/PKA inhibitor) and PSTPIP2 (F-BAR adaptor) were removed:
    # neither is cytokine signalling, both now fall through to "Other".
    exact = c("CD52"),
    regex = c("^IL[0-9]+R", "^JAK[0-9]", "^STAT[0-9]", "^SOCS[0-9]", "^NFKB")
  ),
  "Growth factors" = list(
    exact = c("AREG", "EREG", "HBEGF", "TGFA", "VEGFA"),
    regex = c("^FGF[0-9]", "^PDGF")
  ),
  "Interferon" = list(
    exact = c("IRF1", "IRF7", "IRF8", "GBP1", "GBP2"),
    regex = c("^IFI", "^ISG[0-9]", "^MX[12]$", "^OAS[123]")
  ),
  "S100 / Ca2+\n& annexins" = list(
    exact = c("CAPG", "CALM1", "CALM2"),
    regex = c("^S100", "^ANXA")
  ),
  "Lysosomal" = list(
    # ^LGALS removed: galectins are secreted immunomodulatory lectins,
    # neither lysosomal nor antimicrobial.
    exact = c("LYZ", "CD68", "GPNMB", "APOE", "APOC1", "LIPA", "ACP5", "NPC2"),
    regex = c("^CTS[A-Z]")
  ),
  "Cytoskeleton\n& adhesion" = list(
    exact = c("VIM", "JAML", "CD44", "ACTB", "TMSB4X", "TMSB10", "CORO1A"),
    regex = c("^ACT[BG]", "^TUB[AB]", "^MYO[0-9]", "^ITG", "^SEL[EPL]$")
  ),
  "Ribosomal\n& translation" = list(
    exact = character(0),
    regex = c("^RP[LS]", "^EEF[12]", "^EIF[0-9]")
  ),
  "Metabolism" = list(
    exact = c("GAPDH", "LDHA", "SLC2A1", "PKM", "ENO1"),
    regex = c("^MT-", "^NDUF", "^COX[0-9]", "^ATP5", "^HK[123]$")
  )
)

# Catch-all bucket. Always rendered as the last panel.
MODULE_OTHER <- "Other"

# Module names become facet panel identities, so they must be unique: a duplicate
# would misassign the second entry's genes to "Other" and emit the panel twice,
# which fails inside DotPlot with "factor level [n] is duplicated".
if (anyDuplicated(names(GENE_MODULES)) > 0) {
  stop("Duplicated module name(s) in GENE_MODULES: ",
       paste(unique(names(GENE_MODULES)[duplicated(names(GENE_MODULES))]), collapse = ", "))
}
if (MODULE_OTHER %in% names(GENE_MODULES)) {
  stop("GENE_MODULES must not define a module named '", MODULE_OTHER,
       "' - it is the reserved catch-all bucket.")
}


#' Assign gene symbols to functional modules.
#'
#' Two passes, first match wins in each: literal symbols first, then regex.
#' Anything unmatched falls into MODULE_OTHER.
#'
#' @param genes character vector of HGNC symbols
#' @param modules ordered named list, see GENE_MODULES
#' @return character vector, same length and order as `genes`
assign_gene_modules <- function(genes, modules = GENE_MODULES) {
  assigned <- rep(NA_character_, length(genes))

  # Pass 1: literal symbol match.
  for (mod in names(modules)) {
    exact <- modules[[mod]]$exact
    if (length(exact) == 0) next
    hit <- is.na(assigned) & genes %in% exact
    assigned[hit] <- mod
  }

  # Pass 2: regex match on whatever is left.
  for (mod in names(modules)) {
    patterns <- modules[[mod]]$regex
    if (length(patterns) == 0) next
    for (pat in patterns) {
      hit <- is.na(assigned) & grepl(pat, genes)
      assigned[hit] <- mod
    }
  }

  assigned[is.na(assigned)] <- MODULE_OTHER
  assigned
}


#' Build the named list that Seurat::DotPlot turns into facet panels.
#'
#' Seurat 5 facets on `features` when it is a named list: panel order is the
#' list order, within-panel order is the vector order. It also does
#' `names(feature.groups) <- features`, so a gene must not appear twice -
#' guaranteed here because assign_gene_modules() gives each gene one module.
#'
#' @param genes character vector, already in the desired within-panel order
#' @param modules_assigned module per gene, as returned by assign_gene_modules()
#' @param modules ordered named list, used to order the panels
#' @return named list of character vectors; empty modules dropped, "Other" last
build_module_feature_list <- function(genes, modules_assigned,
                                      modules = GENE_MODULES) {
  stopifnot(length(genes) == length(modules_assigned))
  keep <- !duplicated(genes)
  genes <- genes[keep]
  modules_assigned <- modules_assigned[keep]

  panel_order <- c(names(modules), MODULE_OTHER)
  present <- panel_order[panel_order %in% modules_assigned]

  feature_list <- lapply(present, function(mod) genes[modules_assigned == mod])
  names(feature_list) <- present
  feature_list
}
