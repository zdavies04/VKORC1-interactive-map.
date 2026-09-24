# =============================================================================
# VKORC1 Anticoagulant Rodenticide Resistance — Interactive Global Map
# =============================================================================
# Author:  [Zak Davies]
# Data:   Published meta-analysis
# Run:     shiny::runApp("vkorc1_mutation_map.R")
#          OR open in RStudio and click "Run App"
# Deploy:  rsconnect::deployApp("vkorc1_mutation_map.R")
#
# IMPORTANT: Keep Vkorc1_clean.csv in the SAME folder as this script.
#
# ── 0. INSTALL / LOAD PACKAGES ───────────────────────────────────────────────

required_packages <- c(
  "shiny", "leaflet", "dplyr", "readr", "RColorBrewer",
  "htmltools", "shinyWidgets", "DT", "leaflet.extras", "stringr"
)
installed  <- rownames(installed.packages())
to_install <- required_packages[!required_packages %in% installed]
if (length(to_install) > 0) install.packages(to_install)

library(shiny);        library(leaflet);      library(dplyr)
library(readr);        library(RColorBrewer); library(htmltools)
library(shinyWidgets); library(DT);           library(leaflet.extras)
library(stringr)


# ── 1. LOAD DATA ─────────────────────────────────────────────────────────────

df <- read_csv("Vkorc1_clean.csv", show_col_types = FALSE)

# Stable per-site ID (survives filtering) used to link table rows -> map popups
df$Row_ID <- seq_len(nrow(df))

# Standardise spelling and abbreviation variants without merging genuinely
# distinct taxa. Keep the reported wording for audit and export.
normalise_species <- function(x) {
  original <- trimws(as.character(x))
  key <- tolower(gsub("\\s+", " ", original))
  key <- gsub("\\.", "", key)

  case_when(
    grepl("macedonicus", key) ~ "Mus macedonicus",
    grepl("domesticus", key) ~ "Mus musculus domesticus",
    grepl("castaneus", key) ~ "Mus musculus castaneus",
    key == "m m musculus" ~ "Mus musculus musculus",
    key %in% c("m musculus", "mus musculus") ~ "Mus musculus",
    TRUE ~ original
  )
}

df <- df %>%
  mutate(
    Species_Reported = Species,
    Species = normalise_species(Species)
  )

# Derive Resistance_Level directly from the mutations/evidence encoded at a site.
# Site-level evidence is aggregated from the mutation-level classifications below.
classify_resistance <- function(mutations, rt = NULL) {
  x <- as.character(mutations)
  if (is.na(x) || trimws(x) == "" || grepl("^no mutations$", trimws(x), ignore.case = TRUE)) {
    return("None")
  }

  d <- classify_mutations_detailed(x)
  if (nrow(d) == 0) return("None")

  levels <- unique(na.omit(as.character(d$level)))
  if (length(levels) == 0) return("Unknown")
  if (identical(levels, "Confirmed")) return("Confirmed")
  if (identical(levels, "Suspected")) return("Suspected")
  if (identical(levels, "Unknown")) return("Unknown")
  "Mixed evidence"
}
set.seed(42)   # reproducible jitter

mutation_aliases <- c(
  "Arg12Trp"  = "R12W",
  "Ala48Thr"  = "A48T",
  "Gly2Ala"   = "G2A",
  "Ala14Thr"  = "A14T",
  "Leu20Ile"  = "L20I",
  "Ala21Thr"  = "A21T",
  "Ala26Pro"  = "A26P",
  "Ala26Ser"  = "A26S",
  "Ala26Thr"  = "A26T",
  "Ala32Val"  = "A32V",
  "Asp36His"  = "D36H",
  "Glu37Gly"  = "E37G",
  "Ala48Thr"  = "A48T",
  "Ser52Tyr"  = "S52Y",
  "Phe55Val"  = "F55V",
  "Arg58Gly"  = "R58G",
  "Arg58Trp"  = "R58W",
  "Trp59Gly"  = "W59G",
  "Trp59Ser"  = "W59S",
  "Trp59Arg"  = "W59R",
  "Trp59Cys"  = "W59C",
  "Trp59Leu"  = "W59L",
  "Arg61Leu"  = "R61L",
  "Arg61Gln"  = "R61Q",
  "Arg61Trp"  = "R61W",
  "Ala72Val"  = "A72V",
  "Gln78His"  = "Q78H",
  "Cys85Arg"  = "C85R",
  "Phe87Leu"  = "F87L",
  "Ile104Val" = "I104V",
  "Val114Phe" = "V114F",
  "Ala115Thr" = "A115T",
  "Val118Leu" = "V118L",
  "Leu124Met" = "L124M",
  "Leu124Gln" = "L124Q",
  "Leu128Ser" = "L128S",
  "Ser149Ile" = "S149I",
  "Ser149Asn" = "S149N",
  "Gln151His" = "Q151H",
  "Glu155Lys" = "E155K",
  "Lys157Asn" = "K157N",
  "Tyr139Cys" = "Y139C",
  "Tyr139Phe" = "Y139F",
  "A26Ser"    = "A26S",
  "A26Thr"    = "A26T"
)

# Human-readable names used in the expanded table dropdown/details.
mutation_display_names <- c(
  "G2A"  = "Gly2Ala",
  "R12W" = "Arg12Trp",
  "A14T" = "Ala14Thr",
  "L20I" = "Leu20Ile",
  "A21T" = "Ala21Thr",
  "A26P" = "Ala26Pro",
  "A26S" = "Ala26Ser",
  "A26T" = "Ala26Thr",
  "A32V" = "Ala32Val",
  "D36H" = "Asp36His",
  "E37G" = "Glu37Gly",
  "A48T" = "Ala48Thr",
  "S52Y" = "Ser52Tyr",
  "F55V" = "Phe55Val",
  "R58G" = "Arg58Gly",
  "R58W" = "Arg58Trp",
  "W59G" = "Trp59Gly",
  "W59S" = "Trp59Ser",
  "W59R" = "Trp59Arg",
  "W59C" = "Trp59Cys",
  "W59L" = "Trp59Leu",
  "R61L" = "Arg61Leu",
  "R61Q" = "Arg61Gln",
  "R61W" = "Arg61Trp",
  "A72V" = "Ala72Val",
  "Q78H" = "Gln78His",
  "C85R" = "Cys85Arg",
  "F87L" = "Phe87Leu",
  "I104V" = "Ile104Val",
  "V114F" = "Val114Phe",
  "A115T" = "Ala115Thr",
  "V118L" = "Val118Leu",
  "L124M" = "Leu124Met",
  "L124Q" = "Leu124Gln",
  "L128S" = "Leu128Ser",
  "S149I" = "Ser149Ile",
  "S149N" = "Ser149Asn",
  "Q151H" = "Gln151His",
  "E155K" = "Glu155Lys",
  "K157N" = "Lys157Asn",
  "Y139C" = "Tyr139Cys",
  "Y139F" = "Tyr139Phe"
)

normalise_mutation_names <- function(x) {
  out <- as.character(x)
  for (nm in names(mutation_aliases)) {
    out <- gsub(nm, mutation_aliases[[nm]], out, ignore.case = TRUE, perl = TRUE)
  }
  out
}


# Create ONE dropdown value for the complete mutation combination reported in
# the CSV.  A multi-mutant is any mutation followed by one or more mutations
# using either "+" or "/" (for example A26S+L128S or A26S/L128S).
# The dropdown/table display is normalised to "+" so the same combination is
# not split into separate single-mutation entries.
make_primary_mutation <- function(x) {
  x <- as.character(x)
  if (is.na(x) || trimws(x) == "") return("Other/Complex")
  if (grepl("^no mutations$", trimws(x), ignore.case = TRUE)) return("No Mutations")

  # VKORC1^spr is a named genotype containing these four mutations.
  # Expand it BEFORE tokenisation so it behaves exactly like the explicit
  # combination R12W+A26S+A48T+R61L.
  x <- gsub(
    "VKORC1\\^?spr",
    "R12W+A26S+A48T+R61L",
    x,
    ignore.case = TRUE,
    perl = TRUE
  )

  # Normalise full amino-acid names to canonical one-letter notation first.
  x <- normalise_mutation_names(x)

  # Remove genotype/sample-count annotations before extracting mutations.
  x <- gsub("\\b[0-9]+\\s*(?:HOM|HET)\\b", "", x,
            ignore.case = TRUE, perl = TRUE)
  x <- gsub("\\b(?:HOM|HET)\\b", "", x,
            ignore.case = TRUE, perl = TRUE)

  tokens <- toupper(str_extract_all(x, "\\b[A-Z][0-9]+(?:[A-Z])?\\b")[[1]])
  tokens <- tokens[tokens %in% names(mutation_display_names)]
  tokens <- unique(tokens)

  if (length(tokens) == 0) return("Other/Complex")

  # A multi-mutant requires an actual separator between mutation tokens.
  has_mutation_separator <- grepl(
    "(?:[A-Z][0-9]+[A-Z])\\s*(?:\\+|/|,|;|\\band\\b)\\s*(?:[A-Z][0-9]+[A-Z])",
    x, ignore.case = TRUE, perl = TRUE
  )

  if (length(tokens) >= 2 && has_mutation_separator) {
    # Sort by biological residue number, then retain deterministic ordering.
    nums <- as.integer(str_extract(tokens, "\\d+"))
    tokens <- tokens[order(nums, toupper(tokens), method = "radix")]
    return(paste(tokens, collapse = "+"))
  }

  tokens[[1]]
}

df <- df %>%
  mutate(
    # Primary_Mutation is the COMPLETE mutation combination used by the
    # dropdown/filter.  Multi-mutants therefore appear once, e.g.
    # A26S+L128S+Y139C, rather than as repeated single mutations.
    Primary_Mutation = vapply(Mutations, make_primary_mutation, character(1)),
    # The visible table intentionally says "Multiple" for a multi-mutant.
    # The complete combination remains available in Primary_Mutation and
    # in the expanded mutation details.
    Mutation_Display = if_else(
      str_detect(Primary_Mutation, fixed("+")),
      "Multiple",
      Primary_Mutation
    ),
    DOI_URL  = if_else(str_starts(DOI, "http"), DOI, NA_character_),
    Lat_j    = Latitude  + runif(n(), -0.06, 0.06),
    Lon_j    = Longitude + runif(n(), -0.06, 0.06)
  )

# Extract year from Date_Collected
extract_year <- function(d) {
  m <- str_extract(d, "\\d{4}")
  if (!is.na(m)) as.integer(m) else NA_integer_
}
df$Year <- sapply(df$Date_Collected, extract_year)
df$Year[is.na(df$Year)] <- 2020L


# ── 1b. PER-MUTATION CLASSIFICATION (Confirmed / Suspected / Unknown) ───────
# These lookup tables are the source of truth for mutation-level evidence.
# Site-level Resistance_Level is then derived from the complete set of
# mutation-level classifications at that site.

mutation_level_colours <- c(
  "Confirmed" = "#e74c3c",
  "Suspected" = "#e67e22",
  "Unknown"   = "#95a5a6"
)

# Single mutations with a known classification (uppercase, no spaces/punct.)
single_mutation_levels <- c(
  # Mutations listed in the source data but with no resistance evidence
  # classification encoded in this app are explicitly retained as Unknown.
  "G2A"  = "Unknown",
  "A14T" = "Unknown",
  "L20I" = "Unknown",
  "A21T" = "Unknown",
  "A26P" = "Unknown",
  "A32V" = "Unknown",
  "D36H" = "Unknown",
  "S52Y" = "Unknown",
  "F55V" = "Unknown",
  "R58W" = "Unknown",
  "W59R" = "Unknown",
  "W59C" = "Unknown",
  "W59L" = "Unknown",
  "R61Q" = "Unknown",
  "R61W" = "Unknown",
  "A72V" = "Unknown",
  "Q78H" = "Unknown",
  "C85R" = "Unknown",
  "F87L" = "Unknown",
  "I104V" = "Unknown",
  "V114F" = "Unknown",
  "A115T" = "Unknown",
  "V118L" = "Unknown",
  "S149I" = "Unknown",
  "S149N" = "Unknown",
  "Q151H" = "Unknown",
  "E155K" = "Unknown",
  "K157N" = "Unknown",
  "L128S"  = "Confirmed",
  "Y139C"  = "Confirmed",
  "Y139F"  = "Confirmed",
  "A26S"   = "Suspected",
  "A26T"   = "Suspected",
  "E37G"   = "Suspected",
  "A48T"   = "Suspected",
  "R58G"   = "Suspected",
  "W59G"   = "Suspected",
  "W59S"   = "Suspected",
  "R61L"   = "Suspected",
  "L124M"  = "Suspected",
  "L124Q"  = "Suspected"
)

# Full-name aliases found in the CSV are normalised to standard one-letter
# amino-acid notation before combo matching/tokenisation.
# Multi-mutation combinations that must be recognised as a single evidence unit.
# They are checked longest/most-specific first, before generic token extraction.
combo_patterns <- list(
  list(
    label = "Arg12Trp+Ala26Ser+Ala48Thr+Arg61Leu (VKORC1^spr)",
    level = "Confirmed",
    regex = "R12W\\s*(?:\\+|/)\\s*A26S\\s*(?:\\+|/)\\s*A48T\\s*(?:\\+|/)\\s*R61L"
  ),
  list(
    label = "A26T+L128S",
    level = "Suspected",
    regex = "A26T\\s*(?:\\+|/)\\s*L128S"
  ),
  list(
    label = "A26S+L128S",
    level = "Suspected",
    regex = "A26S\\s*(?:\\+|/)\\s*L128S"
  ),
  list(
    label = "W59G+L124M",
    level = "Suspected",
    regex = "W59G\\s*(?:\\+|/)\\s*L124M"
  ),
  list(
    label = "W59G+L128S",
    level = "Suspected",
    regex = "W59G\\s*(?:\\+|/)\\s*L128S"
  )
)

# Splits a site's raw Mutations string into evidence units with a level.
# Known combinations above stay together; other mutation tokens are classified
# individually, so e.g. L128S+Y139C contributes two confirmed mutations.
classify_mutations_detailed <- function(mutations_string) {
  x <- as.character(mutations_string)
  if (is.na(x) || trimws(x) == "" || grepl("^no mutations$", trimws(x), ignore.case = TRUE)) {
    return(data.frame(mutation = character(0), level = character(0), stringsAsFactors = FALSE))
  }

  remaining <- normalise_mutation_names(x)
  found <- list()

  # 1) pull out known multi-mutation combos so they are not split apart
  for (pat in combo_patterns) {
    mm <- regmatches(remaining, gregexpr(pat$regex, remaining, ignore.case = TRUE, perl = TRUE))[[1]]
    if (length(mm) > 0 && !(length(mm) == 1 && identical(mm[[1]], "-1"))) {
      for (i in seq_along(mm)) {
        if (identical(mm[[i]], "")) next
        found[[length(found) + 1]] <- data.frame(
          mutation = pat$label,
          level = pat$level,
          stringsAsFactors = FALSE
        )
      }
      remaining <- gsub(pat$regex, " , ", remaining, ignore.case = TRUE, perl = TRUE)
    }
  }

  # 2) extract individual mutation tokens from anything left. This intentionally
  # ignores counts/HET/HOM text and handles comma, slash, plus and mixed formats.
  tokens <- str_extract_all(remaining, "[A-Z][0-9]+[A-Z]")[[1]]
  tokens <- toupper(tokens)

  if (length(tokens) > 0) {
    for (p in tokens) {
      lvl <- unname(single_mutation_levels[p])
      if (is.na(lvl) || length(lvl) == 0) lvl <- "Unknown"
      found[[length(found) + 1]] <- data.frame(
        mutation = p,
        level = lvl,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(found) == 0) {
    return(data.frame(mutation = "No classified mutations", level = "Unknown", stringsAsFactors = FALSE))
  }
  do.call(rbind, found)
}

# Build the mutation entries shown in the table dropdown.
# IMPORTANT: commas/semicolons/"and" separate reported entries, but + and /
# stay inside the same entry so multi-mutants are displayed as combinations.
extract_mutation_entries <- function(mutations_string) {
  x <- as.character(mutations_string)
  if (is.na(x) || trimws(x) == "" || grepl("^no mutations$", trimws(x), ignore.case = TRUE)) {
    return(character(0))
  }

  # Normalise mutation aliases first.
  x <- normalise_mutation_names(x)

  # IMPORTANT: remove genotype/sample-count annotations before looking for
  # mutations. The source data can contain entries such as "3 HOM", "28 HOM",
  # "4 HET", etc. These are not mutations and must never appear in the
  # mutation dropdown.
  x <- gsub("\\b[0-9]+\\s*(?:HOM|HET)\\b", "", x, ignore.case = TRUE, perl = TRUE)
  x <- gsub("\\b(?:HOM|HET)\\b", "", x, ignore.case = TRUE, perl = TRUE)

  parts <- str_split(x, "[,;]|\\s+and\\s+", simplify = FALSE)[[1]]
  parts <- trimws(parts)
  parts <- parts[parts != ""]

  vapply(parts, function(part) {
    # Expand VKORC1^spr into its four constituent mutations so the
    # displayed entry is the full biological combination.
    part <- gsub(
      "VKORC1\\^?spr",
      "R12W+A26S+A48T+R61L",
      part,
      ignore.case = TRUE,
      perl = TRUE
    )

    # Extract only valid mutation tokens. These are restricted to the
    # standard one-letter amino-acid notation used by single_mutation_levels
    # and mutation_display_names. This prevents counts, HOM/HET annotations,
    # and other text from becoming dropdown entries.
    tokens <- toupper(str_extract_all(part, "\\b[A-Z][0-9]+(?:[A-Z])?\\b")[[1]])
    tokens <- tokens[tokens %in% names(mutation_display_names)]
    tokens <- unique(tokens)

    # If the part contains no recognised mutation, discard it rather than
    # displaying arbitrary source text as a mutation.
    if (length(tokens) == 0) return(NA_character_)

    token_names <- vapply(tokens, function(tok) {
      unname(mutation_display_names[tok])
    }, character(1))

    if (length(token_names) == 1) token_names[[1]] else paste(token_names, collapse = "+")
  }, character(1)) |>
    stats::na.omit() |>
    as.character() |>
    unique()
}

# Same as extract_mutation_entries() above, but returns each entry as its
# canonical one-letter-code combination (e.g. "R12W+A26S+A48T+R61L") instead
# of the human-readable display name. This is the SOURCE OF TRUTH for the
# Mutations dropdown and for filtering: unlike Primary_Mutation (which can
# incorrectly merge multiple, separately reported entries at one site into a
# single fabricated combination), this respects the comma/semicolon/"and"
# boundaries between distinct reported entries, so a site listing
# "Trp59Leu+Tyr139Phe, Leu128Ser+Tyr139Cys" yields TWO separate entries
# rather than one fake four-mutation combo.
extract_mutation_entries_codes <- function(mutations_string) {
  x <- as.character(mutations_string)
  if (is.na(x) || trimws(x) == "" || grepl("^no mutations$", trimws(x), ignore.case = TRUE)) {
    return(character(0))
  }

  x <- normalise_mutation_names(x)
  x <- gsub("\\b[0-9]+\\s*(?:HOM|HET)\\b", "", x, ignore.case = TRUE, perl = TRUE)
  x <- gsub("\\b(?:HOM|HET)\\b", "", x, ignore.case = TRUE, perl = TRUE)

  parts <- str_split(x, "[,;]|\\s+and\\s+", simplify = FALSE)[[1]]
  parts <- trimws(parts)
  parts <- parts[parts != ""]

  vapply(parts, function(part) {
    part <- gsub(
      "VKORC1\\^?spr",
      "R12W+A26S+A48T+R61L",
      part,
      ignore.case = TRUE,
      perl = TRUE
    )

    tokens <- toupper(str_extract_all(part, "\\b[A-Z][0-9]+(?:[A-Z])?\\b")[[1]])
    tokens <- tokens[tokens %in% names(mutation_display_names)]
    tokens <- unique(tokens)

    if (length(tokens) == 0) return(NA_character_)

    # Order by residue number for a stable, canonical combo key.
    nums <- as.integer(str_extract(tokens, "\\d+"))
    tokens <- tokens[order(nums, toupper(tokens), method = "radix")]

    paste(tokens, collapse = "+")
  }, character(1)) |>
    stats::na.omit() |>
    as.character() |>
    unique()
}

# Does a site (identified by its list of reported entry-codes) match ANY of
# the (possibly several) selected dropdown values?
#  - no values selected               -> nothing matches
#  - a selected value is a full combo (contains "+") -> must match one of the
#    site's reported entries EXACTLY (so picking "R12W+A26S+A48T" never
#    matches a site whose only entry is "R12W+A26S+A48T+R61L")
#  - a selected value is a single mutation code -> matches if that code is
#    one of the components of ANY of the site's reported entries
mutation_entry_matches <- function(entry_codes, filter_mutation) {
  # No mutations selected -> nothing matches (consistent with the other
  # multi-select filters: deselecting everything shows an empty result).
  if (length(filter_mutation) == 0) return(FALSE)

  # A site with NO reported mutations only matches if the special "None"
  # option ("__NONE__") is among the selected values.
  if (length(entry_codes) == 0) return("__NONE__" %in% filter_mutation)

  # TRUE if the site matches ANY of the (possibly several) selected values,
  # so a multi-select acts as an OR across the chosen mutations/combos.
  any(vapply(filter_mutation, function(fm) {
    if (fm %in% entry_codes) return(TRUE)
    if (!grepl("+", fm, fixed = TRUE)) {
      components <- unlist(strsplit(entry_codes, "+", fixed = TRUE))
      return(fm %in% components)
    }
    FALSE
  }, logical(1)))
}

# Per-site list of correctly-separated reported mutation entries (one-letter
# codes), e.g. list("R12W+A26S+A48T", "R12W+A26S+A48T+R61L") for two sites
# that report those as distinct entries. This is what the Mutations dropdown
# and its filter are built from.
df$Entry_Codes <- lapply(df$Mutations, extract_mutation_entries_codes)


# Number of different reported mutation entries at a site. A multi-mutant is
# one displayed entry, e.g. Leu128Ser+Tyr139Cys counts as one entry.
summarise_levels <- function(mutations_string) {
  entries <- extract_mutation_entries(mutations_string)
  as.integer(length(unique(entries)))
}

# Assign the site-level category now that the mutation classifier is defined.
df$Resistance_Level <- vapply(
  df$Mutations,
  function(x) classify_resistance(mutations = x),
  character(1)
)

# ── 2. COLOUR PALETTE ────────────────────────────────────────────────────────
# None is reserved for sites explicitly reporting no mutations. Unknown means
# mutation(s) are present but none have a known resistance classification.
# Mixed evidence means the site contains more than one evidence level.
resistance_levels  <- c("None", "Suspected", "Confirmed", "Mixed evidence", "Unknown")
resistance_colours <- c(
  "None"           = "#27ae60",
  "Suspected"      = "#e67e22",
  "Confirmed"      = "#e74c3c",
  "Mixed evidence" = "#8e44ad",
  "Unknown"        = "#95a5a6"
)
get_colour <- function(level) unname(resistance_colours[level])

# ── 3. POPUP BUILDER ─────────────────────────────────────────────────────────

make_popup <- function(row) {
  doi_html <- if (!is.na(row$DOI_URL)) {
    sprintf('<a href="%s" target="_blank" style="color:#2980b9;">%s</a>',
            htmlEscape(row$DOI_URL), htmlEscape(row$DOI_URL))
  } else "Not available"

  sprintf(
    '<div style="font-family:Segoe UI,Arial,sans-serif;font-size:13px;
                 max-width:320px;line-height:1.55;">
       <b style="font-size:14px;color:#2c3e50;">%s</b><br/>
       <span style="color:#7f8c8d;font-size:11px;"><i>%s</i></span>
       <hr style="margin:5px 0;border-color:#ecf0f1;"/>
       <b>Country:</b> %s &nbsp;|&nbsp; <b>Location:</b> %s<br/>
       <b>Sample size:</b> %s &nbsp;|&nbsp; <b>Site:</b> %s<br/>
       <b>Date collected:</b> %s<br/>
       <b>Resistance level:</b>
         <span style="color:%s;font-weight:bold;">%s</span><br/>
       <b>Resistance types:</b> %s<br/>
       <b>All mutations:</b> <code style="font-size:11px;">%s</code><br/>
       <b>Mutation classification:</b> %s<br/>
       <b>%% with mutation:</b> %s<br/>
       <b>DOI:</b> %s
     </div>',
    htmlEscape(row$Mutations),
    htmlEscape(row$Species),
    htmlEscape(row$Country),
    htmlEscape(row$Location),
    htmlEscape(as.character(row$Sample_Size)),
    htmlEscape(as.character(row$Sample_Site)),
    htmlEscape(as.character(row$Date_Collected)),
    get_colour(row$Resistance_Level),
    htmlEscape(row$Resistance_Level),
    htmlEscape(as.character(row$Resistance_Type)),
    htmlEscape(row$Mutations),
    summarise_levels(row$Mutations),
    htmlEscape(as.character(row$Percent_With_Mutation)),
    doi_html
  )
}


# Order mutation dropdown by amino-acid residue number rather than alphabetically.
# Multi-mutants are ordered by the first residue number, then subsequent residues.
sort_mutations <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & x != ""]

  mutation_number_key <- function(m) {
    nums <- suppressWarnings(as.integer(str_extract_all(m, "\\d+")[[1]]))
    if (length(nums) == 0) return("999999")
    paste(sprintf("%06d", nums), collapse = ".")
  }

  keys <- vapply(x, mutation_number_key, character(1))
  x[order(keys, toupper(x), method = "radix")]
}

# Human-readable labels for the Mutations dropdown. The underlying
# value remains canonical (e.g. A26T or A26T+L128S), so filtering is exact.
primary_mutation_label <- function(x) {
  if (x %in% c("No Mutations", "Other/Complex")) return(x)
  tokens <- strsplit(x, "\\+", fixed = FALSE)[[1]]
  labels <- vapply(tokens, function(tok) {
    if (tok %in% names(mutation_display_names)) unname(mutation_display_names[tok])
    else tok
  }, character(1))
  paste(labels, collapse = "+")
}

# ── 4. UI ─────────────────────────────────────────────────────────────────────

ui <- fluidPage(

  tags$head(tags$style(HTML("
    body { font-family:'Segoe UI',Arial,sans-serif; background:#f0f2f5; }
    .title-bar {
      background:linear-gradient(135deg,#1a252f,#2c3e50);
      color:white; padding:14px 20px; margin-bottom:14px;
      border-radius:8px;
    }
    .title-bar h3 { margin:0; font-size:17px; }
    .title-bar p  { margin:2px 0 0; font-size:11px; opacity:.75; }
    .sidebar-panel {
      background:#fff; border-radius:8px;
      padding:15px; box-shadow:0 2px 8px rgba(0,0,0,.08);
    }
    .sidebar-panel h4 { color:#2c3e50; font-weight:700; margin:0 0 12px; font-size:14px; }
    .filter-label { font-weight:600; color:#555; font-size:12px; display:block; margin-bottom:3px; }
    .stats-bar {
      background:#2c3e50; color:white;
      padding:8px 14px; border-radius:6px;
      margin-bottom:10px; font-size:12px;
    }
    #map_wrap { border-radius:8px; overflow:hidden; box-shadow:0 2px 12px rgba(0,0,0,.12); }
    .section-title { color:#2c3e50; font-weight:700; font-size:14px; margin:14px 0 6px; }
    .export-btn { width:100%; margin-top:6px; font-size:12px; }
    hr.thin { margin:10px 0; border-color:#ecf0f1; }
    .loc-link {
      color:#2980b9; cursor:pointer; font-weight:600;
      text-decoration:underline dotted;
    }
    .loc-link:hover { color:#e74c3c; }
  "))),

  div(class = "title-bar",
    tags$h3("🧬  VKORC1 Anticoagulant Resistance — Global Mutation Map"),
    uiOutput("title_stats")
  ),

  sidebarLayout(

    sidebarPanel(width = 3,
      div(class = "sidebar-panel",
        h4("🔍 Filters"),

        tags$label("Country", class = "filter-label"),
        pickerInput("filter_country", label = NULL,
                    choices  = sort(unique(df$Country)),
                    selected = sort(unique(df$Country)),
                    multiple = TRUE,
                    options  = pickerOptions(
                      actionsBox = TRUE,
                      liveSearch = TRUE,
                      selectedTextFormat = "count > 3",
                      countSelectedText = "{0} countries selected",
                      noneSelectedText = "No countries selected"
                    )),

        tags$br(),
        tags$label("Mutations", class = "filter-label"),
        pickerInput("filter_mutation", label = NULL,
                    multiple = TRUE,
                    options  = pickerOptions(
                      actionsBox = TRUE,
                      liveSearch = TRUE,
                      selectedTextFormat = "count > 3",
                      countSelectedText = "{0} mutations selected",
                      noneSelectedText = "No mutations selected"
                    ),
                    choices = {
  # Build the dropdown from BOTH:
  #   1) every known single mutation in mutation_display_names, and
  #   2) every complete mutation combination ACTUALLY REPORTED as its own
  #      entry in the CSV (via Entry_Codes, which respects the
  #      comma/semicolon/"and" boundaries between distinct reported
  #      entries at a site instead of merging them together).
  #
  # This is deliberately NOT based on Primary_Mutation, which could
  # fabricate combinations that were never reported by merging separate
  # entries listed at the same site (e.g. a site reporting both
  # "Trp59Leu+Tyr139Phe" and "Leu128Ser+Tyr139Cys" would otherwise produce
  # a non-existent "Trp59Leu+Leu128Ser+Tyr139Cys+Tyr139Phe" entry).
  #
  # A single mutation is only offered as its own dropdown entry if it is
  # actually REPORTED STANDALONE somewhere in the data (i.e. one of its
  # Entry_Codes has no "+"). A mutation that only ever occurs inside a
  # multi-mutant combo (e.g. Ala32Val, which is only ever reported as part
  # of Ala32Val+Leu128Ser or Ala32Val+Tyr139Cys) is not offered as a
  # standalone option - it can still be selected via its combo entries.
  vals <- unique(unlist(df$Entry_Codes))
  vals <- vals[!is.na(vals) & vals != ""]
  vals <- sort_mutations(vals)

  # IMPORTANT: setNames(object, nm) makes `nm` the DISPLAYED label and
  # `object` the SUBMITTED value (input$filter_mutation). The displayed
  # label must be the human-readable three-letter name, while the
  # submitted value must stay as the one-letter code so it matches
  # Entry_Codes exactly in filtered_data() below.
  mutation_labels <- vapply(vals, primary_mutation_label, character(1))

  # Prepend a "None" option so sites that explicitly report no mutations
  # (empty Entry_Codes) can be included in the map/table.
  c(
    stats::setNames("__NONE__", "None (no mutations)"),
    stats::setNames(
      vals,             # submitted value  -> one-letter code(s), e.g. "R12W" or "R12W+A26S+A48T+R61L"
      mutation_labels   # displayed label  -> three-letter name(s), e.g. "Arg12Trp"
    )
  )
},
                    selected = {
                      vals <- unique(unlist(df$Entry_Codes))
                      vals <- vals[!is.na(vals) & vals != ""]
                      c("__NONE__", sort_mutations(vals))
                    }),

        tags$br(),
        tags$label("Resistance Level", class = "filter-label"),
        selectInput("filter_resistance", label = NULL,
                    choices = c("All resistance levels" = "__ALL__", resistance_levels),
                    selected = "__ALL__"),

        tags$br(),
        tags$label("Species", class = "filter-label"),
        selectInput("filter_species", label = NULL,
                    choices = c("All taxa" = "__ALL__", sort(unique(df$Species))),
                    selected = "__ALL__"),

        tags$br(),
        tags$label("Year Range", class = "filter-label"),
        sliderInput("filter_year", label = NULL,
                    min = min(df$Year, na.rm = TRUE),
                    max = max(df$Year, na.rm = TRUE),
                    value = range(df$Year, na.rm = TRUE),
                    step = 1, sep = "",
                    animate = animationOptions(interval = 900, loop = FALSE)),

        tags$hr(class = "thin"),
        materialSwitch("cluster_toggle", "Cluster nearby points",
                       value = TRUE, status = "primary"),
        tags$hr(class = "thin"),
        downloadButton("export_csv", "⬇ Export Filtered Data",
                       class = "export-btn btn-sm btn-outline-secondary")
      )
    ),

    mainPanel(width = 9,
      uiOutput("stats_bar"),
      div(id = "map_wrap", leafletOutput("map", height = "560px")),
      p(class = "section-title", "📋 Filtered Data Table   (click a location to view it on the map)"),
      DTOutput("data_table")
    )
  )
)


# ── 5. SERVER ─────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Filtered reactive dataset
  filtered_data <- reactive({
    req(input$filter_resistance, input$filter_species, input$filter_year)
    df %>% filter(
      Country %in% input$filter_country,
      vapply(Entry_Codes, mutation_entry_matches, logical(1),
             filter_mutation = input$filter_mutation),
      input$filter_resistance == "__ALL__" | Resistance_Level == input$filter_resistance,
      input$filter_species == "__ALL__" | Species == input$filter_species,
      Year             >= input$filter_year[1],
      Year             <= input$filter_year[2]
    )
  })

  # Live counter in the title bar — always reflects the actual CSV, not a
  # hardcoded number. Recomputed from df (the full, unfiltered dataset)
  # every time the app starts / the data changes.
  output$title_stats <- renderUI({
    tags$p(sprintf(
      "Interactive meta-analysis  ·   %s sample sites · %s countries · %s unique mutations · Years %d\u2013%d   ·   Click a location in the table to jump to it on the map",
      format(nrow(df), big.mark = ","),
      n_distinct(df$Country),
      n_distinct(unlist(df$Entry_Codes)),
      min(df$Year, na.rm = TRUE),
      max(df$Year, na.rm = TRUE)
    ))
  })

  # Stats bar
  output$stats_bar <- renderUI({
    fd <- filtered_data()
    if (nrow(fd) == 0) return(div(class = "stats-bar", "Showing 0 sample sites"))
    div(class = "stats-bar",
        sprintf("Showing %d sample sites  ·  %d countries  ·  %d unique mutations  ·  Years %d – %d",
                nrow(fd), n_distinct(fd$Country),
                n_distinct(unlist(fd$Entry_Codes)),
                min(fd$Year, na.rm = TRUE), max(fd$Year, na.rm = TRUE)))
  })

  # Base map
  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
      
      addTiles(
        urlTemplate = paste0(
          "https://{s}.basemaps.cartocdn.com/rastertiles/light_all/{z}/{x}/{y}.png?key=",
          "cb1_3tuy_1_ec5589503bed11492c7ed6e2"
        ),
        attribution = '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a>, &copy; <a href="https://carto.com/attributions">CARTO</a>',
        options = tileOptions(
          opacity = 0.92,
          maxZoom = 20
        )
      ) %>%
      
      setView(lng = 10, lat = 45, zoom = 3) %>%
      
      addLegend(
        position = "bottomright",
        colors   = unname(resistance_colours),
        labels   = names(resistance_colours),
        title    = "Resistance Level",
        opacity  = 0.9
      )
  })

  # Update markers
  observe({
    fd <- filtered_data()
    leafletProxy("map") %>% clearMarkerClusters() %>% clearMarkers()
    if (nrow(fd) == 0) return()

    popups <- sapply(seq_len(nrow(fd)), function(i) make_popup(fd[i, ]))

    if (isTRUE(input$cluster_toggle)) {
      leafletProxy("map") %>%
        addCircleMarkers(
          data = fd, lat = ~Lat_j, lng = ~Lon_j,
          layerId = ~as.character(Row_ID),
          radius = 7, color = "#2c3e50", weight = 1,
          fillColor = ~sapply(Resistance_Level, get_colour),
          fillOpacity = 0.85,
          popup = popups,
          label = ~paste0(Primary_Mutation, " · ", Country),
          clusterOptions = markerClusterOptions(
            maxClusterRadius    = 45,
            showCoverageOnHover = FALSE,
            # Fixed colour regardless of how many points are grouped, so a
            # cluster bubble is never mistaken for a resistance-level colour
            iconCreateFunction  = JS("
              function(cluster){
                var n=cluster.getChildCount();
                var bg='#34495e';
                return new L.DivIcon({
                  html:'<div style=\"background:'+bg+';color:#fff;border-radius:50%;'+
                       'width:36px;height:36px;display:flex;align-items:center;'+
                       'justify-content:center;font-weight:700;font-size:12px;'+
                       'border:2px solid #fff;box-shadow:0 1px 5px rgba(0,0,0,.3);\">'+n+'</div>',
                  className:'',iconSize:[36,36]});
              }
            ")
          )
        )
    } else {
      leafletProxy("map") %>%
        addCircleMarkers(
          data = fd, lat = ~Lat_j, lng = ~Lon_j,
          layerId = ~as.character(Row_ID),
          radius = 7, color = "#2c3e50", weight = 1,
          fillColor = ~sapply(Resistance_Level, get_colour),
          fillOpacity = 0.85,
          popup = popups,
          label = ~paste0(Primary_Mutation, " · ", Country)
        )
    }
  })

  # Clicking a location in the table pans/zooms to it and opens its popup
  observeEvent(input$table_location_click, {
    clicked_id <- suppressWarnings(as.numeric(input$table_location_click))
    req(!is.na(clicked_id))
    clicked <- df[df$Row_ID == clicked_id, ]
    if (nrow(clicked) == 0) return()
    clicked <- clicked[1, ]

    leafletProxy("map") %>%
      clearPopups() %>%
      setView(lng = clicked$Lon_j, lat = clicked$Lat_j, zoom = 8) %>%
      addPopups(lng = clicked$Lon_j, lat = clicked$Lat_j,
                popup = make_popup(clicked), layerId = "table_click_popup")
  })

  # Data table
  output$data_table <- renderDT({
    fd <- filtered_data() %>%
      mutate(DOI_Link = if_else(!is.na(DOI_URL),
               sprintf('<a href="%s" target="_blank">🔗</a>', DOI_URL), "—"))

    # Clickable location -> triggers table_location_click in JS below
    fd$Location_Link <- sprintf(
      '<span class="loc-link" data-rowid="%s" title="Click to view on map">%s</span>',
      fd$Row_ID, htmlEscape(fd$Location)
    )

    # Build the dropdown from the reported mutation entries, preserving
    # multi-mutants exactly as combinations. For example, a site containing
    # Leu128Ser, Tyr139Cys and Leu128Ser+Tyr139Cys shows all three entries.
    mut_detail_list <- lapply(fd$Mutations, function(x) {
      entries <- extract_mutation_entries(x)
      if (length(entries) == 0) {
        return(data.frame(mutation = character(0), level = character(0), stringsAsFactors = FALSE))
      }

      data.frame(
        mutation = unique(entries),
        level = vapply(unique(entries), function(entry) {
          # Classify every component of a multi-mutant for its colour.
          tokens <- toupper(str_extract_all(normalise_mutation_names(entry), "[A-Z][0-9]+[A-Z]")[[1]])
          lvls <- unname(single_mutation_levels[tokens])
          lvls <- lvls[!is.na(lvls)]
          if (length(lvls) == 0) "Unknown" else if (all(lvls == "Confirmed")) "Confirmed" else if (all(lvls == "Suspected")) "Suspected" else "Unknown"
        }, character(1)),
        stringsAsFactors = FALSE
      )
    })

    # Numeric value so DT sorts 1, 2, 3, 10... correctly.
    fd$Levels_Summary <- vapply(fd$Mutations, summarise_levels, integer(1))

    n_mut <- vapply(mut_detail_list, nrow, integer(1))
    fd$Detail <- ifelse(n_mut > 0,
                         '<span style="color:#2980b9;font-weight:bold;">&#9654;</span>', "")

    fd$Mutation_Details <- vapply(mut_detail_list, function(d) {
      if (nrow(d) == 0) return("")
      items <- sprintf(
        '<li style="margin-bottom:3px;">
           <span style="display:inline-block;width:10px;height:10px;border-radius:50%%;
                        background:%s;margin-right:6px;"></span>%s
           <span style="color:%s;font-size:10px;font-weight:600;"> (%s)</span>
         </li>',
        unname(mutation_level_colours[d$level]),
        htmlEscape(d$mutation),
        unname(mutation_level_colours[d$level]),
        d$level
      )
      paste0(
        '<b>Mutations at this site:</b>',
        '<ul style="margin:4px 0 0 18px;padding:0;list-style:none;">',
        paste0(items, collapse = ""),
        '</ul>'
      )
    }, character(1))

    fd <- fd %>%
      select(Detail, Country, Location_Link, Species, Mutation_Display, Resistance_Level,
             Levels_Summary, Sample_Size, Date_Collected, Percent_With_Mutation, DOI_Link,
             Mutation_Details)

    datatable(
      fd, escape = FALSE, rownames = FALSE,
      class   = "stripe hover compact",
      options = list(
        pageLength = 15, scrollX = TRUE,
        columnDefs = list(
          list(className = "dt-center", targets = "_all"),
          list(orderable = FALSE, className = "details-control",
               targets = 0, width = "24px"),
          list(visible = FALSE, targets = 11)   # hide the raw Mutation_Details column
        )
      ),
      colnames = c("", "Country","Location","Species","Mutation","Resistance",
                   "Number of different mutations","N","Date","% Mutant","DOI","Mutation_Details"),
      callback = JS("
        // Expand/collapse per-site mutation breakdown
        table.on('click', 'td.details-control', function() {
          var td  = $(this);
          var tr  = td.closest('tr');
          var row = table.row(tr);
          if (row.child.isShown()) {
            row.child.hide();
            td.html('&#9654;');
          } else {
            var details = row.data()[11];
            if (details && details !== '') {
              row.child(
                '<div style=\"padding:10px 18px;background:#f8f9fa;' +
                'border-left:3px solid #2980b9;\">' + details + '</div>'
              ).show();
              td.html('&#9660;');
            }
          }
        });

        // Clicking a location jumps to it on the map
        table.on('click', '.loc-link', function(e) {
          e.stopPropagation();
          var id = $(this).data('rowid');
          Shiny.setInputValue('table_location_click', id, {priority: 'event'});
        });
      ")
    ) %>%
      formatStyle("Resistance_Level",
                  backgroundColor = styleEqual(names(resistance_colours),
                                               unname(resistance_colours)),
                  color     = styleEqual(names(resistance_colours),
                                         c("white","white","white","white","white")),
                  fontWeight = "bold") %>%
      formatStyle("Detail", cursor = "pointer")
  })

  # Export
  output$export_csv <- downloadHandler(
    filename = function() paste0("VKORC1_filtered_", Sys.Date(), ".csv"),
    content  = function(file) write_csv(filtered_data(), file)
  )
}


# ── 6. HELPERS ────────────────────────────────────────────────────────────────

# Split a site's raw Mutations string into individual mutation entries.
# Handles comma/semicolon/slash separated lists and " and "-joined lists.
# Returns character(0) for blank / "No Mutations" cells.
# (kept for backwards compatibility / export use; the table now uses
#  classify_mutations_detailed() above, which also assigns a level to
#  each entry and keeps known combo mutations together.)
split_mutations <- function(x) {
  x <- as.character(x)
  if (is.na(x) || trimws(x) == "" || grepl("no mutations", x, ignore.case = TRUE)) {
    return(character(0))
  }
  parts <- str_split(x, "[,;/]| and ")[[1]]
  parts <- trimws(parts)
  parts[parts != ""]
}

# ── 7. RUN ────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
