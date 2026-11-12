## ============================================================================
## 00 — Expected composition of the mock community by in-silico PCR
##
## The mock is the ZymoBIOMICS Microbial Community DNA Standard D6305, added at
## the PCR step (no extraction). Its passport gives the genomic DNA composition:
## 12% for each of 8 bacteria, 2% for each of 2 yeasts. Equal DNA mass is not an
## equal number of amplicon templates: a smaller genome gives more genome copies
## per ng, and each genome carries several rRNA operons. So the expected share of
## member i among the amplicons of one marker is
##
##   templates_i = DNA_mass_i / genome_size_i  x  amplifiable_copies_i
##   share_i     = templates_i / sum(templates of the marker's members)
##
## (Zymo manual, Table 1 footnote 1: the same formula.) This script finds the
## amplifiable copies by in-silico PCR: it searches the project's primers in
## every genome of the Zymo reference archive, extracts each amplicon, counts
## copies per genome and distinct sequences (intragenomic variants), and applies
## the read-length limit of merged ITS1 reads.
##
## Sources:
##   data/ZymoBIOMICS.STD.refseq.v3.zip   genomes (v3, July 2023)
##   data/ZymoBIOMICS_D6300_instruction_manual.pdf   Table 1 (theoretical
##     composition) and Table 2 (genome size, ploidy, rRNA copies); the D6305
##     DNA standard has the same genomic DNA composition
##   workflow/config_CFM_v2.yaml          primers, read length, filtering
##   workflow/scripts/dada2_denoise.R      mergePairs(minOverlap = 12)
##
## Output (results/<version>/tables/):
##   00_mock_expected_amplicons.csv    one row per distinct amplicon per member
##   00_mock_expected_amplicons.fasta  the amplicon inserts (between the primers)
##   00_mock_expected_composition.csv  expected share per member, marker, read
##                                     set (merged / R1-only) and basis
## Runtime about 1 min (whole genomes, both strands).
## ============================================================================

# ---- 0. Packages ---------------------------------------------------------
library(here)
library(tidyverse)
library(Biostrings)   # DNA sequences: readDNAStringSet(), matchPattern(), reverseComplement()
library(DECIPHER)     # Clusterize(): groups amplicons at 97% identity, like the OTU97 step

# Data version: RESULTS_DIR (the tables go to results/<version>/tables)
source(here("scripts", "utils_paths.R"))

# ---- 1. Parameters -------------------------------------------------------
# Primers: the locus-specific part only (workflow/config_CFM_v2.yaml). The labels
# 515fA, 806r3, ITS1fA, ITS2r4 in the metadata name 8-nt inline tags in front of
# the primer; the tags are removed before primer trimming and do not enter PCR
# specificity, so one in-silico PCR per marker covers all mocks.
PRIMERS <- list(
  "16S"  = c(fwd = "GTGYCAGCMGCCGCGGTAA",    rev = "GGACTACNVGGGTWTCTAAT"),   # 515F / 806R
  "ITS1" = c(fwd = "CTTGGTCATTTAGAGGAAGTAA", rev = "GCTGCGTTCTTCATCGATGC")    # ITS1F / ITS2
)
MAX_MISMATCH <- 3      # mismatches allowed per primer site (see section 4)
MAX_INSERT   <- 1500   # longest insert between the primers that counts as an amplicon

# Merged reads: DADA2 joins R1 and R2 only if they overlap by >= 12 nt, so the
# longest insert that can merge is R1 + R2 - 12. After the 8-nt tag and the
# primer are trimmed, 2 x 250 reads leave R1 about 220 nt (ITS1F, 22 nt) and
# R2 about 222 nt (ITS2, 20 nt; inferred from the read length, not measured).
R1_LEN      <- 220
R2_LEN      <- 222
MIN_OVERLAP <- 12
MAX_MERGED_INSERT <- R1_LEN + R2_LEN - MIN_OVERLAP   # 430 nt
# R1-only ITS1 branch: R1 cut to 220 nt, no merging, so no length limit
R1_ONLY_LEN <- 220

# ---- 2. Passport: DNA composition, genome size, rRNA copies --------------
# DNA %: theoretical genomic DNA composition (manual Table 1).
# Genome size (Mb, haploid) and 16S/18S copies per haploid genome: manual Table 2.
# The yeast copy numbers were estimated by Zymo from shotgun read depth.
# Family_Genus: the label each member gets in our taxonomy (as in script 03),
# for SILVA 138 (v1) and SILVA 144 / UNITE (v2).
passport <- tribble(
  ~member,                     ~species,                    ~group,     ~dna_pct, ~table_mb, ~table_copies, ~label_v1,                                         ~label_v2,
  "Pseudomonas_aeruginosa",    "Pseudomonas aeruginosa",    "bacteria", 12, 6.792,   4, "f__Pseudomonadaceae | g__Pseudomonas",            "f__Pseudomonadaceae | g__Pseudomonas",
  "Escherichia_coli",          "Escherichia coli",          "bacteria", 12, 4.875,   7, "f__Enterobacteriaceae | g__Escherichia-Shigella", "f__Enterobacteriaceae | g__Escherichia-Shigella",
  "Salmonella_enterica",       "Salmonella enterica",       "bacteria", 12, 4.760,   7, "f__Enterobacteriaceae | g__Unassigned",           "f__Enterobacteriaceae | g__Unassigned",
  "Lactobacillus_fermentum",   "Lactobacillus fermentum",   "bacteria", 12, 1.905,   5, "f__Lactobacillaceae | g__Lactobacillus",          "f__Lactobacillaceae | g__Limosilactobacillus",
  "Enterococcus_faecalis",     "Enterococcus faecalis",     "bacteria", 12, 2.845,   4, "f__Enterococcaceae | g__Enterococcus",            "f__Streptococcaceae | g__Enterococcus",
  "Staphylococcus_aureus",     "Staphylococcus aureus",     "bacteria", 12, 2.730,   6, "f__Staphylococcaceae | g__Staphylococcus",        "f__Staphylococcaceae | g__Staphylococcus",
  "Listeria_monocytogenes",    "Listeria monocytogenes",    "bacteria", 12, 2.992,   6, "f__Listeriaceae | g__Listeria",                   "f__Listeriaceae | g__Listeria",
  "Bacillus_subtilis",         "Bacillus subtilis",         "bacteria", 12, 4.045,  10, "f__Bacillaceae | g__Bacillus",                    "f__Bacillaceae | g__Bacillus",
  "Saccharomyces_cerevisiae",  "Saccharomyces cerevisiae",  "yeast",     2, 12.1,  109, "f__Saccharomycetaceae | g__Saccharomyces",        "f__Saccharomycetaceae | g__Saccharomyces",
  "Cryptococcus_neoformans",   "Cryptococcus neoformans",   "yeast",     2, 18.9,   60, "f__Cryptococcaceae | g__Cryptococcus",            "f__Cryptococcaceae | g__Cryptococcus"
)
# Which members each marker targets: 16S the bacteria, ITS1 the yeasts
TARGETS <- c("16S" = "bacteria", "ITS1" = "yeast")

# Manual Table 1, "16S only" column: used below to check the formula
TABLE1_16S <- c(Pseudomonas_aeruginosa = 4.2, Escherichia_coli = 10.1, Salmonella_enterica = 10.4,
                Lactobacillus_fermentum = 18.4, Enterococcus_faecalis = 9.9, Staphylococcus_aureus = 15.5,
                Listeria_monocytogenes = 14.1, Bacillus_subtilis = 17.4)

# ---- 3. Read the genomes from the Zymo archive ---------------------------
# Each genome is a FASTA file of contigs: one per chromosome or plasmid for the
# complete bacterial genomes, many short ones for the draft yeast assemblies.
# The archive is unpacked into a temporary folder, never into the repository.
# (The archive also holds 16S/18S files from older assemblies; they are not
# used, because the ITS1 region has to come from the genome anyway.)
zip_file  <- here("data", "ZymoBIOMICS.STD.refseq.v3.zip")
in_zip    <- unzip(zip_file, list = TRUE)$Name
genome_in <- str_subset(in_zip, "/Genomes/.+\\.(fna|fasta)$")
unzip(zip_file, files = genome_in, exdir = tempdir())

genomes <- tibble(path = file.path(tempdir(), genome_in)) |>
  mutate(member = str_remove(basename(path), "_complete_genome\\..*$"),
         dna    = map(path, \(p) readDNAStringSet(p) |> setNames(nm = NULL)),
         assembly_bp = map_dbl(dna, \(g) sum(width(g))))
stopifnot(setequal(genomes$member, passport$member))

# ---- 4. In-silico PCR ----------------------------------------------------
# For one primer and one strand of every contig, find the binding sites.
# matchPattern(..., fixed = "subject"): the IUPAC codes of the primer (Y, M, N,
# V, W) match any base they stand for, while an N in the genome does not match
# everything. max.mismatch allows up to MAX_MISMATCH mismatches.
# For each site we also record where the mismatch nearest the primer's 3' end
# sits (counted from that end; 0 = no mismatch), because 3'-end mismatches
# block polymerase extension most strongly.
#   three_prime_right = TRUE for the forward primer (its 3' end is on the right);
#   the reverse primer is searched as its reverse complement, whose 3' end is on
#   the left.
find_sites <- function(contigs, primer, three_prime_right) {
  primer <- DNAString(primer)
  imap(as.list(contigs), \(contig, i) {
    hits <- matchPattern(primer, contig, max.mismatch = MAX_MISMATCH, fixed = "subject")
    if (length(hits) == 0) return(NULL)
    # mismatch() lists, for each hit, the positions in the primer that do not match
    mm_pos <- mismatch(primer, hits, fixed = "subject")
    tibble(contig = i, start = start(hits), end = end(hits),
           mismatches = lengths(mm_pos),
           mm_from_3prime = map_int(mm_pos, \(p) {
             if (length(p) == 0) return(0L)
             as.integer(min(if (three_prime_right) length(primer) - p + 1 else p))
           }))
  }) |>
    list_rbind() |>
    (\(res) if (nrow(res) == 0) NULL else res)()   # no site on this strand: NULL, not an empty table
}

# Amplicons of one marker in one genome. The genome is scanned on both strands:
# on each strand the forward primer binds as written and the reverse primer as
# its reverse complement, downstream of it. Each forward site is paired with
# the nearest reverse site within MAX_INSERT nt; the insert is the sequence
# between the two primers, which is what the reads hold after primer trimming.
in_silico_pcr <- function(contigs, marker) {
  fwd    <- PRIMERS[[marker]][["fwd"]]
  rev_rc <- as.character(reverseComplement(DNAString(PRIMERS[[marker]][["rev"]])))
  map(c("plus", "minus"), \(strand) {
    s <- if (strand == "plus") contigs else reverseComplement(contigs)
    f <- find_sites(s, fwd,    three_prime_right = TRUE)
    r <- find_sites(s, rev_rc, three_prime_right = FALSE)
    if (is.null(f) || is.null(r)) return(NULL)
    inner_join(f, r, by = "contig", suffix = c("_fwd", "_rev"), relationship = "many-to-many") |>
      filter(start_rev > end_fwd, start_rev - end_fwd - 1 <= MAX_INSERT) |>
      slice_min(start_rev, n = 1, by = c(contig, start_fwd)) |>
      mutate(strand = strand,
             insert = as.character(subseq(s[contig], start = end_fwd + 1, end = start_rev - 1)))
  }) |> list_rbind()
}

# Every member x marker. A copy counts as amplifiable when both primer sites have
# at most MAX_MISMATCH mismatches, with no extra rule for the 3' end: a strict
# 3'-end rule would call Cryptococcus non-amplifiable (all its ITS2 sites have
# one mismatch at the 3'-terminal base), yet it gives 85-99.9% of the ITS1 mock
# reads. The mismatches are kept in the table as possible causes of bias.
amplicons <- genomes |>
  select(member, dna) |>
  expand_grid(marker = names(PRIMERS)) |>
  mutate(hits = map2(dna, marker, in_silico_pcr)) |>
  select(member, marker, hits) |>
  unnest(hits) |>
  left_join(select(passport, member, group), by = "member") |>
  mutate(off_target = group != TARGETS[marker])   # e.g. 16S primers on yeast 18S

# ---- 5. Copies and intragenomic variants ---------------------------------
# A genome can carry several rRNA operons whose amplicons are identical or
# differ by a few bases. Here identical inserts are counted together: copies =
# how many operons give that sequence. The share of each variant within its
# member is copies / all copies of the member. For the yeasts the rDNA array
# is collapsed in the draft assembly, so these variant proportions are rough.
variants <- amplicons |>
  summarise(copies = n(),
            mismatches_fwd = min(mismatches_fwd), mm_3prime_fwd = min(mm_from_3prime_fwd),
            mismatches_rev = min(mismatches_rev), mm_3prime_rev = min(mm_from_3prime_rev),
            .by = c(member, group, marker, off_target, insert)) |>
  mutate(insert_len   = nchar(insert),
         amplicon_len = insert_len + nchar(PRIMERS[marker] |> map_chr("fwd")) +
                                     nchar(PRIMERS[marker] |> map_chr("rev")),
         # Merged reads exist only for inserts the read pair can overlap (ITS1);
         # the 253-nt 16S insert always merges
         mergeable    = insert_len <= MAX_MERGED_INSERT,
         r1_prefix    = str_sub(insert, 1, R1_ONLY_LEN)) |>
  arrange(marker, off_target, member, desc(copies)) |>
  mutate(variant = row_number(), copies_member = sum(copies),
         variant_share_in_member = copies / copies_member, .by = c(member, marker))

# ---- 6. OTU97 check ------------------------------------------------------
# The ITS1 data are clustered into OTUs at 97% identity (vsearch). Clusterize()
# does the same kind of grouping here (its identity definition differs slightly),
# to check that the intragenomic variants of one yeast fall into one OTU, both
# for the full insert (merged reads) and for the 220-nt R1 prefix (R1-only).
its1 <- variants |> filter(marker == "ITS1", !off_target)
otu97 <- function(seqs) {
  Clusterize(DNAStringSet(setNames(seqs, paste(its1$member, its1$variant))),
             cutoff = 0.03, processors = 1, verbose = FALSE)$cluster
}
its1 <- its1 |> mutate(otu97_full = otu97(insert), otu97_r1 = otu97(r1_prefix))
variants <- variants |>
  left_join(select(its1, member, marker, variant, otu97_full, otu97_r1), by = c("member", "marker", "variant"))
cat("ITS1 OTU97 per member (full insert / R1 prefix):\n")
print(its1 |> summarise(variants = n(), otu97_full = n_distinct(otu97_full),
                        otu97_r1 = n_distinct(otu97_r1), .by = member))

# ---- 7. Expected shares --------------------------------------------------
# templates_i = DNA % / genome size (bp) x amplifiable copies; the unit
# conversion (bp per gram) is the same for every member and cancels.
# Three bases, written side by side so the effect of each choice is visible:
#   default        bacteria: genome size and copies from the archive genomes
#                  (complete assemblies); yeasts: manual Table 2, because the
#                  draft yeast assemblies are larger than the haploid genome
#                  (duplicate haplotype contigs) and collapse the rDNA array
#   zymo_table     manual Table 2 for every member (reproduces Table 1)
#   assembly_only  archive genomes for every member, yeasts included
# Read sets: merged reads (all markers) and R1-only (ITS1). For merged ITS1,
# only the copies with a mergeable insert count.
copies_member <- variants |>
  filter(!off_target) |>
  summarise(copies_insilico   = sum(copies),
            mergeable_fraction = sum(copies[mergeable]) / sum(copies),
            .by = c(member, marker))

composition <- passport |>
  inner_join(copies_member, by = "member") |>
  left_join(select(genomes, member, assembly_bp), by = "member") |>
  expand_grid(basis = c("default", "zymo_table", "assembly_only")) |>
  mutate(
    use_table   = basis == "zymo_table" | (basis == "default" & group == "yeast"),
    genome_bp   = if_else(use_table, table_mb * 1e6, assembly_bp),
    copies      = if_else(use_table, table_copies, copies_insilico),
    source      = if_else(use_table, "manual Table 2", "archive genome, in-silico PCR")
  ) |>
  expand_grid(read_set = c("merged", "R1_only")) |>
  filter(marker == "ITS1" | read_set == "merged") |>
  mutate(templates = dna_pct / genome_bp * copies *
                     if_else(read_set == "merged", mergeable_fraction, 1)) |>
  mutate(expected_share = templates / sum(templates), .by = c(marker, read_set, basis)) |>
  select(marker, read_set, basis, member, species, dna_pct, genome_bp, copies, source,
         mergeable_fraction, templates, expected_share, label_v1, label_v2) |>
  arrange(marker, read_set, basis, desc(expected_share))

# Check: the zymo_table basis must reproduce the manual's Table 1 "16S only"
check <- composition |> filter(marker == "16S", basis == "zymo_table")
stopifnot(all(abs(100 * check$expected_share - TABLE1_16S[check$member]) < 0.1))

# ---- 8. Save -------------------------------------------------------------
dir.create(file.path(RESULTS_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
variants |>
  select(member, group, marker, off_target, variant, copies, copies_member, variant_share_in_member,
         insert_len, amplicon_len, mergeable, mismatches_fwd, mm_3prime_fwd, mismatches_rev,
         mm_3prime_rev, otu97_full, otu97_r1, insert, r1_prefix) |>
  write_csv(file.path(RESULTS_DIR, "tables", "00_mock_expected_amplicons.csv"))
variants |>
  (\(v) DNAStringSet(setNames(v$insert, paste(v$member, v$marker, paste0("variant", v$variant),
                                               paste0("copies", v$copies), sep = "|"))))() |>
  writeXStringSet(file.path(RESULTS_DIR, "tables", "00_mock_expected_amplicons.fasta"))
write_csv(composition, file.path(RESULTS_DIR, "tables", "00_mock_expected_composition.csv"))

composition |>
  filter(basis == "default") |>
  transmute(marker, read_set, member, expected_pct = round(100 * expected_share, 1)) |>
  print(n = Inf)
cat("00 done.\n")
