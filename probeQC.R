# Check if required base R packages are installed and install those that are not
pkgs_r <- c('BiocManager', 'gridExtra', 'tidyverse')
pkgs_bioc <- c('GenomicRanges', 'GenomicFeatures', 'BSgenome.Dmelanogaster.UCSC.dm6',
               'rtracklayer', 'ggbio')
pkgs_r_load <- vapply(pkgs_r, library, logical(1), logical.return = TRUE, character.only = TRUE)
print(pkgs_r_load)
if (!all(pkgs_r_load)) install.packages(pkgs_r[pkgs_r_load], Ncpus = 4)
# Do the same for Bioconductor packages
pkgs_bioc_load <- vapply(pkgs_bioc, library, logical(1), logical.return = TRUE, character.only = TRUE)
print(pkgs_bioc_load)
if (!all(pkgs_bioc_load)) BiocManager::install(pkgs_bioc[pkgs_bioc_load])

# Remove 'chr' from dm6 chromosome names
seqlevels(BSgenome.Dmelanogaster.UCSC.dm6) <- sub("chr", "",
                                                  seqlevels(BSgenome.Dmelanogaster.UCSC.dm6))
# Disable scientific notation for numbers
options(scipen=9999999999999999)
# Read the table with all the probes generated in PaintSHOP 
probes_all <- read_table("data/PaintSHOP-full-probe-file.txt",
                         col_names = T)
# Check how many probes are there for each target region
table(probes_all$target)
# Function for filtering multimapping probes. Threshold = # of probe secondary
# mappings to loci ouside of target region
filter_probes <- function(probe_table, thr = 1) {
  # Loop over each unique target in the probe table
  lapply(unique(probe_table$target), function(pset){
    # Filter the probe file
    probeset <- probe_table %>% filter(target == pset)
    # Extract the target information and make the GRanges object
    target <- str_split(unique(probeset$target), c("_|-"))[[1]]
    tagret <- GRanges(sub("chr", "", target[1]), IRanges(start = as.integer(target[2]),
                                                         end = as.integer(target[3])))
    # Create the DNAStringSet object
    probes_fa <- DNAStringSet(probeset$sequence)
    names(probes_fa) <- paste(probeset$chrom, probeset$start, sep = "_")
    # Write the sequences to a fasta file
    writeXStringSet(unique(probes_fa),
                    paste0("fasta/probes_",
                           pset,
                           ".fasta"))
    # Run the script for aligning the fasta onto dm6.48 genome
    system(paste0("mamba run -n align ./map_probes_dm6.sh ",
                  "fasta/probes_", pset, ".fasta"))
    # Load resulting bed file containing each probe mappings
    probe_bed <- import.bed(paste0("probes_",
                                   pset,
                                   "_dm6.bed"))
    # Identify probes mapping outside of the target region
    probes_outside <- subsetByOverlaps(probe_bed, tagret, ignore.strand = T,
                                       invert = T)
    # Identify 'bad' probes mapping more than 'thr' times outside the target region
    bad_probes <- names(table(probes_outside$name)[table(probes_outside$name) > thr])
    # Filter out the bad probes
    probe_f <- probe_bed[!(probe_bed$name %in% bad_probes)]
    # Keep only major chromosomes in probe file
    seqlevels(probe_f, pruning.mode='coarse') <- c("2L",
                                                   "2R",
                                                   "3L",
                                                   "3R",
                                                   "4",
                                                   "X",
                                                   "Y")
    # Set the chromosome information for the filtered probe GRanges
    seqinfo(probe_f) <- seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)[c("2L",
                                                                   "2R",
                                                                   "3L",
                                                                   "3R",
                                                                   "4",
                                                                   "X",
                                                                   "Y")]
    # Add metadata containing target region info
    probe_f@metadata$pset <- pset
    # Return filtered out probes
    probe_f
  })
}
# Filter out probes
probes_f <- filter_probes(probes_all, thr = 100000000)
# Make karyogram plots showing probe distribution on all the major chroms
karyo_plots <- lapply(probes_f, function(probes){
  probe_f2 <- probes
  pset <- probes@metadata$pset
  # Add some additional info to probe GRanges seqlevels so that it could
  # be displayed on the plot
  seqlevels(probe_f2) <- paste(seqlevels(probe_f2),
                               "\n",
                               table(seqnames(probe_f2)),
                               "probes mapped,\n",
                               sapply(seqlevels(probe_f2), 
                                      function(chr) sum(width(GenomicRanges::reduce(probe_f2[seqnames(probe_f2) == chr])))),
                               "bp")
  # Generate a karyogram plot using 'ggbio' package
  p <- autoplot(seqinfo(probe_f2), layout = 'karyogram')+
    layout_karyogram(probe_f2, geom = 'rect', alpha = 0.01)+
    theme(legend.position = "none",
          text=element_text(size=8))+
    # Set target region as the plot title
    ggtitle(pset)
  p@ggplot
})
# Save karyogram plots for all the target regions as a multipage pdf
ggplot2::ggsave("plots/all_probes_karyo.pdf",
                marrangeGrob(karyo_plots,
                             nrow = 3, ncol = 1, top = NULL),
                width = 6, height = 10)
# Number of probes for each region
sapply(probes_f, function(gr) length(unique(gr$name)))

# Make plots displaying probe distributions across the target sites over 1kb bins
bin1kcov <- lapply(probes_f, function(probe){
  # Prepare a GRanges object containing probe set target coordinates
  target <- str_split(unique(probe@metadata$pset), c("_|-"))[[1]]
  tagret <- GRanges(sub("chr", "", target[1]), IRanges(start = as.integer(target[2]),
                                                       end = as.integer(target[3])))
  # print(seqlevels(tagret))
  # Divide the chromosome with the target regions onto 1 kbp bins
  bins <- tileGenome(seqlengths = seqlengths(BSgenome.Dmelanogaster.UCSC.dm6),
                     tilewidth = 1000, cut.last.tile.in.chrom = T)
  bins <- keepSeqlevels(bins, c("2L",
                                "2R",
                                "3L",
                                "3R",
                                "4",
                                "X",
                                "Y"), pruning.mode = 'coarse')
  # Set the strand as '+' because it doesn't matter to this particular application
  strand(probe) <- "+"
  # Find overlaps between 1kbp bins and probes
  ovs <- findOverlaps(bins, probe, ignore.strand = T)
  # Find number of probes falling for each bin
  bins$nprobes = 0
  bins$nprobes[unique(ovs@from)] <- table(ovs@from)
  # leave only bins that have at least one probe mapping to them
  p <- autoplot(bins,
                geom = "col", aes(y=nprobes), col = "black", main = tagret)+
    facet_wrap(~seqnames, ncol = 1, scales = "free_x")+
    ylab("Number of probes per 1 kb")+
    theme_bw()
  p@ggplot
  # Create a DataTrack from bins using number of probes as the y-variable
  
})

ggplot2::ggsave("plots/all_probes_counts_barps.pdf",
                marrangeGrob(bin1kcov,
                             nrow = 1, ncol = 1, top = NULL),
                width = 6, height = 10)

