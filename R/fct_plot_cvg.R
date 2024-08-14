#' plot_cvg
#'
#' @description A  function get coverage plot of intergenic regions
#'
#' @return The return value, if any, from executing the function.
#'
#' @noRd
#'

#' @param gff_file_path Filepath of a gff annotation file
#' @param ori_sRNA_biotype A string indicating how the biotype of pre-annotated ncRNA, which can be found in the attribute column.In case if the user does not know how the sRNA is annotated, it can be set as "unknown". In this case, all RNAs apart from tRNAs and rRNAs will be removed from the selection.
#' @return A list of GRanges object containing the intergenic coordinates for the plus and minus strand

igr_strand_specific <- function(gff_file_path, ori_sRNA_biotype){

  # the major features dataframe is constructed the same way in baerhunter
  gff <- read.delim(gff_file_path, header = FALSE, comment.char = "#")
  major_f <- gff[grepl("Parent", gff[,9], ignore.case = TRUE)==FALSE & gff[,3]!='chromosome' & gff[,3]!='biological_region' & grepl(ori_sRNA_biotype, gff[,9], ignore.case = TRUE)==FALSE & gff[,3]!='region' & gff[,3]!='sequence_feature',]
  colnames(major_f ) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  major_f$phase <- as.integer(major_f$phase)
  major_f$seqid <- rep("AL123456.3", nrow(major_f))
  granges_major_f <-GenomicRanges::makeGRangesFromDataFrame(major_f ,
                                             keep.extra.columns=TRUE,
                                             ignore.strand=FALSE,
                                             seqinfo=NULL,
                                             seqnames.field=c("seqnames", "seqname",
                                                              "chromosome", "chrom",
                                                              "chr", "chromosome_name",
                                                              "seqid"),
                                             start.field="start",
                                             end.field=c("end", "stop"),
                                             strand.field="strand",
  )
  # calculate gaps in the genome not covered by the features
  gaps_ranges <- gaps(granges_major_f)

  #get the gap ranges on the plus and minus  strand
  gaps_plus <- gaps_ranges[strand(gaps_ranges) == "+"]
  gaps_minus <- gaps_ranges[strand(gaps_ranges) == "-"]

  return(list(gaps_plus = gaps_plus, gaps_minus = gaps_minus))
}

#' @param gff_file_path An rle list containing coverage
#' @return A coverage dataframe

expand_Rle_df <- function(rle_list) {

  #expands the list to get coverage values for each individual position in the gap ranges
  coverage_vectors <- lapply(rle_list, as.integer)

  positions <- unlist(lapply( coverage_vectors, seq_along))

  # Combine coverage values into a single vector
  coverage_values <- unlist(coverage_vectors )

  # Create dataframe
  coverage_df <- data.frame(Position = positions, Coverage = coverage_values)

  return(coverage_df)
}

#' @param bam_file File path of a bam file
#' @param bam_name Name of the bam file
#' @param gaps_plus A Granges onject containing the intergenic coordinates for the plus strand
#' @param gaps_minus A Granges onject containing the intergenic coordinates for the minus strand
#' @param paired_end_data A boolean indicating if the reads are paired-end.
#' @param strandedness A string outlining the type of the sequencing library: stranded, or reversely stranded
#' @return a list of coverage dataframe and percentiles dataframe

process_bam_file <- function(bam_file, bam_name, gaps_plus, gaps_minus, paired_end_data = FALSE, strandedness  = "unstranded" ) {


  if (paired_end_data == FALSE) {
    aln <- readGAlignments(bam_file)
    if (strandedness == "reversely_stranded") {
      aln <- aln[strand(aln) == ifelse(strand(aln) == "+", "-", "+")]
    }
  } else {
    if (strandedness == "stranded") {
      aln <- readGAlignmentPairs(bam_file, strandMode = 1)
    } else if (strandedness == "reversely_stranded") {
      aln <- readGAlignmentPairs(bam_file, strandMode = 2)
    } else {
      aln <- readGAlignmentPairs(bam_file)
    }
  }

  pcvg <- coverage(aln[strand(aln) == "+"])
  pcvg_igr_p <- pcvg[gaps_plus]

  mcvg <- coverage(aln[strand(aln) == "-"])
  mcvg_igr_m <- mcvg[gaps_minus]

  m_cvg_df <- expand_Rle_df(mcvg_igr_m)
  m_cvg_df$Strand <- "-"

  p_cvg_df <- expand_Rle_df(pcvg_igr_p)
  p_cvg_df$Strand <- "+"

  coverage_df<- rbind(m_cvg_df, p_cvg_df)
  coverage_df$File <- bam_name

  #remove all those with zero coverage
  coverage_df <- coverage_df %>% filter(Coverage != 0)

  #get the percentiles
  percentiles <- quantile(coverage_df$Coverage, probs = c(0.25, 0.5, 0.75, 0.8, 0.85, 0.9), names = TRUE)
  percentiles_df <- data.frame(
    File = bam_name,
    TwentyFifth = percentiles['25%'],
    Fifty = percentiles['50%'],
    SeventyFifth = percentiles['75%'],
    Eightieth = percentiles['80%'],
    EightyFifth = percentiles['85%'],
    Ninetieth = percentiles['90%'],
    row.names = NULL
  )

  colnames(percentiles_df) <- c("File", "P25", "Median", "P75", "P80", "P85", "P90")


  return(list(coverage_df = coverage_df, percentiles_df = percentiles_df))
}

