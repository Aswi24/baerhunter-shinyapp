#' summary
#'
#' @description summary function
#'
#' @return A list of 4 summary dataframes
#' @import rtracklayer
#' @import dplyr
#' @import IRanges
#'
#' @param baerhunter_annot  GFF3 genome annotation file produced by baerhunter
#' @param original_annot The available original annotation file
#' @param refseq_name Name of the reference sequence in the annotation file
#' @return A dataframe with detected sRNAs with the nearby genes as per naming conventions
#'
#'

# function to name annotate sRNAs according to Lamichhane et al. 2012
annot_sRNAs <- function(baerhunter_annot, original_annot,refseq_name){

  ref_granges <- import(original_annot)
  # filter for type='gene'
  ref_granges <- ref_granges[mcols(ref_granges)$type=="gene"]
  ref_genes   <- as.data.frame(ref_granges)

  # isolate only those rows with baerhunter sRNA annotations
  gff <- read.delim(baerhunter_annot, header = FALSE, comment.char = "#")
  df_sRNA <- gff[grepl("putative_sRNA", gff$V3), ]
  df_sRNA <- df_sRNA[, c("V3", "V4", "V5", "V7")]

  # rename columns
  names(df_sRNA) <- c("pred_sRNA", "start", "stop", "strand")
  df_sRNA$start <- as.numeric(df_sRNA$start)
  df_sRNA$stop <- as.numeric(df_sRNA$stop)

  fwd_SRNA_df <- df_sRNA[df_sRNA$strand=="+",]
  rev_SRNA_df <- df_sRNA[df_sRNA$strand=="-",]

  srna_fwd_Granges<-GenomicRanges::GRanges(seqnames = refseq_name,
                            ranges   = IRanges(fwd_SRNA_df$start,
                                               end = fwd_SRNA_df$stop),
                            strand   = S4Vectors::Rle(strand(fwd_SRNA_df$strand)))
  pos_df<-as.data.frame(srna_fwd_Granges)
  srna_rev_Granges<-GenomicRanges::GRanges(seqnames = refseq_name,
                            ranges   = IRanges(rev_SRNA_df$start,
                                               end = rev_SRNA_df$stop),
                            strand   = S4Vectors::Rle(strand(rev_SRNA_df$strand)))
  rev_df<-as.data.frame(srna_rev_Granges)


  if (nrow(pos_df) != 0){
    pos_Granges <- GenomicRanges::GRanges(seqnames = refseq_name,
                           ranges   = IRanges(pos_df$start,
                                              end = pos_df$stop),
                           strand   = S4Vectors::Rle(strand(pos_df$strand))
    )
    ## does srna overlap any gene?
    #ignore strand=T, because want overlap on either strand
    # if srna is on positive strand:
    for (i in 1:nrow(pos_df)){
      fwd_overlaps <- findOverlaps(
        pos_Granges[i],
        ref_granges,
        type = "any",
        minoverlap = 1,
        maxgap = -1,
        select = "last",  #want most downstream gene overlap for naming
        ignore.strand = TRUE
      )
      overlap_gene<-ref_genes$ID[fwd_overlaps]
      pos_df$ov_orf[i]<-overlap_gene
      if (is.na(overlap_gene)==FALSE){
        overlap_name <- overlap_gene
        pos_df$srna_name[i] <- paste("Overlapping", overlap_name, sep=" ")
      }else{
        ## use follow for the genes not already named
        # follow is locus that 'is followed' by sRNA (comes before sRNA==upstream)
        ug<-GenomicRanges::follow(pos_Granges[i], ref_granges, ignore.strand=T)
        follow_locus <- ref_genes$ID[ug]
        follow_name <- follow_locus
        pos_df$srna_name[i] <- paste("Upstream", follow_name, sep=" ")
      }
    }


    # if srna is on negative strand
    if (nrow(rev_df) != 0){
      rev_Granges <- GenomicRanges::GRanges(seqnames = refseq_name,
                             ranges = IRanges(rev_df$start,
                                              end = rev_df$stop),
                             strand = S4Vectors::Rle(strand(rev_df$strand))
      )

      for (i in 1:nrow(rev_df)){
        rev_overlaps <- findOverlaps(
          rev_Granges[i],
          ref_granges,
          type = "any",
          minoverlap = 1,
          maxgap = -1,
          select = "first",
          ignore.strand = TRUE
        )
        overlap_gene<-ref_genes$ID[rev_overlaps]
        rev_df$ov_orf[i]<-overlap_gene
        if (is.na(overlap_gene)==FALSE){
          overlap_name <- overlap_gene
          rev_df$srna_name[i] <- paste("Overlapping", overlap_name, sep=" ")
        }else{
          ## use follow for the genes not already named
          # follow is gene that 'is followed' by sRNA (comes before sRNA==upstream)
          # this is reversed for negative strand: use precede
          ug<-GenomicRanges::precede(rev_Granges[i], ref_granges, ignore.strand=T)
          follow_locus <- ref_genes$ID[ug]
          follow_name <- follow_locus
          rev_df$srna_name[i] <- paste("Upstream", follow_name, sep=" ")
        }
      }
    }

    # join back together
    srna_df <- rbind(pos_df, rev_df)
    srna_df <- srna_df[order(srna_df$start),]
    srna_df <- srna_df[, c(1, 2, 3,5,7)]
    new_names <- c("ncRNA", "start", "stop", "strand", "type")
    colnames(srna_df)<- new_names
    srna_df[, 1] <- "putative_sRNA"
    srna_df
    return(srna_df)
  }
}

#' @param baerhunter_annot  GFF3 genome annotation file produced by baerhunter
#' @param original_annot The available original annotation file
#' @param refseq_name Name of the reference sequence in the annotation file
#' @return A dataframe with detected UTRs with the nearby genes as per naming conventions

annot_utrs <- function(baerhunter_annot, original_annot,refseq_name){

  ref_granges <- import(original_annot)
  # filter for type='gene'
  ref_granges <- ref_granges[mcols(ref_granges)$type=="gene"]
  ref_genes   <- as.data.frame(ref_granges)

  #assign utrs
  gff <- read.delim(baerhunter_annot, header = FALSE, comment.char = "#")
  UTR_df <- gff[grepl("putative_UTR", gff$V3), ]
  UTR_df <- UTR_df[, c("V3", "V4", "V5", "V7")]

  # rename columns
  names(UTR_df) <- c("pred_UTR", "start", "stop", "strand")
  UTR_df$start <- as.numeric(UTR_df$start)
  UTR_df$stop <- as.numeric(UTR_df$stop)



  # create grange object of utr coordinates
  utr_Granges<-GenomicRanges::GRanges(seqnames = refseq_name,
                       ranges   = IRanges(UTR_df$start,
                                          end = UTR_df$stop),
                       strand   = S4Vectors::Rle(strand(UTR_df$strand))
  )

  # searching whether utr co-locates with gene?
  # precede is gene that "is preceded" by the UTR (comes after UTR==downstream)
  pg<-GenomicRanges::precede(utr_Granges, ref_granges,ignore.strand=F)
  utr_precede_genes<-ref_genes$ID[pg]
  # precede doesn't pay attention to circular genome, use first gene in gff
  utr_precede_genes[is.na(utr_precede_genes)] <- ref_genes$ID[1]
  # follow is gene that 'is followed' by UTR (comes before UTR==upstream)
  ug<-GenomicRanges::follow(utr_Granges, ref_granges, ignore.strand=F)
  utr_follow_genes<-ref_genes$ID[ug]
  utr_follow_genes[is.na(utr_follow_genes)] <- ref_genes$ID[nrow(ref_genes)]
  ## nearest(utrs, refseq, select="all", ignore.strand = F) will give me nearest gene feature.
  # select=all means returns both nearest features in case of a tie, arbitrary chooses one.
  # see how many ties
  a<-nearest(utr_Granges, ref_granges, select="all", ignore.strand=F)
  b<-as.matrix(a)
  #show frequency of hits for each utr (1 hit or 2 hits, 2=tie)
  c<-as.data.frame(table(b[,1]))
  utr_tied<-c$Freq
  # see nearest gene (arbitrarily choose in ties)
  ng<-nearest(utr_Granges, ref_granges, select="arbitrary", ignore.strand=F)
  utr_nearest_genes<-ref_genes$ID[ng]

  #include preceding (downstream gene), following (upstream gene), and nearest genes in df
  # (and indicate if there is a tie)
  UTR_df$downstream <-utr_precede_genes
  UTR_df$upstream   <-utr_follow_genes
  UTR_df$nearest    <-utr_nearest_genes
  UTR_df$tied       <-utr_tied




  # subset utr_granges by strand
  utr_fwd_Granges<-utr_Granges[strand(utr_Granges) == "+"]
  fwd_df<-as.data.frame(utr_fwd_Granges)
  utr_rev_Granges<-utr_Granges[strand(utr_Granges) =="-"]
  rev_df<-as.data.frame(utr_rev_Granges)
  utr_df<-rbind(fwd_df, rev_df)

  for (i in 1:nrow(UTR_df)){
    if (UTR_df$strand[i] == "+"){
      UTR_df$downstream_start[i]  <-
        ref_genes$start[match(UTR_df$downstream[i], ref_genes$ID)]
      # consider edge case of final utr on strand
      UTR_df$dist_to_start[i] <- ifelse(UTR_df$downstream_start[i] - UTR_df$stop[i] > 0,
                                        UTR_df$downstream_start[i] - UTR_df$stop[i], NA)
    }else{
      # for minus strand, need end coordinate of downstream gene
      UTR_df$downstream_start[i]  <-
        ref_genes$end[match(UTR_df$downstream[i], ref_genes$ID)]
      UTR_df$dist_to_start[i] <- ifelse(UTR_df$start[i] - UTR_df$downstream_start[i] > 0,
                                        UTR_df$start[i] - UTR_df$downstream_start[i], NA) }
  }
  # assign 5', 3' or btwn to each utr

  for (i in 1:nrow(UTR_df)){
    # if tied (same dist between two genes, located between two genes)
    if (UTR_df$tied[i] == 2){
      UTR_df$utr[i] <- paste("UTR BTWN genes;downstream", UTR_df$downstream[i], sep=" ")
      #not tied (closer to either upstream or downstream gene)
    }else{
      if (UTR_df$nearest[i] == UTR_df$upstream[i]) {
        # check distance from end of UTR to start of downstream gene
        if (is.na(UTR_df$dist_to_start[i]) | UTR_df$dist_to_start[i] > 40){
          UTR_df$utr[i] <- paste("3UTR", UTR_df$upstream[i], sep=" ")
        }else{
          UTR_df$utr[i] <- paste("UTR BTWN genes;downstream", UTR_df$downstream[i], sep=" ") }
      }else{
        # if nearest==downstream, likely 5' UTR
        UTR_df$utr[i] <- paste("5UTR", UTR_df$downstream[i], sep=" ") }
    }
  }
  UTR_df <- UTR_df[, c(1, 2, 3, 4, 11)]
  UTR_df[, 1] <- "putative_UTR"
  new_names <-c("ncRNA", "start", "stop", "strand", "type")
  colnames(UTR_df)<- new_names
  return(UTR_df)
}

#' @param ori_ranges  Start and end ranges of ncRNA genes found in the original annotation file
#' @param bh_ranges Start and end ranges of ncRNA genes found in the annotation file produced by baerhunter
#' @return A dataframe with those ranges that overlap between the 2


find_overlaps <- function (ori_ranges, bh_ranges){
  overlaps <- findOverlaps(ori_ranges, bh_ranges)
  overlapping_ori_range <- ori_ranges[S4Vectors::queryHits(overlaps)]
  overlapping_bh_range <- bh_ranges[S4Vectors::subjectHits(overlaps)]
  overlapping_ori_df <- as.data.frame(overlapping_ori_range)
  overlapping_bh_df <- as.data.frame(overlapping_bh_range)
  # Adding prefixes to column names
  colnames(overlapping_ori_df) <- paste("Original", colnames(overlapping_ori_df), sep="_")
  colnames(overlapping_bh_df) <- paste("BaerHunter", colnames(overlapping_bh_df), sep="_")

  # Combining the data frames side by side
  overlap_df <- cbind(overlapping_bh_df, overlapping_ori_df)
  return( overlap_df)
}

#' @param baerhunter_annot  GFF3 genome annotation file produced by baerhunter
#' @param original_annot The available original annotation file
#' @param refseq_name Name of the reference sequence in the annotation file
#' @return A list of dataframes - a df of UTRs and SRNAs; overlapping annotations detected by baerhunter,
#' summary of counts of the categories of nCRNAs detected; the ncRNAs in original annotation file
#' undetected by baerhunter


ncrna_annot_compare <- function(original_annot, baerhunter_annot,refseq_name) {

  if (!file.exists(original_annot)) {
    print("original annotation file does not exist! Please provide a valid file path")
    stop(paste("File not found:", original_annot))
  }

  if (!file.exists(baerhunter_annot)) {
    print(" baerhunter_annot file does not exist! Please provide a valid file path")
    stop(paste("File not found:", baerhunter_annot))
  }

  ori_gff <- read.delim(original_annot, header = FALSE, comment.char = "#")
  bh_gff <- read.delim(baerhunter_annot, header = FALSE, comment.char = "#")



  sRNA_df <- annot_sRNAs(baerhunter_annot, original_annot,refseq_name)
  UTR_df <- annot_utrs(baerhunter_annot, original_annot,refseq_name)
  sRNA_UTR_df <- dplyr::bind_rows(UTR_df, sRNA_df)


  search_term_ori <- "UTR|[^trm]RNA"
  # select for annotations of UTR and sRNA according to baerhunter's nomenclature
  search_term_bh <- "putative_(UTR|sRNA)"

  ori_ncrna <-  ori_gff [grepl(search_term_ori, ori_gff [,3], ignore.case = TRUE)==TRUE ,]

  # look for ncRNAs in baerhunter annotation
  bh_ncrna <- bh_gff[grepl(search_term_bh, bh_gff[,3], ignore.case = TRUE)==TRUE ,]

  ori_ranges <- IRanges(start = ori_ncrna[,4], end = ori_ncrna[,5], names = ori_ncrna[,3])
  bh_ranges <- IRanges(start =  bh_ncrna[,4], end =  bh_ncrna[,5], names =  bh_ncrna[,3])

  overlap_df <- find_overlaps(ori_ranges,bh_ranges)

  ori_df <- as.data.frame(ori_ranges)
  colnames(ori_df) <- c("Original_start", "Original_end","Original_width","Original_names")

  unique_df <- overlap_df %>%
    mutate(location = paste(Original_start, Original_end, sep = "-")) %>%
    distinct(location, .keep_all = TRUE)

  undetected_annot <- anti_join(ori_df, overlap_df, by = c("Original_start", "Original_end", "Original_names"))
  undetected_annot <- rename(undetected_annot, undetected_ncRNAs = Original_names)


  sRNA_UTR_df$category <- NA


  sRNA_UTR_df$category[grepl('overlapping', tolower(sRNA_UTR_df$type)) & sRNA_UTR_df$strand == '-'] <- "Antisense RNA "
  sRNA_UTR_df$category[grepl('overlapping', tolower(sRNA_UTR_df$type)) & sRNA_UTR_df$strand == '+'] <- "Between UTRs"
  sRNA_UTR_df$category[grepl('upstream', tolower(sRNA_UTR_df$type)) & sRNA_UTR_df$strand == '-'] <- "Intergenic RNA (minus strand)"
  sRNA_UTR_df$category[grepl('upstream', tolower(sRNA_UTR_df$type)) & sRNA_UTR_df$strand == '+'] <- "Intergenic RNA (plus strand)"
  sRNA_UTR_df$category[grepl('^3UTR', sRNA_UTR_df$type, ignore.case = TRUE)] <- "3' UTRs"
  sRNA_UTR_df$category[grepl('^5UTR', sRNA_UTR_df$type, ignore.case = TRUE)] <- "5' UTRs"
  sRNA_UTR_df$category[grepl('BTWN', sRNA_UTR_df$type, ignore.case = TRUE)] <- "Unassigned UTRs"

  # Produce a final counts dataframe
  counts_df <- as.data.frame(table(sRNA_UTR_df$category))
  colnames(counts_df) <- c("Type", "Count")
  # count total original annotations

  total_original_annotations <- nrow(ori_ncrna)
  counts_df <- rbind(counts_df, data.frame(Type = "Total original annotations",
                                           Count = total_original_annotations,
                                           stringsAsFactors = FALSE))

  #count annotations produced by baerhunter
  annotations_produced <- nrow(bh_ncrna)
  counts_df <- rbind(counts_df, data.frame(Type = "Annotations produced by baerhunter",
                                           Count = annotations_produced,
                                           stringsAsFactors = FALSE))

  # count overlaps of baerhunter annotations with the original
  counts_overlaps <- nrow(overlap_df)
  counts_df <- rbind(counts_df, data.frame(Type = "Overlapping annotations",
                                           Count = counts_overlaps,
                                           stringsAsFactors = FALSE))

  num_unique_locations <- nrow(unique_df)
  counts_df <- rbind(counts_df, data.frame(Type = "Overlapping original annotations",
                                           Count = num_unique_locations,
                                           stringsAsFactors = FALSE))




  return(list(sRNA_UTR_df, overlap_df,  counts_df, unique_df,undetected_annot))
}


