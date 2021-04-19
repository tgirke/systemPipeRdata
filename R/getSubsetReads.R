############################
## Subsetting fastq data ##
############################
getSubsetReads <- function(args,
                           geneList = NULL, gr = NULL, MappingRegion = 1:100000,
                           sample_range = 90000:100000,
                           truncate_refs = TRUE,
                           id_read_number = TRUE,
                           annotation = "data/tair10.gff",
                           reference = "data/tair10.fasta",
                           annot_outname = "tair10_sub.gff",
                           ref_outname = "tair10_sub.fasta",
                           outdir = "data/subset/", silent = FALSE) {
  ## Validation
  if (all(class(args) != "SYSargs" & class(args) != "SYSargs2")) stop("Argument 'args' needs to be assigned an object of class 'SYSargs' OR 'SYSargs2")
  if (!file.exists(annotation)) stop("Please provide a valid annotation path.")
  if (!file.exists(reference)) stop("Please provide a valid reference genome path.")
  ## check and create output directory
  if (!dir.exists(outdir)) {
    dir.create(outdir)
    if (!(silent)) cat(paste0("Creating '", outdir, "' directory"), "\n")
  }
  ## obtain bam PATH
  outpaths <- systemPipeR::subsetWF(args, slot = "output", subset = 1, index = 1)
  ## Create `txdb` for annontation
  if (!is.null(geneList)) {
    suppressWarnings(txdb <- GenomicFeatures::makeTxDbFromGFF(file = annotation))
    eByg <- GenomicFeatures::exonsBy(txdb, by = "gene")
    gr_shortlist <- eByg[names(eByg) %in% geneList]
    gr <- BiocGenerics::unlist(range(gr_shortlist)) # filter by names first before unlist
  }
  for (i in seq_along(systemPipeR::outpaths)) {
    print(names(systemPipeR::targets(args))[i])
    ## Get sequence information
    bam <- outpaths[i]
    si <- GenomicRanges::seqinfo(Rsamtools::BamFile(bam))
    # find the GRanges based on `geneList` OR `gr` OR `MappingRegion`
    if (!is.null(geneList)) {
      gr_target <- GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(gr), ranges = GenomicRanges::ranges(gr))
      if (is.null(gr_target)) stop("Please provide valid gene names.")
    } else if (!is.null(gr)) {
      gr_target <- GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(gr), ranges = GenomicRanges::ranges(gr))
    } else if (!is.null(MappingRegion)) {
      gr_target <- GenomicRanges::GRanges(GenomicRanges::seqnames(si), IRanges::IRanges(MappingRegion))
      gr_target
    } else {
      stop("One of the arguments 'geneList', 'gr' or 'chrHead' must not be empty.")
    }
    ## TODO: double check here. Need more testing
    gr_outer <- range(gr_target) # obtain the outer boundary for the following subsetting
    ## Reference genome and annotatio
    ## TODO: need more testing
    if(truncate_refs){
      ## subsetting
      .subsetRef(reference, truncate_refs, gr, gr_outer, ref_outname)
      ## subsetting gff annotation file
      .subsetAnno(anno, truncate_refs, gr, gr_outer, outdir)
    }
    ## faster option from 'readGAlignments'
    ## get GR ids
    param <- Rsamtools::ScanBamParam(which = gr_target, what = c("qname", "pos"))
    aligns <- unname(Rsamtools::scanBam(bam, param = param))
    # aligns <- unname(aligns)
    elts <- stats::setNames(Rsamtools::bamWhat(param), Rsamtools::bamWhat(param))
    aligns <- lapply(elts, function(elt) BiocGenerics::unlist(lapply(aligns, "[[", elt)))
    names(aligns) <- elts
    keepids <- aligns$qname[!is.na(aligns$pos)]
    keepids
    # whether the reads id will be sampled
    if (!is.null(sample_range)) {
      myN <- sample(sample_range, 1)
      keepids <- sample(IRanges::unique(keepids), myN, replace = TRUE)
    } else {
      keepids <- keepids
    }
    ## Filter fastq and writeout subset
    filter <- function(x) x[as.character(ShortRead::id(x)) %in% keepids]
    if (!is.null(systemPipeR::infile1(args)[i])) {
      if (id_read_number) keepids <- paste(keepids, "/1", sep = "")
      ShortRead::filterFastq(files = systemPipeR::infile1(args)[i], destinations = file.path(outdir, basename(systemPipeR::infile1(args)[i])), filter = filter)
      if (!silent) {
        cat(
          "\t", "Data subset written to the file:",
          file.path(outdir, BiocGenerics::basename(systemPipeR::infile1(args)[i])), "\n"
        )
      }
      gc()
    }
    if (!is.null(systemPipeR::infile2(args)[i])) {
      if (id_read_number) keepids <- IRanges::gsub("/1", "/2", keepids)
      ShortRead::filterFastq(files = systemPipeR::infile2(args)[i], destinations = file.path(outdir, BiocGenerics::basename(systemPipeR::infile2(args)[i])), filter = filter)
      if (!silent) {
        cat(
          "\t", "Data subset written to the file:",
          file.path(outdir, basename(systemPipeR::infile2(args)[i])), "\n"
        )
      }
      gc()
    }
  }
}
## Usage:
# getSubsetReads(args, MappingRegion = 1:900, sample_range = 800:900, outdir = "data/subset/", silent = FALSE)
# getSubsetReads(args, MappingRegion = 1:900, sample_range = NULL, outdir = "data/subset/", silent = FALSE)

#############################
## .subsetRef function ##
#############################

.subsetRef <- function(reference, truncate_refs, gr, gr_outer, ref_outname) {
  refpath <- normalizePath(reference)
  if (truncate_refs){
    sref <- getSeq(FaFile(refpath), gr)
    srefl <- split(sref, names(sref))
    subsetRef <- DNAStringSet(lapply(srefl, unlist))
    writeXStringSet(subsetRef, filepath = ref_outname)
  } else {
    sref <- getSeq(FaFile(refpath), gr_outer)
    writeXStringSet(sref, filepath = ref_outname)
  }
}

#############################
## .subsetAnno function ##
#############################
.subsetAnno <- function(annotation, truncate_refs, gr, gr_outer, annot_outname) {
  anno <- rtracklayer::import(annotation)
  if (truncate_refs){
    annosub <- anno[anno %within% gr]
    rdc <- reduce(annosub, with.revmap = TRUE, ignore.strand = TRUE)
    rmp <- mcols(rdc)$revmap
    grl <- relist(annosub[unlist(rmp)], rmp)
    widthList <- split(width(rdc), as.character(seqnames(rdc))) # base:split since width(rdc) is a vector
    new_end <- unlist(unname(lapply(widthList, cumsum)))
    subsetGFF <- unlist(shift(grl, new_end - end(rdc)))           
  } else {
    annosub <- anno[anno %within% gr_outer]
    grl <- split(annosub, seqnames(annosub), drop = TRUE)
    subsetGFF <- unlist(shift(grl, 1 - start(gr_outer)))
  }
  con <- file(annot_outname, open = "w")
  rtracklayer::export(subsetGFF, con, format = "gff3")
  close(con)
}
