
setMethod("getGenePositions", c("ensemblGenome","character"), function(object, by, force = FALSE, ...)

{

    # - - - - - - - - - - - - - - - - - - - - - - - - - - #

    if(!is.logical(force))

        stop("[getGenePositions.ensemblGenome] force must be logical!")



    # Copy of table will be in ev -> positions

    # need only once be calculated.

    if(exists("genes", where=object@ev, inherits=FALSE) & !force)

        return(object@ev$genes)



    # - - - - - - - - - - - - - - - - - - - - - - - - - - #

    # Differing by 'gene_id' makes sense for Ensembl

    # - - - - - - - - - - - - - - - - - - - - - - - - - - #

    if(missing(by))

        by <- "gene_id"

    if(!is.character(by))

        stop("'by' must be character!")





    if(!exists("gtf", where=object@ev, inherits=FALSE))

        return(NULL)

    if(!is.data.frame(object@ev$gtf))

        stop("gtf-table must be data.frame!")



    if(is.na(match("gene_name", names(object@ev$gtf))))

        stop("No 'gene_name' data found!")



    # - - - - - - - - - - - - - - - - - - - - - - - - - - #

    # Different gene-identifications

    # - - - - - - - - - - - - - - - - - - - - - - - - - - #

    if(by=="gene_id")

    {

        # Get (sorted) gene_id's

        # (as.numeric(genes)==1:n! asc sorted!)

        if(is.factor(object@ev$gtf$gene_id))

          genes <- factor(levels(object@ev$gtf$gene_id))

        else

          genes <- sort(unique(object@ev$gtf$gene_id))





        n <- length(genes)

        # Point back into source table

        mtc <- match(genes, object@ev$gtf$gene_id)



        # Min start position (table has same order as genes!)

        mig <- summaryBy(start~gene_id, data=object@ev$gtf, FUN=min)

        # Max end   position (table has same order as genes!)

        mxg <- summaryBy(end~gene_id, data=object@ev$gtf, FUN=max)

    }

    else if(by=="gene_name")

    {

        # Get (sorted) gene_id's

        # (as.numeric(genes)==1:n! asc sorted!)

        genes <- factor(levels(object@ev$gtf$gene_name))

        n <- length(genes)

        # Point back into source table

        mtc <- match(genes, object@ev$gtf$gene_name)



        # Min start position (table has same order as genes!)

        mig <- summaryBy(start~gene_name, data=object@ev$gtf, FUN=min)

        # Max end   position (table has same order as genes!)

        mxg <- summaryBy(end~gene_name, data=object@ev$gtf, FUN=max)

    }

    else

        stop("[getGenePositions.ensemblGenome] by must be 'gene_id' or 'gene_name'!")



    # - - - - - - - - - - - - - - - - - - - - - - - - - - #

    # Assemble result

    # - - - - - - - - - - - - - - - - - - - - - - - - - - #

    if(is.na(match("gene_biotype", names(object@ev$gtf))))

    {

        res <- data.frame(id = 1:n,  gene_id = object@ev$gtf$gene_id[mtc],

                          gene_name = object@ev$gtf$gene_name[mtc],

                          seqid = object@ev$gtf$seqid[mtc],

                          start = mig[, 2],

                          end = mxg[, 2],

                          strand = object@ev$gtf$strand[mtc])



    }else{

        res <- data.frame(id = 1:n, gene_id = object@ev$gtf$gene_id[mtc],

                          gene_name = object@ev$gtf$gene_name[mtc],

                          seqid = object@ev$gtf$seqid[mtc],

                          start = mig[, 2],

                          end=mxg[, 2],

                          strand = object@ev$gtf$strand[mtc],

                          gene_biotype = object@ev$gtf$gene_biotype[mtc])

    }



    message("[getGenePositions.ensemblGenome] Adding 'start_codon' and 'stop_codon' positions.")

    strt <- extractFeature(object, "start_codon")@ev$gtf

    mtc <- match(res$gene_id, strt$gene_id)

    stap <- strt$start[mtc]

    stam <- strt$end[mtc]

    res$start_codon <- ifelse(res$strand=='+', stap, stam)

    stpp <- extractFeature(object, "stop_codon")@ev$gtf

    mtc <- match(res$gene_id, stpp$gene_id)

    sttp <- stpp$start[mtc]

    sttm <- stpp$end[mtc]

    res$stop_codon <- ifelse(res$strand=='+', sttp, sttm)



    res <- res[order(res$seqid, res$start), ]

    assign("genes", res, envir=object@ev)

    return(invisible(res))

})

