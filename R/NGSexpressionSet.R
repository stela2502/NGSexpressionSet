library('StefansExpressionSet')
require('StefansExpressionSet')
library('DESeq')
require('DESeq')
library('Rsubread')
require('Rsubread')

#' @name NGSexpressionSet
#' @title NGSexpressionSet
#' @docType package
#' @description  An S4 class to visualize Expression data from a NGS experiment.
#' @slot data a data.frame containing the expression values for each gene x sample (gene = row)
#' @slot samples a data.frame describing the columnanmes in the data column
#' @slot annotation a data.frame describing the rownames in the data rows
#' @slot outpath the default outpath for the plots and tables from this package
#' @slot name the name for this package (all filesnames contain that)
#' @slot rownamescol the column name in the annotation table that represents the rownames of the data table
#' @slot sampleNamesCol the column name in the samples table that represents the colnames of the data table
#' @slot stats the stats list for all stats created in the object
#' @importClassesFrom StefansExpressionSet StefansExpressionSet
#' @exportClass NGSexpressionSet
setClass( 
		Class='NGSexpressionSet', 
		representation = representation ( 
			##	NGS = 'binary'
		),
		contains='StefansExpressionSet',
		prototype(outpath ='', name = 'NGSexpressionSet', 
				rownamescol=NA_character_, 
				sampleNamesCol=NA_character_, 
				stats=list() )
)


#' @name read.bams
#' @aliases read.bams,NGSexpressionSet-method
#' @rdname read.bams-methods
#' @docType methods
#' @description  Use most of'StefansExpressionSet'functionallity with minor changes to NGS data (like normalizing)
#' @description  This package is mainly meant to plot the data - all important quantification or gene
#' @description  annotation is performed using DEseq functionallity. Read a set of bam files and perform
#' @description  the quantification (better do that without using this function!)
#' @param bamFiles a file containing one file name per line
#' @param annotation the gene level annotation to use
#' @title description of function read.bams
#' @export 
setGeneric('read.bams', ## Name
	function ( bamFiles, annotation, GTF.featureType='exon', GTF.attrType = "gene_id", isPairedEnd = FALSE, nthreads = 2) { ## Argumente der generischen Funktion
		standardGeneric('read.bams') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('read.bams', signature = c ('NGSexpressionSet') ,
	definition = function ( bamFiles, annotation, GTF.featureType='exon', GTF.attrType = "gene_id", isPairedEnd = FALSE, nthreads = 2) {
	if (file.exists(bamFiles)){
		bamFiles <- readLines(bamFiles)
	}
	counts <- featureCounts(files =bamFiles,annot.ext = annotation ,isGTFAnnotationFile = TRUE,GTF.featureType = GTF.featureType,
		GTF.attrType = GTF.attrType,allowMultiOverlap=T, isPairedEnd =isPairedEnd , nthreads = nthreads)
	counts.tab <- cbind(counts$annotation,counts$counts)  # combine the annotation and counts
	counts.tab
})

extends("NGSexpressionSet", 'StefansExpressionSet')

#' @name NGSexpressionSet
#' @aliases NGSexpressionSet,NGSexpressionSet-method
#' @rdname NGSexpressionSet-methods
#' @docType methods
#' @description  create a new NGSexpressionSet object This object is mainly used for plotting the
#' @description  data
#' @param dat data frame or matrix containing all expression data
#' @param Samples A sample description table
#' @param Analysis If the samples table contains a Analysis column you can subset the data already here to the Analysis here (String). This table has to contain a column filenames that is expected to connect the sample table to the dat column names.
#' @param name The name of this object is going to be used in the output of all plots and tables - make it specific
#' @param namecol The samples table column that contains the (to be set) names of the samples in the data matrix
#' @param namerow This is the name of the gene level annotation column in the dat file to use as ID 
#' @param outpath Where to store the output from the analysis
#' @param annotation The annotation table from e.g. affymetrix csv data
#' @param newOrder The samples column name for the new order (default 'Order')
#' @return A new NGSexpressionSet object
#' @title description of function NGSexpressionSet
#' @export 
setGeneric('NGSexpressionSet', ## Name
	function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = NULL) { ## Argumente der generischen Funktion
		standardGeneric('NGSexpressionSet') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('NGSexpressionSet', signature = c ('data.frame'),
	definition = function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = NULL) {
	x <- StefansExpressionSet( dat, Samples,  Analysis = Analysis, name=name, namecol=namecol, namerow= namerow, usecol= usecol, outpath =  outpath)
	as(x,'NGSexpressionSet')
} )




#' @name normalize
#' @aliases normalize,NGSexpressionSet-method
#' @rdname normalize-methods
#' @docType methods
#' @description  normalize the expression data (sample wise)
#' @param x The NGSexpressionSet
#' @param readCounts The number of reads from each bam file or another value you want to normalize the data to
#' @param to_gene_length FALSE whether or not to normalize the data to gene length
#' @param geneLengthCol the column in the annotation data.frame to (in addition) normalize the genes to (e.g. trancript length)
#' @return the normalized data set (original data stored in NGS$raw
#' @title description of function normalize
#' @export 
setGeneric('normalize', ## Name
	function ( object , ..., readCounts=NULL, to_gene_length=FALSE, geneLengthCol='transcriptLength' ) { ## Argumente der generischen Funktion
		standardGeneric('normalize') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('normalize', signature = c ('NGSexpressionSet'),
		definition = function (  object, readCounts=NULL, to_gene_length=FALSE, geneLengthCol='transcriptLength' ) {
			if ( ! object@snorm ){
				if ( is.null( readCounts ) ) {
					readCounts <- as.vector( DESeq::estimateSizeFactorsForMatrix ( object@data) )
				}
				browser()
				object@samples$SizeFactor <- readCounts
				object@raw <- object@data
				object@data =  data.frame(t(apply(object@data,1, function(a) { a / readCounts } ) ))
				colnames(object@data) = colnames(object@raw)
				rownames(object@data) = rownames(object@raw)
				if (to_gene_length){
					for ( i in 1:nrow(object@data)) {
						object@data[i,] <- object@data[i,]/ (object@annotation[i,geneLengthCol ] / 1000)
					}
				}
				
				object@snorm=TRUE
			}
			object
		})


#' @name removeBatch
#' @aliases removeBatch,NGSexpressionSet-method
#' @rdname removeBatch-methods
#' @docType methods
#' @description  merge the groups in the dataset into one single column The sample description is
#' @description  restricted to the first entry in each group. Most variables in the sample description
#' @description  table might be useless after this collape
#' @param dataObj the NGSexpressionSet
#' @param what merge by 'median','mean','sd' or 'sum'
#' @return the collapsed NGSexpressionSet
#' @title description of function removeBatch
#' @export 
setGeneric('removeBatch', ## Name
	function ( x, phenotype ) { ## Argumente der generischen Funktion
		standardGeneric('removeBatch') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('removeBatch', signature = c ('NGSexpressionSet'),
	definition = function ( x, phenotype ) {
	if ( x@batchRemoved==1) {
		return (x)
	}
	browser()
	exprs =  dataOBJ$cds
	null <- which ( exprs == 0)
	exprs[null]<- 1
	log <- log(exprs)
	svseq <-  ComBat(dat=filtered.log , mod=mod1, batch=x@samples[,phenotype], par.prior=TRUE, prior.plots=FALSE )
	svseq <- exp(log)
	svseq[null] = 0
	x@cds <- svseq
	x@batchRemoved = 1
	x
})

#' @name check.depth
#' @aliases check.depth,NGSexpressionSet-method
#' @rdname check.depth-methods
#' @docType methods
#' @description this function checks how many reads are taken from the total dataset using the top 5 percent of genes only
#' @param x the NGSexpressionSet
#' @param gene an optional list of genes you want to select for
#' @param cutoff a cutoff that defines how many reads may be taken from the top 5 percent genes before a sample gets rejected
#' @title description of function check
#' @return a vector of sample names that did not pass the check
#' @export 
setGeneric('check.depth', ## Name
	function (x, genes=NULL, cutoff = 0.77 ) { ## Argumente der generischen Funktion
		standardGeneric('check.depth') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('check.depth', signature = c ('NGSexpressionSet'),
	definition = function (x, genes=NULL, cutoff = 0.77 ) {
	percent5 <-  reads.taken(x, 0.05, genes)
	names(which(percent5$reads.taken > cutoff )) ## obtained experimentally using Jennies HSC dataset
})

#' @name reads.taken
#' @aliases reads.taken,NGSexpressionSet-method
#' @rdname reads.taken-methods
#' @docType methods
#' @description  this check is testing how many percent of the total reads are taken from the first
#' @description  x percent of the genes this function calls reads.taken.NGSexpressionSet internally.
#' @param x the NGSexpressionSet
#' @param genes see reads.taken 
#' @param cutoff = 0.77 a good expression dataset should not use more than 77\% of the reads in the top 5\%
#' @param tf if you supply a list of genes here the same values as for all genes are calulated for this list
#' @return a list of samples which have passed the test
#' @title description of function reads.taken
#' @return a complicated list with al measured values
#' @export 
setGeneric('reads.taken', ## Name
	function ( x, percentile= 0.05, tf=NULL ) { ## Argumente der generischen Funktion
		standardGeneric('reads.taken') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('reads.taken', signature = c ('NGSexpressionSet'),
	definition = function ( x, percentile= 0.05, tf=NULL ) {
	top.genes <- list()
	reads.taken <- vector( 'numeric', ncol(x@data))
	nTF <- vector('numeric',  ncol(x@data)+1)
	percentile= 1- percentile
	for ( i in 1:ncol(x@data) ){
		qua <- quantile(x@data[,i], percentile)
		reads.taken [i] <- sum (x@data[which(x@data[,i] > qua),i] ) / sum(x@data[,i])
		top.genes[[i]] <- rownames(x@data) [ which(x@data[,i] > qua) ]
		if ( ! is.null(tf) ){
		nTF[i] <- length( intersect( tf, str_replace(top.genes[[i]],  "\\.[0-9]*" ,'') ) )
		}
	}
	names( reads.taken) <- colnames(x@data)
	names(top.genes) <- colnames(x@data)
		
	inter <- intersect( top.genes[[1]], top.genes[[2]])
	for (i in 3:ncol(x@data) ) { inter <- intersect( inter, top.genes[[i]]) }
	if ( ! is.null(tf) ) {nTF[ncol(x@data)+1] <- length(intersect(str_replace(inter,  "\\.[0-9]*" ,''), tf ) )}
	reads.taken.intersect <- vector( 'numeric', ncol(x@data))
	for ( i in 1:ncol(x@data) ){
		reads.taken.intersect [i] <- sum ( x@data[inter ,i] ) / sum(x@data[,i])
	}
	names( reads.taken.intersect) <- colnames(x@data)
	list( reads.taken = reads.taken, top.genes = top.genes, intersect = inter,reads.taken.intersect = reads.taken.intersect, nTF = nTF )
})
	

#' @name preprocess
#' @aliases preprocess,NGSexpressionSet-method
#' @rdname preprocess-methods
#' @docType methods
#' @description  calculates how many reads are consumed by the top genes performes a set of preprocess
#' @description  function used by DEseq Do not use
#' @param x the NGSexpressionSet
#' @param percentile how many of the top genes to use (0.05)
#' @param tf further subselect the top genes using these genes
#' @param x the NGSexpressionSet
#' @return a list of values - please read the source code of this function
#' @title description of function preprocess
#' @export 
setGeneric('preprocess', ## Name
	function (x) { ## Argumente der generischen Funktion
		standardGeneric('preprocess') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('preprocess', signature = c ('NGSexpressionSet'),
	definition = function (x) {
	if ( ! exists(where=x, 'vsdFull')){
		condition <- as.factor(x@samples$GroupName)
		print ( condition )
		x@cds <- newCountDataSet(x@data, condition)
		x@cds <- estimateSizeFactors(x@cds) 
		sizeFactors(x@cds)
		x@cds <- estimateDispersions(x@cds) 
		x@vsdFull = varianceStabilizingTransformation( x@cds )
	}
	x
})


#' @name anayse_working_set
#' @aliases anayse_working_set,NGSexpressionSet-method
#' @rdname anayse_working_set-methods
#' @docType methods
#' @description  zscore the ngs dataset
#' @param m the NGSexpressionSet
#' @return The z scored NGSexpressionSet
#' @title description of function anayse_working_set
#' @export 
setGeneric('anayse_working_set', ## Name
	function ( dataOBJ, name,  p_values = c ( 0.1,  1e-2 ,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8 ),geneNameCol= "mgi_symbol", batchCorrect=NULL) { ## Argumente der generischen Funktion
		standardGeneric('anayse_working_set') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('anayse_working_set', signature = c ('NGSexpressionSet'),
	definition = function ( dataOBJ, name,  p_values = c ( 0.1,  1e-2 ,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8 ),geneNameCol= "mgi_symbol", batchCorrect=NULL) {
	dataOBJ <- preprocess.NGSexpressionSet ( dataOBJ )
	if ( ! is.null(batchCorrect) ){
		removeBatch.NGSexpressionSet( dataOBJ ,batchCorrect )
	}
	png(file=paste(dataOBJ@name,"_",version,'_interesting_samples_clean_up.png',sep=''), width=800, height=800) 
	plot(hclust(as.dist(1-cor(dataOBJ@data))),main=paste(dataOBJ@name,version, ' samples')) 
	dev.off()
	dataOBJ <- do_comparisons.NGSexpressionSet ( dataOBJ, geneNameCol=geneNameCol)
	#dataOBJ@expr <- exprs(dataOBJ@vsdFull)
	dataOBJ <- get_gene_list.NGSexpressionSet(dataOBJ,p_values, geneNameCol=geneNameCol)	
	dataOBJ
})


#' @name createStats
#' @aliases createStats,NGSexpressionSet-method
#' @rdname createStats-methods
#' @docType methods
#' @description  calculate staistics on all possible groupings using the DEseq nbinomTest test Both
#' @description  together create the group 'A vs. B'
#' @param x the NGSexpressionSet
#' @param condition the grouping column in the samples data.frame
#' @param files write the statistics tables (FALSE)
#' @param A the first component to analyze (Group A)
#' @param B the second component to analyze (Group B)
#' @return the NGSexpressionSet with a set of ststs tables
#' @title description of function createStats
#' @export 
setGeneric('createStats', ## Name
	function (x, condition, files=F, A=NULL, B=NULL) { ## Argumente der generischen Funktion
		standardGeneric('createStats') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('createStats', signature = c ('NGSexpressionSet'),
	definition = function (x, condition, files=F, A=NULL, B=NULL) {
	if ( nrow(x@data) < 2e+4 ) {
	    stop ( "Please calculate the statistics only for the whole dataset!" )
	}
	if ( length( grep ( condition, colnames(x@samples))) > 0 ) {
		condition = factor( x@samples[,condition] )
	}
	if ( exists( 'raw', where=x) ){
		cds <- newCountDataSet(x@raw, condition)
	}else {
		cds <- newCountDataSet(x@data, condition)
	}
	cds <- estimateSizeFactors(cds)
	# sizeFactors(cds)
	cds <- estimateDispersions(cds)
	vsdFull = varianceStabilizingTransformation( cds )
	res <- list()
	ln= 0
	na <- c()
	conditions <- as.vector(unique(condition))
	if ( ! is.null(A) && ! is.null(B)) {
		ln = ln +1
		na[ln] = paste( A, B ,sep=' vs. ')
		res[[ln]] <- nbinomTest(cds, A, B )
	}
	else {
		for ( i in 1:(length(conditions)-1) ){
			for ( a in (i+1):length(conditions) ){
				ln = ln +1
				print ( paste( conditions[i], conditions[a], sep='_vs_') )
				na[ln] = paste( conditions[i], conditions[a],sep=' vs. ')
				res[[ln]] <- nbinomTest(cds, conditions[i], conditions[a])
				#res[[ln]] <- res[[res[[1]]+1]][which(res[[res[[1]]+1]]$padj < 0.2),]
			}
		}
	}
	names(res) <- na
	if ( exists( 'stats', where=x )) {
		x@stats <- c( x@stats, res)
	}else {
		x@stats <- res
	}
	if ( files ) {
		writeStatTables( x )
	}
	x
})
