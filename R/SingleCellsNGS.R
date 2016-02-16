#' @name SingleCellsNGS
#' @title SingleCellsNGS
#' @docType package
#' @description  An S4 class to visualize Expression data from a single cell NGS experiment.
#' @slot data a data.frame containing the expression values for each gene x sample (gene = row)
#' @slot samples a data.frame describing the columnanmes in the data column
#' @slot annotation a data.frame describing the rownames in the data rows
#' @slot outpath the default outpath for the plots and tables from this package
#' @slot name the name for this package (all filesnames contain that)
#' @slot rownamescol the column name in the annotation table that represents the rownames of the data table
#' @slot sampleNamesCol the column name in the samples table that represents the colnames of the data table
#' @slot stats the stats list for all stats created in the object
setClass( 
		Class='SingleCellsNGS', 
		representation = representation ( 
		##	NGS = 'binary'
		),
		contains='NGSexpressionSet',
		prototype(outpath ='', name = 'SingleCellsNGS', 
				rownamescol=NA_character_, 
				sampleNamesCol=NA_character_, 
				stats=list() )
)

extends("SingleCellsNGS","NGSexpressionSet")

#' @name SingleCellsNGS
#' @aliases SingleCellsNGS,SingleCellsNGS-method
#' @rdname SingleCellsNGS-methods
#' @docType methods
#' @description  create a new SingleCellsNGS object This object is mainly used for plotting the
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
#' @return A new SingleCellsNGS object
#' @title description of function SingleCellsNGS
#' @export 
setGeneric('SingleCellsNGS', ## Name
		function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = NULL) { ## Argumente der generischen Funktion
			standardGeneric('SingleCellsNGS') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('SingleCellsNGS', signature = c ('data.frame'),
		definition = function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = NULL) {
			x <- StefansExpressionSet( dat, Samples,  Analysis = Analysis, name=name, namecol=namecol, namerow= namerow, usecol= usecol, outpath =  outpath)
			as(x,'SingleCellsNGS')
		} )

#' @name normalize
#' @aliases normalize,NGSexpressionSet-method
#' @rdname normalize-methods
#' @docType methods
#' @description  
#' normalize the expression data by subsampling as described in PMID 24531970
#' @param x The SingleCellsNGS object
#' @param reads the required read depth
#' @param name the name of the new object
#' @return the normalized data set (original data stored in slot 'raw'
#' @title description of function normalize
#' @export 
setGeneric('normalize', ## Name
		function ( object , ..., reads= 600, name='normalized' ) { ## Argumente der generischen Funktion
			standardGeneric('normalize') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)
setMethod('normalize', signature = c ('SingleCellsNGS'),
		definition = function (  object, ..., reads=600, name='normalized' ) {
			if ( length( object@samples$counts ) == 0 ) {
				object@samples$counts <- apply( object@data, 2, sum)
			}
			object <- drop.samples( object, object@samples[which(object@samples$counts < reads), object@sampleNamesCol ] 
					, name=name )
			if ( ! object@snorm ){
				object@raw <- object@data
				object@snorm = TRUE
			}
			## resample the data
			n <- nrow(object@raw)
			object@data[] <- 0
			for ( i in 1:ncol(object@raw) ) {
				d <- sample(rep ( 1:n, object@raw[,i]) , reads)
				t <- table(d)
				object@data[ as.numeric(names(t)),i] <- as.numeric(t)
			}
			object
		}
)



#' @name z.score
#' @aliases z.score,SingleCellsNGS-method
#' @rdname z.score-methods
#' @docType methods
#' @description  Use most of'StefansExpressionSet'functionallity with minor changes to NGS data (like normalizing)
#' @description  This package is mainly meant to plot the data - all important quantification or gene
#' @description  annotation is performed using DEseq functionallity. Read a set of bam files and perform
#' @description  the quantification (better do that without using this function!)
#' @param x the SingleCellsNGS object
#' @title description of function z.score
#' @export 
setMethod('z.score', signature = c ('SingleCellsNGS'),
		definition = function ( m ) {
			
			if ( ! m@zscored ) {
				m@raw <- m@data
				ma  <- as.matrix(m@data)
				i = 0
				ret <- t(
						apply(ma,1, function (x) {
									i = i+1
									n <- which(x==0)
									if ( length(x) - length(n) > 1 ){
										x[-n] <- scale(as.vector(t(x[-n])))
									}
									else {
										x[] = -20
									}
									x[n] <- -20
									x}
						)
				)
				m@data <- data.frame(ret)
				m@zscored = TRUE
			}
			m
		})


#' @name gg.heatmap.list
#' @aliases gg.heatmap.list,SingleCellsNGS-method
#' @rdname gg.heatmap.list-methods
#' @docType methods
#' @description  uses ggplot2 to plot heatmaps
#' @param dat the StefansExpressionSet object
#' @param glist a list of probesets to plot (or all)
#' @param colrs a list of colors for the sample level boxes (or rainbow colors)
#' @param groupCol the column group in the samples table that contains the grouping strings
#' @param colCol the column group in the samples table that contains the color groups
#' @title description of function gg.heatmap.list
#' @export 
setMethod('gg.heatmap.list', signature = c ( 'SingleCellsNGS') ,
		definition = function (dat,glist=NULL, colrs=NULL, groupCol='GroupID', colCol=NULL) {
			
			if ( ! is.null(glist) ) {
				isect <- reduce.Obj ( dat, glist)
			}else {
				isect <- dat
			}
			if ( is.null(colCol)){
				colCol <- groupCol
			}
			if ( is.null(colrs) ){
				colrs = rainbow( length(unique(isect@samples[,colCol])))
			}
			if ( ! isect@gnorm ) {isect <- z.score(isect)}
			dat.ss <- melt.StefansExpressionSet ( isect, probeNames=isect@rownamescol, groupcol=groupCol,colCol=colCol)
			#dat.ss <- dat[which(is.na(match(dat$Gene.Symbol,isect))==F),]
			colnames(dat.ss) <- c( 'Gene.Symbol', 'Sample', 'Expression', 'Group', 'ColorGroup' )
			dat.ss$z <- ave(dat.ss$Expression, dat.ss$Gene.Symbol, FUN = function (x)  {
						n <- which(x==0)
						if ( length(x) - length(n) > 1 ){
							x[-n] <- scale(as.vector(t(x[-n])))
						}
						else {
							x[] = -20
						}
						x[n] <- -20
						x
					}
			)
			samp.cast <- reshape2::dcast(dat.ss,Gene.Symbol~Sample,mean,value.var="z")
			samp.mat <- as.matrix(samp.cast[,2:ncol(samp.cast)])
			ord.genes <-
					as.vector(samp.cast[hclust(dist(samp.mat),method="ward.D")$order,1])
			dat.ss$Gene.Symbol <- with(dat.ss,factor(Gene.Symbol,levels =
									unique(as.character(ord.genes))))
			dat.ss$Sample <- with(dat.ss,factor(Sample,levels =
									unique(as.character(Sample))))
			dat.ss$Group <- with(dat.ss,factor(Group,levels =
									unique(as.character(Group))))
			dat.ss$colrss <- colrs[dat.ss$Group]
			ss <-dat.ss[which(dat.ss$Gene.Symbol==dat.ss$Gene.Symbol[1]),]
			brks= c( -20.1, as.vector(quantile(dat.ss$z[which(dat.ss$z != -20)],seq(0,1,by=0.1)) ))
			brks = unique(brks)
			print ( brks )
			brks[length(brks)] = brks[length(brks)] + 0.1
			dat.ss$z <- cut( dat.ss$z, breaks= brks)
			
			
			p = ( ggplot2::ggplot(dat.ss, ggplot2::aes(x=SampleName,y=ProbeName))
						+ ggplot2::geom_tile(ggplot2::aes(fill=z)) 
						+ ggplot2::scale_fill_manual( values =  c( 'gray', gplots::bluered(length(brks) -2  )) ) 
						+ ggplot2::theme(
								legend.position= 'bottom',
								axis.text.x=ggplot2::element_blank(),
#axis.ticks.x=element_line(color=ss$colrss),
								axis.ticks.length=ggplot2::unit(0.00,"cm")
						)+ labs( y='') )
			if ( ncol(dat.ss) == 5 ){
				p <- p + ggplot2::facet_grid( ColorGroup ~ Group,scales="free", space='free')
			}else if ( ncol(dat.ss) == 4 ) {
				p <- p + ggplot2::facet_grid( . ~ Group,scales="free", space='free')
			}
			list ( plot = p, not.in = setdiff( glist, rownames(isect@data)) )
			
#			list ( plot = ggplot2::ggplot(dat.ss, ggplot2::aes(x=Sample,y=Gene.Symbol))
#							+ ggplot2::geom_tile(ggplot2::aes(fill=z))
#							+ ggplot2::scale_fill_manual( values = c( 'gray', gplots::bluered(length(brks) -2  )) ) 
#							+ ggplot2::theme(
#									legend.position= 'bottom',
#									axis.text.x=ggplot2::element_blank(),
#									axis.ticks.x=ggplot2::element_line(color=ss$colrss),
#									axis.ticks.length=ggplot2::unit(0.6,"cm")
#							)
#							+ ggplot2::labs( y=''),
#					not.in = setdiff( glist, rownames(isect@data)) )
		})
