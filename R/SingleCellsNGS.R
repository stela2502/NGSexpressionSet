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
		function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = '') { ## Argumente der generischen Funktion
			standardGeneric('SingleCellsNGS') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('SingleCellsNGS', signature = c ('data.frame'),
		definition = function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = '') {
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
										if (length(n) == 0 ){
											x <-  scale(as.vector(t(x)))
										}
										else {
											x[-n] <- scale(as.vector(t(x[-n])))
											x[n] <- -20
										}
										
									}
									else {
										x[] = -20
									}
									x}
						)
				)
				m@data <- data.frame(ret)
				colnames(m@data)<- colnames(m@raw)
				m@zscored = TRUE
			}
			m
		})

#' @name defineHeatmapColors
#' @aliases defineHeatmapColors,SingleCellsNGS-method
#' @rdname defineHeatmapColors-methods
#' @docType methods
#' @description  uses ggplot2 to plot heatmaps
#' @param merged the merged data object with the Expression column that should be colored
#' @param colrs and optional colors vector( gray + bluered for data and rainbow for samples)
#' @title description of function gg.heatmap.list
#' @return a list with the modified merged table and the colors vector
#' @export 
setGeneric('defineHeatmapColors', ## Name
		function ( melted,..., colrs=NULL) { ## Argumente der generischen Funktion
			standardGeneric('defineHeatmapColors') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('defineHeatmapColors', signature = c ( 'data.frame') ,
		definition = function (melted, colrs=NULL ){
			if ( is.factor( melted$Expression )) {
				## here might be some row grouping going on!
				d <- levels(melted$Expression)[melted$Expression]
				prob.id <- which(is.na(as.numeric(d))==T)
				treat.separate <- unique(d[prob.id])
				n <- as.numeric(d[-prob.id])
				brks= c( -20.1, as.vector(quantile(n[which(n != -20)],seq(0,1,by=0.1)) ))
				brks = unique(brks)
				d[-prob.id]  <- brks [cut( n, breaks= brks)]
				melted$Expression <- factor( d, levels= c(brks, treat.separate ) )
				colors= c(
						'gray', 
						gplots::bluered(length(brks) -2  ), ## the expression
						rainbow( length(treat.separate) ) ## the sample descriptions
				)
			}
			else {
				n <- as.numeric(melted$Expression )
				brks= c( -20.1, as.vector(quantile(n[which(n != -20)],seq(0,1,by=0.1)) ))
				brks = unique(brks)
				melted$Expression <- factor( brks [cut( n, breaks= brks)] , levels= c(brks) )
				
				colors= c(
								'gray', 
								gplots::bluered(length(brks) -2  ) ## the expression
				)
			}
			list (melted = melted, colors = colors)
		}
)


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
			
			print ( "You are using the latest version!")
			if ( ! is.null(glist) ) {
				isect <- reduce.Obj ( dat, glist)
			}else {
				isect <- dat
			}
			if ( is.null(colrs) ){
				colrs = rainbow( length(unique(isect@samples[,colCol])))
			}
			if ( ! isect@zscored ) {isect <- z.score(isect)}
			dat.ss <- melt.StefansExpressionSet ( isect, probeNames=isect@rownamescol, groupcol=groupCol,colCol=colCol)
			#dat.ss <- dat[which(is.na(match(dat$Gene.Symbol,isect))==F),]
			colnames(dat.ss) <- c( 'Gene.Symbol', 'Sample', 'Expression', 'Group', 
					paste('ColorGroup', 1:10) )[1:ncol(dat.ss)]
			r <- defineHeatmapColors( dat.ss )
			dat.ss <- r$melted
			ord.genes <- rownames(isect@data)[hclust(dist(isect@data),method="ward.D2")$order]
			if ( ! is.null(colCol) ){
				ord.genes <- c( ord.genes,colCol  )
			}
			dat.ss$Gene.Symbol <- with(dat.ss,factor(Gene.Symbol,levels =
									unique(as.character(ord.genes))))
			dat.ss$Sample <- with(dat.ss,factor(Sample,levels =
									unique(as.character(Sample))))
			dat.ss$Group <- with(dat.ss,factor(Group,levels =
									unique(as.character(Group))))
			dat.ss$colrss <- colrs[dat.ss$Group]
			ss <-dat.ss[which(dat.ss$Gene.Symbol==dat.ss$Gene.Symbol[1]),]
			p = ( ggplot2::ggplot(dat.ss, ggplot2::aes(x=Sample,y=Gene.Symbol))
						+ ggplot2::geom_tile(ggplot2::aes(fill=Expression)) 
						+ ggplot2::scale_fill_manual( values =  r$colors ) 
						+ ggplot2::theme(
								legend.position= 'bottom',
								axis.text.x=ggplot2::element_blank(),
								axis.ticks.length=ggplot2::unit(0.00,"cm")
						)+ ggplot2::labs( y='') )
			if ( ncol(dat.ss) == 6 ){
				p <- p + ggplot2::facet_grid( colrss ~ Group,scales="free", space='free')
			}else if ( ncol(dat.ss) == 5 ) {
				p <- p + ggplot2::facet_grid( . ~ Group,scales="free", space='free')
			}
			list ( plot = p, not.in = setdiff( glist, rownames(isect@data)) )
		})
