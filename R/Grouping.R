

#' @name bestGrouping
#' @aliases 'bestGrouping,SingleCellsNGS-method
#' @docType methods
#' @description The function is using a randomForest classifier with 2000 trees to classify the given data using the given grooping
#' @description All groups that fail to be prediceted using the random forest are deemed ungrouped.
#' @description All groups where less than 50 percent of the total samples geting classified as being from that group fail.
#' @param x the single cells ngs object
#' @param group a vector of sample columns that should be checked (the most complex is used only)
#' @param bestColname the column name to store the best grouping in
#' @param cutoff the cutoff percentage where all groups showing less than this percentacge of remapped samples are dropped
#' @title description of function randomForest
#' @return a distRF object to be analyzed by pamNew
#' @export 
setGeneric('bestGrouping',
		function ( x, group , bestColname='QualifiedGrouping', cutoff=0.5){
			standardGeneric('bestGrouping')
		}
)
setMethod('bestGrouping', signature = c ('SingleCellsNGS'),
		definition = function (x, group, bestColname='QualifiedGrouping' , cutoff=0.5) {
			uObj <- paste( 'predictive RFobj', group )
			rf <- NULL
			if (  is.null( x@usedObj[[uObj]])){
				x@usedObj[[uObj]] <- randomForest( x= t(as.matrix(x@data)), y=factor(x@samples[, group]),ntree=2000 )
			}
			t <- table( observed= x@samples[,group ], predicted = x@usedObj[[uObj]]$predicted )
			i <- 0
			r <- vector('numeric', ncol(t))
			names(r) <- colnames(t)
			for (i in 1:nrow(t)) {
				if ( which(t[i,] == max(t[i,])) == i) {
					r[i]= max(t[i,]) / sum(t[i,])
				}
				else {
					r[i]= 0
				}
			}
			BAD <- which(r < cutoff )
			## remove an optional 'gr. ' from the group ids
			x@samples[,bestColname] <- as.numeric(str_replace_all( x@samples[, group], 'gr. ', ''))
			for ( b in BAD ) {
				x@samples[ which(x@samples[,bestColname] == b), bestColname] <- 0
			}
			
			for (i in 0:(length(table(x@samples[,bestColname]))-1)){
				modify <- which(x@samples[,bestColname] >= i )
				if ( length(modify) == 0 ) { break}
				while( length(which(x@samples[modify,bestColname] == i)) == 0 ){
					x@samples[modify,bestColname] = x@samples[modify, bestColname] -1
				}
			}
			x@samples[,bestColname] <- paste( 'gr.', x@samples[,bestColname])
			x
		}
)

#' @name predict.rf
#' @aliases 'predict.rf,SingleCellsNGS-method
#' @docType methods
#' @description simple prediction of groups using a random forest trained during the bestGrouping process
#' @param x the single cells ngs object
#' @param rf the random forst model to use for the classification
#' @param bestColname the column name to store the results in
#' @title description of function predict.rf
#' @return a SingleCellsNGS object including the results and storing the RF object in the usedObj list (bestColname)
#' @export 
setGeneric('predict.rf',
		function ( x, rf,  bestColname='predicted group using random forest'){
			standardGeneric('predict.rf')
		}
)
setMethod('predict.rf', signature = c ('SingleCellsNGS'),
		definition = function (x, rf, bestColname='predicted group using random forest') {
			predicted2 <-predict( rf, t(as.matrix(workingSet@data)) )
			
		}
)


#' @name rfCluster
#' @aliases 'rfCluster,SingleCellsNGS-method
#' @title rfCluster
#' @name rfCluster-methods
#' @docType methods
#' @description This fucntion uses the RFclust.SGE to create fandomForest based unsupervised clusters on a subset of the data.
#' @description Default is on 200 cells using all (provided) genes with 500 forests and 500 trees per forest for 5 repetitions.
#' @description You are asked to give a k numer of expected clusters (better too many than too little), classifies the total 
#' @description data using the 5 different unsupervised runs and all cluster ids from these runs are merged into the final cluster id.
#' @description This <summaryCol> will be part of the return objects samples table, together with a <usefulCol> where
#' @description all clusters with less than 10 cells have been merged into the 'gr. 0'.
#' @description The final results will be reported as new columns in the samples table containing the 'name'
#' @param x the single cells ngs object
#' @param email your email to use together with the SGE option
#' @param SGE whether to use the sun grid engine to calculate the rf grouping
#' @param rep how many repetitions for the random forest grouping should be run (default = 5)
#' @param slice how many processes should be started for each random forest clustering (default = 30)
#' @param bestColname the column name to store the results in
#' @param k the numer of expected clusters (metter more than to view)
#' @param subset how many cells should be randomly selected for the unsupervised clustering (default = 200)
#' @param name if you want to run multiple RFclusterings on e.g. using different input genes you need to specify a name (default ='RFclust')
#' @param pics create a heatmap for each grouping that has been accessed (in the outpath folder; default = FALSE)
#' @param nforest the numer of forests to grow for each rep (defualt = 500)
#' @param ntree the numer of trees per forest (default = 500)
#' @return a SingleCellsNGS object including the results and storing the RF object in the usedObj list (bestColname)
#' @export 
setGeneric('rfCluster',
		function ( x, rep=5, SGE=F, email, k=16, slice=30, subset=200, pics=F ,nforest=500, ntree=500, name='RFclust'){
			standardGeneric('rfCluster')
		}
)
setMethod('rfCluster', signature = c ('SingleCellsNGS'),
		definition = function ( x, rep=5, SGE=F, email, k=16, slice=30, subset=200, pics=F ,nforest=500, ntree=1000, name='RFclust') {
			summaryCol=paste( 'All_groups', name,sep='_')
			usefulCol=paste ('Usefull_groups',name, sep='_')
			n= paste(x@name, name,sep='_')
			m <- max(k)
			OPATH <- paste( x@outpath,"/",str_replace( x@name, '\\s', '_'), sep='')
			opath = paste( OPATH,"/RFclust.mp",sep='' )
			
			if ( ! dir.exists(OPATH)){
				dir.create( OPATH )
			}
			processed = FALSE
			single_res_col <- paste('RFgrouping',name)
			for ( i in 1:rep) {
				tname = paste(n,i,sep='_')
				
				if ( is.null(x@usedObj[['rfExpressionSets']][[tname]]) ){
					## start the calculations!
					if ( dir.exists(opath)){
						if ( opath == '' ) {
							stop( "Are you mad? Not giving me an tmp path to delete?")
						}
						system( paste('rm -f ',opath,"/*",tname,'*', sep='') )
					}else {
						dir.create( opath )
					}
					total <- ncol(x@data)
					if ( total-subset <= 20 ) {
						stop( paste( 'You have only', total, 'samples in this dataset and request to draw random',subset, "samples, which leaves less than 20 cells to draw on random!") )
					}
					if ( is.null(x@usedObj[['rfExpressionSets']])){
						x@usedObj[['rfExpressionSets']] <- list()
						x@usedObj[['rfObj']][[ i ]] <- list()
					}
					
					if ( length( x@usedObj[['rfExpressionSets']] ) < i  ) {
						x@usedObj[['rfExpressionSets']][[ i ]] <- drop.samples( x, colnames(x@data)[sample(c(1:total),total-subset)], tname )
						x@usedObj[['rfObj']][[ i ]] <- RFclust.SGE ( dat=x@usedObj[['rfExpressionSets']][[ i ]]@data, SGE=SGE, slice=slice, email=email, tmp.path=opath, name= tname )
					}
					names(x@usedObj[['rfExpressionSets']]) [i] <- tname
					names(x@usedObj[['rfObj']]) [i] <- tname
					x@usedObj[['rfObj']][[ i ]] <- runRFclust ( x@usedObj[['rfObj']][[ i ]] , nforest=nforest, ntree=ntree, name=tname )
					if ( SGE){
						print ( "You should wait some time now to let the calculation finish! check: system('qstat -f') -> re-run the function")
					}
					else {
						print ( "You should wait some time now to let the calculation finish! -> re-run the function")
						print ( "check: system( 'ps -Af | grep Rcmd | grep -v grep')")
					}
				}
				else {
					
					## read in the results
					try ( x@usedObj[['rfObj']][[ i ]] <- runRFclust ( x@usedObj[['rfObj']][[ i]] , nforest=nforest, ntree=ntree, name=tname ) )
					if ( ! is.null(x@usedObj[['rfObj']][[ i ]]@RFfiles[[tname]]) ){
						stop( "please re-run this function later - the clustring process has not finished!")
					}
					for ( a in k ){
						x@usedObj[["rfExpressionSets"]][[i]]@samples <- 
								x@usedObj[["rfExpressionSets"]][[i]]@samples[ ,
										is.na(match ( colnames(x@usedObj[["rfExpressionSets"]][[i]]@samples), paste('group n=',a) ))==T 
								]
					}
					groups <- createGroups( x@usedObj[['rfObj']][[i]], k=k, name=tname )
					x@usedObj[['rfExpressionSets']][[i]]@samples <- cbind ( x@usedObj[['rfExpressionSets']][[i]]@samples, groups[,3:(2+length(k))] )
					
					le <- ncol(x@usedObj[['rfExpressionSets']][[i]]@samples)
					colnames(x@usedObj[['rfExpressionSets']][[i]]@samples)[(le-length(k)+1):le] <- paste('group n=',k)
					
					## create the required RF object
					m <- max(k)
					x@usedObj[['rfExpressionSets']][[i]] <- bestGrouping( x@usedObj[['rfExpressionSets']][[i]], group=paste('group n=', m), bestColname = paste('OptimalGrouping',m ,name) )
					## the 'predictive RFobj group n=' object is created by the bestGrouping call
					x@samples[, paste( single_res_col, i) ] <-
							predict( x@usedObj[['rfExpressionSets']][[i]]@usedObj[[paste( 'predictive RFobj group n=',m) ]], t(as.matrix(x@data)) )
					if ( pics ){
						fn <- paste(OPATH,'/heatmap_rfExpressionSets_',i,'.png', sep='')
						png ( file=fn, width=800, height=1600 )
						gg.heatmap.list( x, groupCol=paste( single_res_col , i) )
						dev.off()
						print ( paste('heatmap stored in', fn) )
					}
					print ( paste("Done with cluster",i))
					processed = TRUE
				}
			}
			if ( processed ) {
				try ( {combine <- identifyBestGrouping( x, c( paste(single_res_col, 1:rep)) )} , silent =T)
				if ( all.equal(as.vector(combine$res), rep('', rep)) ) {
					print( 'No really usful grouping of the data obtained - I recommend re-run with more trees/forests and a new name')
					x <- combine$x
				}
				else {
					x <- combine$x
					colnames(x@samples)[which( colnames(x@samples) == names(combine$res)[1] )] <- usefulCol
					if ( pics ){
						fn <- paste(OPATH,'/heatmap_',str_replace( usefulCol, '\\s', '_'),'.png', sep='')
						png ( file=fn, width=800, height=1600 )
						gg.heatmap.list( x, groupCol= usefulCol )
						dev.off()
						print ( paste('heatmap stored in', fn ))
					}
				}
				x@usedObj$combinationAnalysis <- list ( 'initial_significants' = combine$names, 'merged_significants' = combine$res )				
			}
			x		
		}
)

#' @name qualityTest
#' @aliases qualityTest,SingleCellsNGS-method
#' @rdname qualityTest-methods
#' @docType methods
#' @description This function calculates an anova test for the groupCol nd all data.
#' @return A list containing the updated 'x' SingleCellsNGS object and res a vector of pasted gene names that turned out to be significantly different.
#' @param x The SingleCellsNGS object
#' @param groups the group name to clister the data (default=NULL)
#' @param cut the p value cut off (BenjaminHochberg corrected; default=0.05)
#' @param numbers if true return amount of significant genes, not the names
#' @title description of function qualityTest
#' @export 
setGeneric('qualityTest', ## Name
		function (x, groups=NULL, cut=0.05, numbers=F  ) { ## Argumente der generischen Funktion
			standardGeneric('qualityTest') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('qualityTest', signature = c ('SingleCellsNGS'),
		definition = function (x, groups=NULL, cut=0.05, numbers=F ) {
			if (is.null(groups) ) {
				stop( "please give me a vector of columns to check!" )
			}
			res <- vector( 'numeric', length(groups))
			names(res) <- groups
			for ( i in 1:length(groups) ) {
				x <- simpleAnova( x,groups[i] )
				grname <- paste('simpleAnova', groups[i])
				g <-which(x@stats[[grname]][,3] < cut)
				
				if ( numbers ){
					res[i] <- length( x@stats[[grname]][which(x@stats[[grname]][,3]< cut),1] )
				}
				else {
					if ( length( g ) > 0 ) {
						res[i] <- paste(collapse=" ",as.vector(x@stats[[grname]][which(x@stats[[grname]][,3] < cut),1] ))
					}
					else{
						res[i] <- ''
					}
				}
			}
			list ( res = res, x=x)
		} )

#' @name identifyBestGrouping
#' @aliases identifyBestGrouping,SingleCellsNGS-method
#' @rdname identifyBestGrouping-methods
#' @docType methods
#' @description This function compares several groupings to each other. 
#' @description In addition the groupings are combined in a way, that if a sample end up in group 1 in grouing A but in group 2 in groupung B
#' @description it might be something else than a sample ending up in group 1 for both.
#' @description The quality of the groupings will be accessed in two stated (1) there should not be too many small groups and
#' @description (2) the genes should be differentially expressed in these groups. The differential expression is accessed 
#' @description using a straight forward anova approach (excluding the not called genes for a SingleCellsNGS object).
#' @param x the SingleCellsNGS
#' @param groups the colnames (samples) that contain the grouping information (default)
#' @param namePrefix a common name prefix for this analysis
#' @param cut define the p value cut off for the test ( if too view signfican genes are detected at 0.05)
#' @return A list containing the number of significant genes in each possible combination of groups and the SingleCellsNGS object containg all annotations and stats.
#' @title description of function identifyBestGrouping
#' @export 
setGeneric('identifyBestGrouping', ## Name
		function (x, groups, namePrefix='identifyBestGrouping', cut=0.05) { ## Argumente der generischen Funktion
			standardGeneric('identifyBestGrouping') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('identifyBestGrouping', signature = c ('SingleCellsNGS'),
		definition = function (x, groups, namePrefix='identifyBestGrouping', cut=0.05) {
			## I do not want to loose all so get me the most out of this gouping
			ret <- list( 'x' = 0, groups = 0 )
			groupLengthT <- function ( a, g ) {
				sum (which(t) < 10 ) < 20 && length(t) < 20
			}
			groupPaste <- function ( a, g, name ) {
				if ( ! is.na(match (name,colnames(a@samples)) ) ) {
					stop( paste( "The column",name,'already exists - STOP') )
				}
				a@samples[, name ] <- apply( a@samples[, g ],1, function (x) {paste(x,collapse= ' ') } )
				a
			}

			## first oder the group cols by the output of identifyBestGrouping			
			te <- qualityTest ( x, groups, numbers=T , cut=cut)
			x <- te$x
			groups <- groups[ order( te$res, decreasing =T ) ]
			names = c(paste( namePrefix, 'All (',length(groups),')' ))
			x <- groupPaste (x, groups, names)
			for ( i in 2:length(groups) ) {
				## then paste them best to worst together and check where you get better stats
				x <- groupPaste( x, groups[1:i], paste( namePrefix, i,'/',length(groups) ) )
				names<- c(names, paste( namePrefix, i,'/',length(groups) ))
			}
			ret <- qualityTest ( x, names, cut=cut )
			ret$names <- te$res
			ret
} )



#' @name update.measurements
#' @aliases update.measurements,SingleCellsNGS-method
#' @rdname update.measurements-methods
#' @docType methods
#' @description This method is inspired by PMID:24531970; it calculates the expression variance over all samples
#' @description and the fraction of cells expressing for each gene. This data is used by interesting_genes() to identify
#' @description putatively interesting genes in a SingleCellNGS object.
#' @param x the SingleCellsNGS object
#' @title description of function update.measurements
#' @export 
setGeneric('update.measurements', ## Name
		function ( x ) { ## Argumente der generischen Funktion
			standardGeneric('update.measurements') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('update.measurements', signature = c ('SingleCellsNGS'),
		definition = function ( x ) {
			x@annotation$mean_var <- apply( x@data , 1, function (x ) { 
							a<-x[which(x > 0)]
							if ( length(a) < 5 ){
								0
							}
							else {
								mean(a) / var(a) 
							}
						} )
			x@annotation$fraction_expressing <-apply( x@data , 1, function (x ) { length(which(x > 0)) / length(x) })
			x
		} 
)
#' @name interesting_genes
#' @aliases interesting_genes,SingleCellsNGS-method
#' @rdname interesting_genes-methods
#' @docType methods
#' @description This method uses the values created by update.measurements() and selects the most interesting subset of genes.
#' @param x the SingleCellsNGS object
#' @param Min the minimal fraction of cells expresing the gene default=0.4
#' @param Max the maximal fraction of cells expresing the gene default=0.8 
#' @param genes how many genes you want to get back default=100
#' @title description of function interesting_genes
#' @export 
setGeneric('interesting_genes', ## Name
		function (x, Min=0.4, Max=0.8 , genes=100) { ## Argumente der generischen Funktion
			standardGeneric('interesting_genes') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('interesting_genes', signature = c ('SingleCellsNGS'),
		definition = function (x, Min=0.4, Max=0.8 , genes=100) {
			if ( is.null(x@annotation$mean_var) ) {
				x <- update.measurements ( x )
			}
			possible <- which( x@annotation$fraction_expressing > Min & x@annotation$fraction_expressing < Max )
			order.possible <- order(x@annotation$mean_var[possible])
			rownames(x@data)[possible[order.possible][1:genes]]
		} 
)




