

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
			uObj <- paste( 'RFobj', group )
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
#' @docType methods
#' @description This fucntion uses the RFclust.SGE to create fandomForest based unsupervised clusters on a subset of the data.
#' @description Default is on 200 cells using all (provided) genes with 500 forests and 500 trees per forest for 5 repetitions.
#' @description You are asked to give a k numer of expected clusters (better too many than too little), classifies the total 
#' @description data using the 5 different unsupervised runs and all cluster ids from these runs are merged into the final cluster id.
#' @description This <summaryCol> will be part of the return objects samples table, together with a <usefulCol> where
#' @description all clusters with less than 10 cells have been merged into the 'gr. 0'.
#' @param x the single cells ngs object
#' @param email your email to use together with the SGE option
#' @param SGE whether to use the sun grid engine to calculate the rf grouping
#' @param rep how many repetitions for the random forest grouping should be run (default = 5)
#' @param slice how many processes should be started for each random forest clustering (default = 30)
#' @param bestColname the column name to store the results in
#' @param k the numer of expected clusters (metter more than to view)
#' @param subset how many cells should be randomly selected for the unsupervised clustering (default = 200)
#' @param summaryCol the column storing the summary of the random forest clusterings (default= 'Combined_Group')
#' @param usefulCol the column where all summaryCol groups with less than 10 cells have been merged into the group 0 (default= 'Usefull_Grouping')
#' @return a SingleCellsNGS object including the results and storing the RF object in the usedObj list (bestColname)
#' @export 
setGeneric('rfCluster',
		function ( x, rep=5, SGE=F, email, k=16, slice=30, subset=200, summaryCol='Combined_Group', usefulCol='Usefull_Grouping'){
			standardGeneric('rfCluster')
		}
)
setMethod('rfCluster', signature = c ('SingleCellsNGS'),
		definition = function ( x, rep=5, SGE=F, email, k=16, slice=30, subset=200, summaryCol='Combined_Group', usefulCol='Usefull_Grouping' ) {
			opath = paste( x@outpath,"/RFclust.mp/",sep='' )
			n= paste(x@name, 'RFclust',sep='_')
			m <- max(k)
			if ( length(x@usedObj[['rfExpressionSets']]) == 0 ){
				## start the calculations!
				system( paste('rm -Rf',opath) )
				x@usedObj[['rfExpressionSets']] <- list()
				x@usedObj[['rfObj']] <- list()
				total <- ncol(x@data)
				for ( i in 1:rep) {
					name = paste(n,i,sep='_')
					x@usedObj[['rfExpressionSets']][[ i ]] <- drop.samples( x, colnames(x@data)[sample(1:total,total-subset)], name )
					x@usedObj[['rfObj']][[ i ]] <- RFclust.SGE ( dat=x@usedObj[['rfExpressionSets']][[ i ]]@data, SGE=SGE, slice=30, email=email, tmp.path=opath, name= name )
					x@usedObj[['rfObj']][[ i ]] <- runRFclust ( x@usedObj[['rfObj']][[ i ]] , nforest=500, ntree=500, name=name )
				}
				print ( "You should wait some time now to let the calculation finish! check: system('qstat -f') - re-run the function")
			}
			else {
				for ( i in 1:rep) {
					name = paste(n,i,sep='_')
					## read in the results
					x@usedObj[['rfObj']][[ i ]] <- runRFclust ( x@usedObj[['rfObj']][[ i]] , name=paste(n,i,sep='_') )
					if ( ! is.null(x@usedObj[['rfObj']][[ i ]]@RFfiles[[name]]) ){
						stop( "please re-run this function later - the clustring process has not finished!")
					}
				}
				browser()
				for ( i in 1:rep ) {
					name = paste(n,i,sep='_')
					groups <- createGroups( x@usedObj[['rfObj']][[i]], k=k, name=name )
					x@usedObj[['rfExpressionSets']][[i]]@samples <- cbind ( x@usedObj[['rfExpressionSets']][[i]]@samples, groups[,3:(2+length(k))] )
					## create the required RF object
					x@usedObj[['rfExpressionSets']][[i]] <- bestGrouping( x@usedObj[['rfExpressionSets']][[i]], group=paste('group n=', m) )
					x@samples[, paste( 'RFgrouping', i) ] <-
							predict( x@usedObj[['rfExpressionSets']][[i]]@usedObj[[1]], t(as.matrix(x@data)) )
					print ( paste("Done with cluster",i))
				}
				x@samples[,summaryCol ] <- apply( x@samples[, c( paste('RFgrouping', 1:rep))],1,function (x ) { paste( x, collapse=' ') } )
				useful_groups <- names( which(table( x@samples[,summaryCol ] ) > 10 ))
				x@samples[,usefulCol] <- x@samples[,summaryCol ]
				x@samples[is.na(match ( x@samples[,summaryCol], unique(useful_groups)))==T,usefulCol] <- 'gr. 0'
			}
			x		
		}
)


