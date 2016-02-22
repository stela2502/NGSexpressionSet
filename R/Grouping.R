

#' @name bestGrouping
#' @aliases 'bestGrouping,SingleCellsNGS-method
#' @rdname 'bestGrouping-methods
#' @docType methods
#' @description The function is using a randomForest classifier with 2000 trees to classify the given data using the given grooping
#' @description All groups that fail to be prediceted using the random forest are deemed ungrouped.
#' @description All groups where less than 50 percent of the total samples geting classified as being from that group fail.
#' @param x the single cells ngs object
#' @param groups a vector of sample columns that should be checked (the most complex is used only)
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
			rf <- randomForest( x= t(as.matrix(x@data)), y=x@samples[, group],ntree=2000 )
			t <- table( observed= x@samples[,group ], predicted = rf$predicted )
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
			x@samples[,bestColname] <- as.numeric(x@samples[, group])
			for ( b in BAD ) {
				x@samples$QualifiedGrouping[ which(x@samples[,bestColname] == b)] <- 0
			}
			
			for (i in 0:(length(table(x@samples[,bestColname]))-1)){
				modify <- which(x@samples[,bestColname] >= i )
				if ( length(modify) == 0 ) { break}
				while( length(which(x@samples[modify,bestColname] == i)) == 0 ){
					x@samples[modify,bestColname] = x@samples[modify, bestColname] -1
				}
			}
			x
		}
)