## here I need to test the grouping strategy!

old_LTHSC <- update.measurements(old_LTHSC)
rfAnalysis <- list()
i <- 1
start <- 0.9
run = FALSE
if ( run ) {
	rfAnalysis[[i]] <- reduce.Obj (old_LTHSC, interesting_genes( old_LTHSC , Max=start, Min=start-0.1, genes=200),name= paste('old_LTHSC_default_test_200_', start,(start-0.1),sep='_'))
	rfAnalysis[[i]]  <- rfCluster ( rfAnalysis[[i]],rep=2, SGE=F,slice=2, k=10, email='stefan.lang@med.lu.se', pics=T, subset=100, nforest=6, ntree=200, name='test1')
}
