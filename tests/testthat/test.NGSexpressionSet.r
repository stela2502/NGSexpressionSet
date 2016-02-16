#library('NGSexpressionSet')

PMID25158935 <- NGSexpressionSet( PMID25158935exp, PMID25158935samples, Analysis = NULL, name='PMID25158935', namecol='Sample',		namerow= 'GeneID', usecol=NULL , outpath = '')

expect_equal(class(PMID25158935)[1], 'NGSexpressionSet' )

red2 <- reduce.Obj ( PMID25158935, rownames(PMID25158935@data)[1:100], name='minimal' )
expect_true(  all.equal( red2, red) , 'saved data OK' )


expect_equal( class(red)[1], 'NGSexpressionSet' )
expect_equal( dim(red@data), c(100,15) )
expect_equal( red@name, 'minimal' )

dropS <- drop.samples( red, samplenames=red@samples[1:5, red@sampleNamesCol])
expect_equal( dim(dropS@data), c(100,10) )
expect_equal( dim(dropS@samples), c( 10,21) )
expect_equal( dim(dropS@annotation), c( 100,2) )

expect_equal(PMID25158935@outpath, pwd())


#PMID25158935 <- createStats(PMID25158935, condition='GroupName', A='HSC', B='MPP1' )
PMID25158935$stats <- list( 'HSC vs. MPP1' = HSC_MPP1 )

normalized <- normalize(PMID25158935)
expect_equal(round(as.vector(t(normalized$data[4,]))), c(21773, 28980, 17473, 22466, 24822, 31913, 24797, 25424, 26666, 29825, 23307, 27206, 25841, 27314 ,25949))

PMID25158935 <- simpleAnova( PMID25158935 )

