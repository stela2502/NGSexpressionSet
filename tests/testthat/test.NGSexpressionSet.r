library('NGSexpressionSet')

PMID25158935samples <- PMID25158935samples[ - grep ('ERR420370', PMID25158935samples[,'Sample'] ) ,]
colnames(PMID25158935samples)[20] <- 'bam filename'
PMID25158935 <- NGSexpressionSet( PMID25158935exp, PMID25158935samples,  Analysis = NULL, name='PMID25158935', namecol='Sample', namerow= 'GeneID', usecol=NULL , outpath = NULL)

expect_equal(class(PMID25158935), c('NGSexpressionSet', 'ExpressionSet') )

expect_equal(PMID25158935$outpath, pwd())


#PMID25158935 <- createStats(PMID25158935, condition='GroupName', A='HSC', B='MPP1' )
PMID25158935$stats <- list( 'HSC vs. MPP1' = HSC_MPP1 )

normalized <- normalize(PMID25158935)
expect_equal(round(as.vector(t(normalized$data[4,]))), c(21773, 28980, 17473, 22466, 24822, 31913, 24797, 25424, 26666, 29825, 23307, 27206, 25841, 27314 ,25949))

PMID25158935 <- simpleAnova( PMID25158935 )

