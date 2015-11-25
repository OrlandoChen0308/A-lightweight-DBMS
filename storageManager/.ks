ff = function( experiment, query, element="z", keep.scores=FALSE) {
  ## rank data matrix in descending order
  data.matrix <- as( assayDataElement(experiment, element), "matrix" )

  ## subset objects to shared genes
  matched.features <- match( featureNames( experiment ), featureNames(query))
  matched.sets <- query[na.omit(matched.features),]

  ## extract scores for each gene set
  sets.up     <- mclapply( seq(ncol(matched.sets)),
    function( x ) which(members( matched.sets )[ ,x ] == 1 ))

  sets.down     <- mclapply( seq(ncol(matched.sets)),
    function( x ) which(members( matched.sets )[ ,x ] == -1 ))

  ## transform experiment to (reverse) ranks
  rank.matrix <- apply(data.matrix, 2, function(x) { length(x) - rank(x) + 1 } )

  ## calculate connectivity score
  raw.score <- apply( rank.matrix, 2, function( r ) {
    sapply(seq_along( sets.up ), function( n ) {
      .s( r[sets.up[[n]]], r[sets.down[[n]]], length( r ) )
})
    })

  print(raw.score)
}

.ks <- function( V, n ) {
  t <- length( V )

  if( t == 0 )  {
    return( 0 )
  } else {

    if ( is.unsorted( V ) )
      V <- sort( V )
    d <- (1:t) / t - V / n
    a <- max( d )
    b <- -min( d ) + 1 / t
    ifelse( a > b, a, -b )
  }
}

.s <- function( V_up, V_down, n ) {
  ks_up <- .ks( V_up, n )
  ks_down <- .ks( V_down, n )
  ifelse( sign( ks_up ) == sign( ks_down ), 0, ks_up - ks_down )
}

