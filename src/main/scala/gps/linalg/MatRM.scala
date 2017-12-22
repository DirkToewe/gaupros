package gps.linalg

/**
  * Created by Dirk Toewe on 30.07.17.
  */
class MatRM private[linalg](
  override val nRows: Int,
  override val nCols: Int,
  private[linalg] val vals: Vec
) extends Mat
{
  assert( vals.length == nRows*nCols )

  override def T: MatCM = new MatCM(nCols,nRows,vals)

  override def copy = new MatRM(nRows,nCols,vals.clone)

  override def apply( row: Int, col: Int ) =
  {
    if( row < 0 ) throw new IndexOutOfBoundsException
    if( col < 0 ) throw new IndexOutOfBoundsException
    if( row >= nRows ) throw new IndexOutOfBoundsException
    if( col >= nCols ) throw new IndexOutOfBoundsException
    vals( nCols*row + col )
  }

  override def update( row: Int, col: Int, value: Double ) =
  {
    if( row < 0 ) throw new IndexOutOfBoundsException
    if( col < 0 ) throw new IndexOutOfBoundsException
    if( row >= nRows ) throw new IndexOutOfBoundsException
    if( col >= nCols ) throw new IndexOutOfBoundsException
    vals( nCols*row + col ) = value
  }

  override def modify( row: Int, col: Int, modifier: Double => Double ) =
  {
    if( row < 0 ) throw new IndexOutOfBoundsException
    if( col < 0 ) throw new IndexOutOfBoundsException
    if( row >= nRows ) throw new IndexOutOfBoundsException
    if( col >= nCols ) throw new IndexOutOfBoundsException
    val i = nCols*row + col
    vals(i) = modifier{vals(i)}
  }

  override def foreach[@specialized U]( consumer: Mapper[U] ) =
  {
    var i = vals.length
    var r = nRows
    while( r > 0 )
    {
      r -= 1
      var c = nCols
      while( c > 0 )
      {
        c -= 1; i -= 1
        consumer(vals(i),r,c)
      }
    }
  }

  @inline def diagReduce( reducer: (Double,Double) => Double ): Double =
  {
    assert( nRows == nCols )
    var i = nRows*nRows - 1
    var result = vals(i)
    while( i > 0 )
    {
      i -= nRows+1
      result = reducer( result, vals(i) )
    }
    result
  }

  @inline override def mapReduce[M <: R, @specialized R]( mapper: Mapper[M] )( reducer: (R,R) => R ): R =
  {
    var i = vals.length-1
    var r = nRows
    var c = nCols-1
    var result: R = mapper(vals(i), r-1, c)
    while( r > 0 )
    {
      r -= 1
      while( c > 0 )
      {
        c -= 1; i -= 1
        result = reducer( result, mapper(vals(i), r, c) )
      }
      c = nCols
    }
    result
  }

  override def *# ( mat: Mat ) =
  {
    if( mat.nRows != nCols )
      throw new IllegalArgumentException("Matrix shapes do not match.")
    val lenL=     nRows *     nCols
    val lenR= mat.nRows * mat.nCols
    val result = Vec.zeros(nRows * mat.nCols)
    mat match {
      case MatRM(_,rCols,rVals) =>
        var rStart = result.length
        var i = lenL
        while( rStart > 0 )
        {
          val rEnd = rStart; rStart -= rCols
          var j = lenR
          while( j > 0 )
          {
            i -= 1
            var k = rEnd
            while( k > rStart )
            {
              j -= 1; k -= 1
              result(k) += vals(i) * rVals(j)
            }
          }
        }
      case MatCM(_,_,rVals) =>
        var rStart = lenL
        var i = result.length
        while( rStart > 0 )
        {
          val rEnd = rStart; rStart -= nCols
          var j = lenR
          while( j > 0 )
          {
            i -= 1
            var k = rEnd
            while( k > rStart )
            {
              j -= 1; k -= 1
              result(i) += vals(k) * rVals(j)
            }
          }
        }
    }
    new MatRM(nRows,mat.nCols,result)
  }

  override def *# ( vec: Vec ): Vec =
  {
    if( vec.length != nCols )
      throw new IllegalArgumentException
    val result = Vec.zeros(nRows)
    var i =result.length
    var j =  vals.length
    while( i > 0 )
    {
      i -= 1
      var k = vec.length
      while( k > 0 )
      {
        j -= 1; k -= 1
        result(i) += vec(k) * vals(j)
      }
    }
    result
  }

  def += ( mat: Mat ): Unit =
  {
    if( mat.nRows != nRows ) throw new IllegalArgumentException
    if( mat.nCols != nCols ) throw new IllegalArgumentException
    mat match {
      case MatRM(_,_,summands) => vals += summands
    }
  }
}
object MatRM extends MatFactory
{
  override def zeros( nRows: Int, nCols: Int ): MatRM =
  {
    if( nRows < 0 ) throw new IllegalArgumentException
    if( nCols < 0 ) throw new IllegalArgumentException
    if( nRows*nCols / nCols != nRows ) throw new IllegalArgumentException
    if( nCols*nRows / nRows != nCols ) throw new IllegalArgumentException
    val vals = Vec.zeros(nRows*nCols)
    new MatRM(nRows,nCols,vals)
  }

  override def tabulate( nRows: Int, nCols: Int )( tabulator: (Int, Int) => Double ): Mat
    = tabulateRM(nRows,nCols)( i => j => tabulator(i,j) )

  def tabulateRM( nRows: Int, nCols: Int )( tabulator: Int => Int => Double ): Mat =
  {
    if( nRows < 0 ) throw new IllegalArgumentException
    if( nCols < 0 ) throw new IllegalArgumentException
    if( nRows*nCols / nCols != nRows ) throw new IllegalArgumentException
    if( nCols*nRows / nRows != nCols ) throw new IllegalArgumentException
    val vals = Vec.zeros(nRows*nCols)
    var i = vals.length
    var r = nRows
    while( r > 0 )
    {
      r -= 1
      val tab = tabulator(r)
      var c = nCols
      while( c > 0 )
      {
        c -= 1; i -= 1
        vals(i) = tab(c)
      }
    }
    new MatRM(nRows,nCols,vals)
  }

  override def tabulateSym( size: Int )( tabulator: (Int, Int) => Double ): MatRM
    = tabulateSymRM(size)( i => j => tabulator(i,j) )

  def tabulateSymRM( size: Int )( tabulator: Int => Int => Double ): MatRM =
  {
    if( size < 0 ) throw new IllegalArgumentException
    if( size*size / size != size ) throw new IllegalArgumentException
    val vals = Vec.zeros(size*size)
    var i = vals.length
    var r = size
    while( r > 0 )
    {
      r -= 1
      val tab = tabulator(r)
      var c = size
      while( c > 0 )
      {
        c -= 1; i -= 1
        vals(i) = if( r < c ) vals( r + c*size ) else tab(c)
      }
    }
    new MatRM(size,size,vals)
  }

  private[linalg] def unapply( matRM: MatRM )
    = Some(matRM.nRows, matRM.nCols, matRM.vals)
}