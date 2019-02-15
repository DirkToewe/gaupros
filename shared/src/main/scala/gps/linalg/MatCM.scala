/* This file is part of GauProS.
 *
 * GauProS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GauProS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GauProS.  If not, see <https://www.gnu.org/licenses/>.
 */

package gps.linalg

/**
  * Created by Dirk Toewe on 30.07.17.
  */
class MatCM private[linalg](
  override val nRows: Int,
  override val nCols: Int,
  private[linalg] val vals: Vec
) extends Mat
{
  assert( vals.length == nRows*nCols )

  override def T: MatRM = new MatRM(nCols,nRows,vals)

  override def copy = new MatCM(nRows,nCols,vals.clone)

  override def apply( row: Int, col: Int ) =
  {
    if( row < 0 ) throw new IndexOutOfBoundsException
    if( col < 0 ) throw new IndexOutOfBoundsException
    if( row >= nRows ) throw new IndexOutOfBoundsException
    if( col >= nCols ) throw new IndexOutOfBoundsException
    vals( row + col*nRows )
  }

  override def update( row: Int, col: Int, value: Double ) =
  {
    if( row < 0 ) throw new IndexOutOfBoundsException
    if( col < 0 ) throw new IndexOutOfBoundsException
    if( row >= nRows ) throw new IndexOutOfBoundsException
    if( col >= nCols ) throw new IndexOutOfBoundsException
    vals( row + col*nRows ) = value
  }

  override def modify( row: Int, col: Int, modifier: Double => Double ) =
  {
    if( row < 0 ) throw new IndexOutOfBoundsException
    if( col < 0 ) throw new IndexOutOfBoundsException
    if( row >= nRows ) throw new IndexOutOfBoundsException
    if( col >= nCols ) throw new IndexOutOfBoundsException
    val i = row + col*nRows
    vals(i) = modifier{vals(i)}
  }

  override def foreach[@specialized U]( consumer: Mapper[U] ) =
  {
    var i = vals.length
    var c = nCols
    while( c > 0 )
    {
      c -= 1
      var r = nRows
      while( r > 0 )
      {
        r -= 1; i -= 1
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
    var r = nRows-1
    var c = nCols
    var result: R = mapper(vals(i), r, c-1)
    while( c > 0 )
    {
      c -= 1
      while( r > 0 )
      {
        r -= 1; i -= 1
        result = reducer( result, mapper(vals(i), r, c) )
      }
      r = nRows
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
        var cStart = lenL
        var i = lenR
        while( cStart > 0 )
        {
          val cEnd = cStart; cStart -= nRows
          var j = result.length
          while( j > 0 )
          {
            i -= 1
            var k = cEnd
            while( k > cStart )
            {
              j -= 1; k -= 1
              result(j) += vals(k) * rVals(i)
            }
          }
        }
      case MatCM(_,_,rVals) =>
        var cStart = result.length
        var i = lenR
        while( cStart > 0 )
        {
          val cEnd = cStart; cStart -= nRows
          var j = lenL
          while( j > 0 )
          {
            i -= 1
            var k = cEnd
            while( k > cStart )
            {
              j -= 1; k -= 1
              result(k) += vals(j) * rVals(i)
            }
          }
        }
    }
    new MatCM(nRows,mat.nCols,result)
  }

  override def *# ( vec: Vec ): Vec =
  {
    if( vec.length != nCols )
      throw new IllegalArgumentException
    val result = Vec.zeros(nRows)
    var j =vals.length
    var i = vec.length
    while( i > 0 )
    {
      i -= 1
      var k = result.length
      while( k > 0 )
      {
        j -= 1; k -= 1
        result(k) += vals(j) * vec(i)
      }
    }
    result
  }

  def += ( mat: Mat ): Unit =
  {
    if( mat.nRows != nRows ) throw new IllegalArgumentException
    if( mat.nCols != nCols ) throw new IllegalArgumentException
    mat match {
      case MatCM(_,_,summands) => vals += summands
    }
  }
}
object MatCM extends MatFactory
{
  override def zeros( nRows: Int, nCols: Int ): MatCM =
  {
    if( nRows < 0 ) throw new IllegalArgumentException
    if( nCols < 0 ) throw new IllegalArgumentException
    if( nRows*nCols / nCols != nRows ) throw new IllegalArgumentException
    if( nCols*nRows / nRows != nCols ) throw new IllegalArgumentException
    val vals = Vec.zeros(nRows*nCols)
    new MatCM(nRows,nCols,vals)
  }

  override def tabulate( nRows: Int, nCols: Int )( tabulator: (Int, Int) => Double ): MatCM
    = tabulateCM(nRows,nCols)( j => i => tabulator(i,j) )

  def tabulateCM( nRows: Int, nCols: Int )( tabulator: Int => Int => Double ): MatCM =
  {
    if( nRows < 0 ) throw new IllegalArgumentException
    if( nCols < 0 ) throw new IllegalArgumentException
    if( nRows*nCols / nCols != nRows ) throw new IllegalArgumentException
    if( nCols*nRows / nRows != nCols ) throw new IllegalArgumentException
    val vals = Vec.zeros(nRows*nCols)
    var i = vals.length
    var c = nCols
    while( c > 0 )
    {
      c -= 1
      val tab = tabulator(c)
      var r = nRows
      while( r > 0 )
      {
        r -= 1; i -= 1
        vals(i) = tab(r)
      }
    }
    new MatCM(nRows,nCols,vals)
  }

  override def tabulateSym( size: Int )( tabulator: (Int,Int) => Double ): MatCM
    = tabulateSymCM(size){ j => i => tabulator(i,j) }

  def tabulateSymCM( size: Int )( tabulator: Int => Int => Double ): MatCM =
  {
    if( size < 0 ) throw new IllegalArgumentException
    if( size*size / size != size ) throw new IllegalArgumentException
    val vals = Vec.zeros(size*size)
    var i = vals.length
    var c = size
    while( c > 0 )
    {
      c -= 1
      val tab = tabulator(c)
      var r = size
      while( r > 0 )
      {
        r -= 1; i -= 1
        vals(i) = if( r > c )  vals( size*r + c ) else tab(r)
      }
    }
    new MatCM(size,size,vals)
  }

  private[linalg] def unapply( matCM: MatCM )
    = Some(matCM.nRows, matCM.nCols, matCM.vals)
}