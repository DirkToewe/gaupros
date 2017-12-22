package gps.linalg

import LMat.MAX_SIZE

import scala.annotation.tailrec

/** Abstract super trait for lower triangular matrix factories.
  *
  * Created by Dirk Toewe on 19.06.17.
  */
trait LMatFactory
{
  /** Creates a new lower triangular matrix of the specified size using the
    * given function to initialize the entries.
    *
    * @param size    The size of the created matrix, i.e. the number of rows
    *                as well as the number of columns of said matrix.
    * @param idx2val The function used to initialize the lower triangular
    *                entry values.
    * @return        A new lower triangular matrix <code>m</code> with
    *                <code>m<sub>i,j</sub> = idx2val(i,j)</code> for
    *                <code>i≥j</code>.
    */
  def tabulate( size: Int )( idx2val: (Int,Int) => Double ): LMat

  def zeros( size: Int ): LMat

  def apply( values: Double* ): LMat =
  {
    @tailrec
    def binSearch( lo: Int, hi: Int ): Int =
    {
      assert( lo < hi )
      val mid = (lo+hi)/2 ensuring (_ < hi)
      val nVals = mid*(mid+1) / 2
           if( nVals > values.size )                binSearch(lo,mid)
      else if( nVals < values.size ) if( mid > lo ) binSearch(mid,hi)
      else throw new IllegalArgumentException("values do not fill triangular matrix.")
      else mid
    }
    val size = binSearch(-1,MAX_SIZE+1)

    tabulate(size){ (row,col) => values(col + row*(row+1) / 2) }
  }
}
