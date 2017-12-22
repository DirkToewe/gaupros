package gps.linalg

/**
  * Created by Dirk Toewe on 30.07.17.
  */
trait MatFactory
{
  def zeros( nRows: Int, nCols: Int ): Mat

  def tabulate( nRows: Int, nCols: Int )( tabulator: (Int,Int) => Double ): Mat

  def tabulateSym( size: Int )( tabulator: (Int,Int) => Double ): Mat

  def apply( nRows: Int, nCols: Int )( values: Double* ): Mat =
  {
    val vals: IndexedSeq[Double] = values match {
      case vals: IndexedSeq[Double] => vals
      case _ => values.toArray
    }
    tabulate(nRows,nCols){ (i,j) => vals(nCols*i + j) }
  }
}