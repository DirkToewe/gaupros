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