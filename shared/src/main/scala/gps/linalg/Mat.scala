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
  * Created by dtitx on 29.07.17.
  */
trait Mat extends ( (Int,Int) => Double )
{
  assert( nRows*nCols / nCols == nRows )
  assert( nCols*nRows / nRows == nCols )

  /** Returns the number of rows in this matrix.
    * @return the number of rows in this matrix.
    */
  def nRows: Int
  /** Returns the number of columns in this matrix.
    * @return the number of columns in this matrix.
    */
  def nCols: Int

  def update( row: Int, col: Int, value: Double ): Unit

  def modify( row: Int, col: Int, modifier: Double => Double )

  def foreach[@specialized U]( consumer: Mapper[U] ): Unit

  def trace = diagReduce(_ + _)

  @inline def diagReduce( reducer: (Double,Double) => Double ): Double

  def copy: Mat

  @inline def mapReduce[M <: R, @specialized R]( mapper: Mapper[M] )( reducer: (R,R) => R ): R

  /** Returns the matrix product of this matrix and the matrix given.
    * The asymmetry of the infix operator is supposed to remind of its
    * non-commutativity.
    *
    * @param mat The right factor of the matrix product.
    * @return this times mat
    */
  def *# ( mat: Mat ): Mat

  def *# ( vec: Vec ): Vec

  def += ( mat: Mat ): Unit

  /** Returns a transpose <b>view</b> of this matrix.
    *
    * @return a transpose <b>view</b> of this matrix.
    */
  def T: Mat

  // def luDecomp( permutations: Array[Int] ): Unit
  // def luSolve( permutations: Array[Int], rhs: Vec ): Unit

  override def toString: String
    = 0 until nRows map {
        row => 0 until nCols map (
          col => f"${this(row,col)}%12f"
          ) mkString ","
      } mkString ("[ [", " ],\n  [", " ] ]")
}
object Mat extends MatFactory
{
  override def zeros( nRows: Int, nCols: Int ): Mat
    = MatRM.zeros(nRows,nCols)

  override def tabulate( nRows: Int, nCols: Int )( tabulator: (Int, Int) => Double )
    = MatRM.tabulate(nRows,nCols)(tabulator)

  override def tabulateSym( size: Int )( tabulator: (Int, Int) => Double )
    = MatRM.tabulateSym(size)(tabulator)

  /** Returns the inner(-like) product of a vector over the given matrix,
    * i.e. the scalar result of <b>u<sup>T</sup> ⋅ A ⋅ v</b>.
    *
    * @param u
    * @param A
    * @param v
    * @return <b>u<sup>T</sup> ⋅ A ⋅ v</b>
    */
  def inner( u: Vec, A: Mat, v: Vec ): Double =
  {
    if( u.length != A.nRows ) throw new IllegalArgumentException
    if( v.length != A.nCols ) throw new IllegalArgumentException
    A.mapReduce{ _ * u(_) * v(_) }(_ + _)
  }
}