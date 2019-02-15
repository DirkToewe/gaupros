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

/** A super-trait for lower triangular matrix implementations. The main purpose
  * of this implementation is to solve linear equation system with symmetric
  * positive definite matrices. The solver of choice is of course the Cholesky
  * decomposition.
  *
  * The methods of LMat have different prefixes depending on how the data of
  * the triangular matrix is interpreted.
  *
  * <p> If the triangular matrix is interpreted as a lower triangular matrix with
  *     the upper right part of the matrix being filled with zeros, the methods
  *     have <b>tri</b> as a prefix or no prefix at all.
  *
  * <p> If the triangular matrix is interpreted as the lower triangle matrix of a
  *     symmetric matrix, the methods are prefixed by <b>sym</b>.
  *
  * <p> If the triangular matrix is interpreted as the lower triangular matrix result
  *     of a Cholesky decomposition, the methods are prefixed by <b>cholesky</b> with
  *     the exception of `choleskyDecomp` which assumes the matrix to be the lower
  *     triangle of a symmetric positive definite matrix.
  */
trait LMat extends ( (Int,Int) => Double )
{
  import LMat._

  /** The size of this quadratic matrix. size = #rows = #columns.
    */
  def size: Int

  /** Returns the entry value from the (row+1)-th row and (col+1)-th column
    * from this matrix. For the upper lesser triangle, 0 is returned.
    *
    * @param row    The row index.
    * @param col The column index.
    * @return Entry value in the (row+1)-th row and (col+1)-th column.
    */
  @throws[IndexOutOfBoundsException]
  override def apply( row: Int, col: Int ): Double

  /** Assuming this to be the lower triangle of a symmetric matrix A, this matrix
    * returns the entry value A(row,col).
    *
    * @param row The row of the entry that is to be read.
    * @param col The column of the antry that is to be read.
    * @return col > row ? this(col,row) : this(row,col)
    */
  def symGet( row: Int, col: Int ): Double
    = if( col > row )
        this(col,row)
      else
        this(row,col)

  /** Changes the entry value from the (row+1)-th row and (col+1)-th column
    * from this matrix. For the upper lesser triangle, an IndexOutOfBoundsException
    * is thrown.
    *
    * @param row    The row index.
    * @param col The column index.
    * @throws IndexOutOfBoundsException If col > row or either of the indices
    *                                   is out of matrix bounds.
    */
  @throws[IndexOutOfBoundsException]
  def update( row: Int, col: Int, value: Double ): Unit

  /** Changes the entry value from the (row+1)-th row and (col+1)-th column
    * from this matrix using the given function to map the entry value to
    * it's new value. For the upper lesser triangle, an IndexOutOfBoundsException
    * is thrown.
    *
    * @param row    The row index.
    * @param col The column index.
    * @throws IndexOutOfBoundsException If col > row or either of the indices
    *                                   is out of matrix bounds.
    */
  @throws[IndexOutOfBoundsException]
  def modify( row: Int, col: Int, modifier: Double => Double ): Unit

  /** Applies each entry of the lower triangular matrix to the function and writes
    * the return value back as that entry's value.
    *
    * @param modifier The function used to modify each entry in the lower triangle
    *                 of this matrix.
    */
  def modify( modifier: Mapper[Double] ): Unit

//  def foreach[@specialized U]( consumer: Double => U ): Unit

  def foreach[@specialized U]( consumer: Mapper[U] ): Unit

  @inline def mapReduce[M <: R, @specialized R]( mapper: Mapper[M] )( reducer: (R,R) => R ): R

  /** Creates a deep clone of this lower triangular matrix.
    *
    * @return A deep copy of this lower triangular matrix.
    */
  def copy: LMat

  def unary_- : LMat

  /** Returns the determinant of this lower triangular matrix, which is the
    * product of its diagonal values.
    */
  def det: Double

  /** Solves <b>x</b> in <b>this ⋅ x = y</b> in-place, overwriting the right
    * hand side vector in the process. Since this matrix is lower triangular
    * forward substitution is used to solve the equation in
    * <b>O(size<sup>2</sup>)</b> time.
    *
    * @param y The right hand side vector of the equation. The solution is
    *          written back into that same vector.
    */
  def triSolve( y: Vec ): Unit

  /** Solves <b>x</b> in <b>this ⋅ x = y</b> in-place, overwriting the right
    * hand side vector in the process. Since this matrix is lower triangular
    * forward substitution is used to solve the equation.
    *
    * @param y The right hand side matrix of the equation. The solution is
    *          written back into that same matrix.
    */
  def triSolve( y: Mat ): Unit

  /** With solution == this.triSolve(y), the method computes a weight matrix W
    * with Sum[i,j]( W[i,j]*dThis[i,j]/dp ) == d( solution ⋅ weights )/dp under the
    * condition that y and weights is independent of p. <b style="color:#FF0000">The
    * contents of this matrix are overwritten with W.</b>
    *
    * @param solution The solution vector of which a weighted sum is to be derivated.
    * @param weights The weight vector of the weighted sum that is to be derivated.
    *                <b style="color:#FF0000">Content of weights is overwritten in
    *                the process, caching information temporarily.</b>
    */
  @inline def triSolveSumDerive(solution: Vec, weights: Vec ): Unit

  /** Calculates in-place the inverse of this lower triangular matrix.
    */
  def triInvert(): Unit

  /** Returns the diagonal of this matrix.
    *
    * @return The diagonal of this matrix from top left to bottom right.
    */
  def diag(): Vec

  /** Goes through the diagonal entries and calls the consumer function
    * for each of them.
    *
    * @param consumer The callback to be called for each diagonal entry.
    * @tparam U consumer's return type. The return value is ignored.
    */
  @inline def diagForeach[@specialized U]( consumer: Double => U )

  /** Goes through the diagonal entries and calls the consumer function
    * for each of them.
    *
    * @param consumer The callback to be called for each diagonal entry.
    *                 The second argument is the entry index.
    * @tparam U consumer's return type. The return value is ignored.
    */
  @inline def diagForeach[@specialized U]( consumer: (Double,Int) => U )

  @inline def diagModify( modfier: (Double,Int) => Double )

//  def diagReduce( reductor: (Double,Double) => Double )

  /** Performs an in-place cholesky decomposition on this matrix. This method
    * assumes this matrix to be the lower triangle of a symmetric matrix.
    */
  def choleskyDecomp(): Unit

  /** Assuming this matrix is the cholesky decomposition of a matrix A, i.e. A = this *# this<sup>T<sup>,
    * this method computes a lower triangular  weight matrix W such that:
    * d( sum[i,j]( weights[i,j]*this[i,j] )/dp = sum[i,j]( W[i,j]*(dA/dp)[i,j] ).
    *
    * @param weights The weights of the weighted sum that is to be derived. <b style="color:#FF0000">
    *                Contents of weights is overwritten with the result W in the process.</b>
    */
  def choleskySumDerive( weights: LMat ): Unit

  /** Solves <b>x</b> in <b>this ⋅ this<sup>T</sup> ⋅ x = y</b> and writes the
    * result back into y.
    *
    * @param y The right hand side of the linear equation system to be solved.
    *          The array's content is overwritten with the solution.
    */
  def choleskySolve( y: Vec ): Unit

  /** Solves <b>x</b> in <b>this ⋅ this<sup>T</sup> ⋅ x = y</b> and writes the
    * result back into y.
    *
    * @param y The right hand side of the linear equation system to be solved.
    *          The array's content is overwritten with the solution.
    */
  def choleskySolve( y: Mat ): Unit

  /** Assuming this matrix represents the lower triangular result of a Cholesky
    * decomposition, this method returns the determinant of the decomposed matrix.
    *
    * @return <p>det(this ⋅ this<sup>T</sup>) = det(this) * det(this<sup>T</sup>) = det(this)<sup>2</sup>
    */
  final def choleskyDet(): Double =
  {
    // det(L @ L.T) = det(L) * det(L.T) = det(L)²
    val det = this.det
    return det*det
  }

  /** Assuming this matrix represents the lower triagular result of a Cholesky,
    * this method compute the inverse of this matrix in-place.
    */
  def choleskyInvert(): Unit

  /** Assuming this matrix is the lower triangle of a symmetric matrix, the method
    * goes through all elements of the specified row of said matrix and calls
    * a consumer each time.
    *
    * @param row
    * @param consumer
    * @tparam U
    */
  def symRowForeach[@specialized U]( row: Int )( consumer: (Double,Int) => U )

  /** Assuming <b>this</b> is the lower triangle of a symmetric matrix <b>A</b>
    * and <b>b</b> returns the lower triangle of a symmetric matrix <b>B</b>,
    * this methods computes the trace of the matrix product <b>A⋅B</b>.

    * @param b <b>b(i,j)</b> returns <b>B<sub>ij</sub></b> for <b>0 ≤ j ≤ i < a.size</b>
    * @return <b>tr(this⋅B)</b>
    */
  def symMatProdTrace( b: (Int,Int) => Double ): Double

  /** Assuming this is the lower triangle of a symmetric matrix A, returns the lower triangle of the symmetric
    * matrix A⋅A.
    *
    * @return
    */
  def symSquare: LMat

  /** Assuming this is the lower triangle of a symmetric matrix A, this method returns the lower triangle b
    * of a matrix B, where B<sub>i,j</sub> = ∑<sub>k</sub> product(a<sub>i,k</sub>, a<sub>k,j</sub>, i, k, j).
    * The matrix product A⋅A could for example be computed using symMatProdSquareLike(a)( (a_ik,a_kj, i,k,j) => a_ik*a_kj ).
    *
    * @param product
    * @return
    */
  @inline def symSquareLike( product: ProductLike ): LMat

  /** Assuming this is the lower triangle of a symmetric matrix A, this method
    * computes product of A and the given vector.
    *
    * @param vec The right factor of the matrix-vector product.
    * @return <b>A⋅vec</b>
    */
  def symMatMul( vec: Vec ) =
  {
    val result = Vec.zeros(vec.length)
    foreach{
      (Aij,i,j) => result(i) += vec(j) * Aij
        if(i != j) result(j) += vec(i) * Aij
    }
    result
  }

  /** Assuming this matrix is the lower triangle of a symmetric matrix, this
    * methods returns a string representation of said triangular matrix.
    *
    * @return A string representation of this matrix interpreted as a symmetric
    *         matrix's lower triangle.
    */
  def symToString
    = 0 until size map (
        row => 0 until size map {
          col => f"${symGet(row,col)}%12f"
        } mkString ","
      ) mkString ("[ [", " ],\n  [", " ] ]")

  override def toString
    = 0 until size map {
        row => 0 until size map (
          col => f"${this(row,col)}%12f"
        ) mkString ","
      } mkString ("[ [", " ],\n  [", " ] ]")
}
object LMat extends LMatFactory
{
  // FIXME that is not the correct MAX_SIZE, which would be 65535. But using said MAX_SIZE requires changes in the implementation code
  // val MAX_SIZE = 46340
  val MAX_SIZE = 65535

  /** Returns the number of entries in a size*size triangular matrix.
    *
    * @param size The size of the triangular matrix.
    * @return size*(size+1)/2
    */
  def nEntries( size: Int ): Int =
  {
    assert( size <= MAX_SIZE )
    if (size % 2 == 0)
      size / 2 * (size + 1)
    else
      (size + 1) / 2 * size
  }

  override def zeros(size: Int): LMat
    = LMatCM.zeros(size)

  override def tabulate( size: Int )( idx2val: (Int,Int) => Double ): LMat
    = LMatCM.tabulate(size)(idx2val)

  trait ProductLike {
    @inline def apply( A_ik: Double, A_kj: Double, i: Int, k: Int, j: Int ): Double
  }

  /** Given the lower triangle of a symmetric matrix A and the vectors u and v,
    * this matrix computes the inner product <b>u<sup>T</sup>Av</b>
    *
    * @param u The left factor.
    * @param A The lower triangle of the symmetric matrix A, which is the base of the inner product.
    * @param v The right factor.
    * @return <b>u<sup>T</sup>Av</b>
    */
  def symInner( u: Vec, A: LMat, v: Vec ): Double =
  {
    if( u.length != A.size ) throw new IllegalArgumentException
    if( v.length != A.size ) throw new IllegalArgumentException
    A.mapReduce{
      (Aij,i,j) =>
        var        uv  = u(i)*v(j)
        if(i != j) uv += u(j)*v(i)
        Aij * uv
    }(_ + _)
  }
}