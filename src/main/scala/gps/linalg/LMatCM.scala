package gps.linalg

import math.{sqrt => √}
import LMat.{MAX_SIZE, ProductLike, nEntries}

/** A column major 1D-array-based implementation of a lower triangular matrix.
  *
  * Created by Dirk Toewe on 18.06.17.
  */
final class LMatCM private( override val size: Int, private val vals: Vec ) extends LMat
{
  private def index( row: Int, col: Int ): Int
    = if( col%2 == 0 )
        row + col/2 * (2*size - col-1)
      else
        row + (2*size - col-1)/2 * col

  override def apply( row: Int, col: Int ): Double =
  {
    if( row < 0 ) throw new IndexOutOfBoundsException(   "Negative row index.")
    if( col < 0 ) throw new IndexOutOfBoundsException("Negative column index.")
    if( row >= size ) throw new IndexOutOfBoundsException(   "Row index out of bounds.")
    if( col >= size ) throw new IndexOutOfBoundsException("Column index out of bounds.")
    if( col > row )
      return 0
    vals{index(row,col)}
  }

  override def update( row: Int, col: Int, value: Double ): Unit =
  {
    if( col < 0    ) throw new IndexOutOfBoundsException("Negative column index.")
    if( col > row  ) throw new IndexOutOfBoundsException("Upper right lesser triangle is read-only.")
    if( row >= size) throw new IndexOutOfBoundsException("Row index out of bounds.")
    vals{index(row,col)} = value
  }

  override def modify( row: Int, col: Int, modifier: Double => Double ): Unit =
  {
    if( col < 0    ) throw new IndexOutOfBoundsException("Negative column index.")
    if( col > row  ) throw new IndexOutOfBoundsException("Upper right lesser triangle is read-only.")
    if( row >= size) throw new IndexOutOfBoundsException("Row index out of bounds.")
    val i = index(row,col)
    vals(i) = modifier(vals{i})
  }

  override def modify( modifier: Mapper[Double] ): Unit =
  {
    var i = 0; var c = 0
    while( c < size )
    {
      var r = c
      while( r < size )
      {
        vals(i) = modifier(vals(i),r,c)
        r += 1; i += 1
      }
      c += 1
    }
    // assert( i == vals.size )
  }

  override def foreach[@specialized U]( consumer: Mapper[U] ): Unit =
  {
    var i = 0; var c = 0
    while( c < size )
    {
      var r = c
      while( r < size )
      {
        consumer(vals(i),r,c)
        r += 1; i += 1
      }
      c += 1
    }
  }

  @inline override def mapReduce[M <: R, @specialized R]( mapper: Mapper[M] )
                                                        ( reducer: (R, R) => R ): R =
  {
    if( 0 == size )
      throw new IllegalStateException
    var i = vals.length - 1
    var c = size-1
    var result: R = mapper(vals(i), size-1,c)
    while( c > 0 )
    {
      c -= 1
      var r = size
      while( r > c )
      {
        r -= 1; i -= 1
        result = reducer( result, mapper(vals(i), r,c) )
      }
    }
    result
  }

  override def copy    = new LMatCM(size,  vals.clone)
  override def unary_- = new LMatCM(size, -vals)

  override def det: Double =
  {
    var det = 1.0
    var  i = vals.length-1
    var di = size
    while( di > 0 )
    {
      det *= vals(i)
      i -= di; di -= 1
    }
    det
  }

  override def diag(): Vec  =
  {
    val result = Vec.zeros(size)
    var cSize = 1
    var i = vals.length
    while( i > 0 )
    {
      i -= cSize
      result(size-cSize) = vals{i}
      cSize += 1
    }
    result
  }

  @inline override def diagForeach[@specialized U]( consumer: Double => U ): Unit =
  {
    var cSize = 1
    var i = vals.length
    while( i > 0 )
    {
      i -= cSize; cSize += 1
      consumer(vals{i})
    }
  }

  @inline override def diagForeach[@specialized U]( consumer: (Double,Int) => U ): Unit =
  {
    var cSize = 1
    var i = vals.length
    while( i > 0 )
    {
      i -= cSize
      consumer(vals{i}, size-cSize)
      cSize += 1
    }
  }

  @inline override def diagModify( modifier: (Double,Int) => Double ): Unit =
  {
    var cSize = 1
    var i = vals.length
    while( i > 0 )
    {
      i -= cSize
      vals(i) = modifier(vals{i}, size-cSize)
      cSize += 1
    }
  }

  override def triSolve( rhs: Vec ) =
  {
    if( rhs.length != size )
      throw new IllegalArgumentException
    // forward substitution
    var i = 0; var c = 0
    while( c < size )
    {
      rhs(c) /= vals(i)
      var r = c+1; i += 1
      while( r < size )
      {
        rhs(r) -= rhs(c) * vals(i)
        i += 1; r += 1
      }
      c += 1
    }
  }

  override def triSolve( rhs: Mat ) =
  {
    if( rhs.nRows != size )
      throw new IllegalArgumentException
    val len = rhs.nRows * rhs.nCols
    rhs match {
      case MatRM(_,nCols,rVals) =>
        var i = 0; var rEnd = 0
        while( rEnd < len )
        {
          val rStart = rEnd; rEnd += nCols
          var j = rStart
          while( j < rEnd )
          {
            rVals(j) /= vals(i)
            j += 1
          }
          i += 1
          while( j < len )
          {
            var k = rStart
            while( k < rEnd )
            {
              rVals(j) -= rVals(k) * vals(i)
              j += 1; k += 1
            }
            i += 1
          }
        }
      case MatCM(_,_,rVals) =>
        var i = 0
        while( i < len )
        {
          val cEnd = i+size // <- end of the currently solved column in rhs
          var j = 0
          while( i < cEnd )
          {
            rVals(i) /= vals(j)
            var k = i+1; j += 1
            while( k < cEnd )
            {
              rVals(k) -= rVals(i) * vals(j)
              j += 1; k += 1
            }
            i += 1
          }
        }
    }
  }

  @inline override def triSolveSumDerive(solution: Vec, weights: Vec ): Unit =
  {
    // BACKTRACKING TRI_SOLVE
    //
    // β[i] = (y[i] - K[i,0]*β[0] - K[i,1]*β[1] - ... - K[i,i-1]*β[i-1]) / K[i,i]
    //
    // dβ[i]/dp = - dK[i,i]/dp * β[i] / K[i,i]
    //            - ( dK[i,0]/dp * β[0] + dK[i,1]/dp * β[1] + ... + dK[i,i-1]/dp * β[i-1] ) / K[i,i]
    //            - ( K[i,0] * dβ[0]/dp + K[i,1] * dβ[1]/dp + ... + K[i,i-1] * dβ[i-1]/dp ) / K[i,i]
    //
    var i = vals.length; var r = size
    while( i > 0 )
    {
      i -= 1; r -= 1
      var j = size-1
      while( j > r )
      {
        weights(r) -= weights(j) * vals(i)
        vals(i)  =  - weights(j) * solution(r)
        i -= 1; j -= 1
      }
      weights(r) /= vals(i)
      vals(i) = - weights(j) * solution(r)
    }
  }

  override def triInvert(): Unit =
  {
    var cStart = vals.length
    var cSize = 0
    while( cStart > 0 )
    {
      cStart -= 1
      val cEnd = cStart; cStart -= cSize; cSize += 1

      var i = vals.length-1
      var j = cEnd
      while( j > cStart )
      {
        var k = cEnd
        while( k > j )
        {
          vals(k) -= vals(j)*vals(i)
          i -= 1; k -= 1
        }
        vals(j) *= -vals(i)
        i -= 1; j -= 1
      }

      vals(cStart) = 1 / vals(cStart)
      while( i > cStart )
      {
        vals(i) *= vals(cStart)
        i -= 1
      }
    }
  }

  override def choleskyDecomp(): Unit =
  {
    var i=0
    var cSize = size
    while( cSize > 0 )
    {
    // SOLVE COLUMN
      // diagonal element
      if( ! {+0.0 < vals(i)} ) // handles NaN as well
        throw MatrixNotPositiveDefinite()
      vals(i) = √(vals{i})
      // elements below diagonal
      var j = i+1
      while( j < i+cSize )
      {
        vals(j) /= vals(i)
        j += 1
      }
    // ELIMINATE COLUMN FROM REMAINING COLUMNS
      i += 1
      var cs = cSize-1
      while( cs > 0 )
      {
        var k = i
        while( k < i+cs )
        {
          vals(j) -= vals(i) * vals(k)
          j += 1; k += 1
        }
        i += 1; cs -= 1
      }
      cSize -= 1
    }
  }

  override def choleskySolve( rhs: Vec ): Unit =
  {
    // forward substitution
    triSolve(rhs)
    // backward substitution
    var i = vals.length; var r = size
    while( r > 0 )
    {
      r -= 1
      var c = size-1
      while( c > r )
      {
        i -= 1
        rhs(r) -= rhs(c) * vals(i)
        c -= 1
      }
      i -= 1
      rhs(r) /= vals(i)
    }
  }

  override def choleskySolve( rhs: Mat ): Unit =
  {
    // forward substitution
    triSolve(rhs)
    // backward substitution
    val len = rhs.nRows * rhs.nCols
    rhs match {
      case MatRM(_,nCols,rVals) =>
        var rStart = len
        var i = nEntries(size)
        while( rStart > 0 )
        {
          val rEnd = rStart; rStart -= nCols
          var j = len
          while( j > rEnd )
          {
            i -= 1
            var k = rEnd
            while( k > rStart )
            {
              j -= 1; k -= 1
              rVals(k) -= rVals(j) * vals(i)
            }
          }
          i -= 1
          while( j > rStart )
          {
            j -= 1
            rVals(j) /= vals(i)
          }
        }
      case MatCM(_,_,rVals) =>
        val jEnd = nEntries(size)
        var i = len
        while( i > 0 )
        {
          val cEnd = i-1; var j = jEnd // <- cStart: start of currently solved column in rhs
          while( j > 0 )
          {
            i -= 1; j -= 1
            var k = cEnd
            while( k > i )
            {
              rVals(i) -= rVals(k) * vals(j)
              j -= 1; k -= 1
            }
            rVals(i) /= vals(j)
          }
        }
    }
  }

  override def choleskySumDerive( weights: LMat ): Unit =
  {
    // BACKTRACKING CHOLESKY_DECOMPOSTION
    //
    // for a matrix A the cholesky decomposition L with A = L *# L^T can be determined as follows
    //
    // L[i,j < i] = ( A[i,j] - L[i,0]*L[j,0] - L[i,1]*L[j,1] - ... - L[i,j-1]*L[j,j-1] ) / L[j,j]
    //
    // L[i,i] = √( A[i,i] - L[i,0]² - L[i,1]² - ... - L[i,j-1]² )
    //
    // =>
    //
    // dL[i,j < i]/dp = - dL[j,j]/dp * L[i,j] / L[j,j] + dA[i,j]/dp / L[j,j]
    //                  - ( dL[i,0]/dp*L[j,0] + dL[i,1]/dp*L[j,1] + ... + dL[i,j-1]/dp*L[j,j-1] ) / L[j,j]
    //                  - ( L[i,0]*dL[j,0]/dp + L[i,1]*dL[j,1]/dp + ... + L[i,j-1]*dL[j,j-1]/dp ) / L[j,j]
    //
    //
    if( weights.size != size )
      throw new IllegalArgumentException

    weights match {
      case LMatCM(_,rVals) =>
        val len = vals.length
        var i = len
        while( i > 0 )
        {
          val cEnd = i
          var j = len

          while( j > cEnd )
          {
            i -= 1
            var k = cEnd
            while( k > i )
            {
              j -= 1; k -= 1
              rVals(i) -= rVals(j) * vals(k)
              rVals(k) -= rVals(j) * vals(i)
            }
          }

          i -= 1; j -= 1
          while( j > i )
          {
            rVals(j) /= vals(i)
            rVals(i) -= rVals(j) * vals(j)
            j -= 1
          }
          rVals(i) /= 2*vals(i)
        }
    }
  }

  override def choleskyInvert(): Unit =
  {
    // STEP 1: calculate this^-1
    triInvert()
    // STEP 2: with the now inverted this, calculate (this^T ⋅ this), i.e. ◥⋅◣
    var i = 0
    var cSize = size
    while( cSize > 0 )
    {
      val cEnd = i + cSize
      var j = i
      while( i < cEnd )
      {
        var k = i
        var sum = 0.0
        while( k < cEnd )
        {
          sum += vals(j)*vals(k)
          j += 1; k += 1
        }
        vals(i) = sum
        i += 1
      }
      cSize -= 1
    }
  }

  override def symRowForeach[@specialized U]( row: Int )( consumer: (Double,Int) => U ) =
  {
    if( row >= size )
      throw new IndexOutOfBoundsException
    var i = row
    var j = 0
    while( j < row )
    {
      consumer(vals(i),j)
      j += 1; i += size-j
    }
    while( j < size )
    {
      consumer(vals(i),j)
      i += 1; j += 1
    }
  }

  /** Assuming <b>this</b> is the lower triangle of a symmetric matrix <b>A</b>
    * and <b>b</b> returns the lower triangle of a symmetric matrix <b>B</b>,
    * this methods computes the trace of the matrix product <b>A⋅B</b>.

    * @param b <b>b(i,j)</b> returns <b>B<sub>ij</sub></b> for <b>0 ≤ j ≤ i < a.size</b>
    * @return <b>tr(this⋅B)</b>
    */
  @inline
  def symMatProdTrace( b: (Int,Int) => Double ): Double =
  {
    var tr = 0.0
    var i = vals.length
    var c = size
    while( c > 0 )
    {
      c -= 1
      var r = size-1
      while( r > c )
      {
        i -= 1
        tr += 2 * vals(i) * b(r,c)
        r -= 1
      }
      i -= 1
      tr += vals(i) * b(r,c)
    }
    tr
  }

  override def symSquare: LMatCM =
  {
    val result = new LMatCM( size, Vec.zeros(nEntries{size}) )

    var cEnd = 0
    var cSize = size

    while( cSize > 0 ) // <- for each column (left -> right)
    {
      val cStart = cEnd
      var i = cStart
      var j = cStart
      cEnd += cSize; cSize -= 1

      while( j < cEnd )
      {
        var k = j

        while( k < cEnd ) {
                                     result.vals(j) += vals(i)*vals(k)
          if( j != cStart && k > j ) result.vals(k) += vals(i)*vals(j)
          if( i != cStart )          result.vals(i) += vals(k)*vals(j)
          i += 1; k += 1
        }

        j += 1
      }
    }
    result
  }

  @inline override def symSquareLike( product: ProductLike ): LMatCM =
  {
    val result = new LMatCM( size, Vec.zeros(nEntries{size}) )

    var cEnd = 0
    var cSize = size
    var col_jk = 0 // <- column(j) == column(k)

    while( cSize > 0 )
    {
      val cStart = cEnd; cEnd += cSize; cSize -= 1
      var i = cStart
      var j = cStart
      var rj_ci = col_jk // <- row(j) == column(i)

      while( j < cEnd )
      {
        var row_ik = rj_ci // <- row(i) == row(k)
        var k = j

        while( k < cEnd ) {
//           assert( vals(i) == this(row_ik, rj_ci ) )
//           assert( vals(j) == this(rj_ci,  col_jk) )
//           assert( vals(k) == this(row_ik, col_jk) )
                                     result.vals(j) += product( vals(i), vals(k),  rj_ci,  row_ik, col_jk)
          if( j != cStart && k > j ) result.vals(k) += product( vals(i), vals(j),  row_ik, rj_ci,  col_jk)
          if( i != cStart )          result.vals(i) += product( vals(k), vals(j),  row_ik, col_jk, rj_ci )

          row_ik += 1; i += 1; k += 1
        }

        rj_ci += 1; j += 1
      }

      col_jk += 1
    }

    result
  }
}
object LMatCM extends LMatFactory
{
  def tabulate( size: Int )( idx2val: (Int,Int) => Double ): LMatCM =
  {
    val result = zeros(size)
    result.modify( (_,i,j) => idx2val(i,j) )
    result
  }

  override def zeros( size: Int ): LMatCM =
  {
    if( size < 1        ) throw new IllegalArgumentException("size may not be less than 1.")
    if( size > MAX_SIZE ) throw new IllegalArgumentException("Matrix size exceeds limit.")
    new LMatCM( size, Vec.zeros(nEntries{size}) )
  }

  /** A column major way of creating a new lower triangular matrix using a curried function
    * to initialize the entry values.
    *
    * @param size
    * @param colRow2val The tabulation function column -> row -> entry_value.
    * @return
    */
  def tabulateCM( size: Int )( colRow2val: Int => Int => Double ): LMatCM =
  {
    val result = zeros(size)
    var i = result.vals.length
    var c = size
    while( c > 0 )
    {
      var r = size; c -= 1
      val row2val = colRow2val(c)
      while( r > c )
      {
        r -= 1; i -= 1
        result.vals(i) = row2val(r)
      }
    }
    result
  }

  /** Assuming <b>a</b> is the lower triangle of a symmetric matrix <b>A</b>
    * and <b>b</b> is the lower triangle of a symmetric matrix <b>B</b>, this
    * methods computes the trace of the matrix product <b>A⋅B</b>.
    *
    * @param a The lower triangle of a symmetric matrix <b>A</b>.
    * @param b The lower triangle of a symmetric matrix <b>B</b>
    * @return <b>tr(A⋅B)</b>
    */
  def symMatProdTrace( a: LMatCM, b: LMatCM ): Double =
  {
    if( a.size != b.size ) throw new IllegalArgumentException("size of a and b do not match.")
    var tr = 0.0
    var i = 0
    var cSize = a.size
    while( cSize > 0 )
    {
      cSize -= 1
      tr += a.vals(i) * b.vals(i)
      var j = cSize
      while( j > 0 )
      {
        i += 1; j -= 1
        tr += 2 * a.vals(i) * b.vals(i)
      }
      i += 1
    }
    tr
  }

  def unapply( lMat: LMatCM )
    = Some(lMat.size,lMat.vals)
}