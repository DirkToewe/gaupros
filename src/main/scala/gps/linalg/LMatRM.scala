package gps.linalg

import math.{sqrt => √}
import LMat.{MAX_SIZE, nEntries}
import gps.linalg.LMat.ProductLike

/**
  * Created by dtub on 19.06.17.
  */
// FIXME: vals.length is used multiple times in a way that is incompatible with choleskyAddRow
final class LMatRM private( private var _size: Int, private var vals: Vec ) extends LMat
{
  override def size =_size

  override def apply( row: Int, col: Int ): Double =
  {
    if( row < 0 ) throw new IndexOutOfBoundsException(   "Negative row index.")
    if( col < 0 ) throw new IndexOutOfBoundsException("Negative column index.")
    if( row >= size ) throw new IndexOutOfBoundsException(   "Row index out of bounds.")
    if( col >= size ) throw new IndexOutOfBoundsException("Column index out of bounds.")
    if( col > row )
      return 0
    vals{col + nEntries(row)}
  }

  override def update( row: Int, col: Int, value: Double ): Unit =
  {
    if( col < 0    ) throw new IndexOutOfBoundsException("Negative column index.")
    if( col > row  ) throw new IndexOutOfBoundsException("Upper right lesser triangle is read-only.")
    if( row >= size) throw new IndexOutOfBoundsException("Row index out of bounds.")
    vals{col + nEntries(row)} = value
  }

  override def modify( row: Int, col: Int, modifier: Double => Double ): Unit =
  {
    if( col < 0    ) throw new IndexOutOfBoundsException("Negative column index.")
    if( col > row  ) throw new IndexOutOfBoundsException("Upper right lesser triangle is read-only.")
    if( row >= size) throw new IndexOutOfBoundsException("Row index out of bounds.")
    val i = col + nEntries(row)
    vals(i) = modifier(vals{i})
  }

  override def modify( modifier: Mapper[Double] ): Unit =
  {
    var i = 0; var r=0
    while( r < size )
    {
      var c = 0
      while( c <= r )
      {
        vals(i) = modifier(vals(i),r,c)
        c += 1; i += 1
      }
      r += 1
    }
  }

  override def foreach[@specialized U]( consumer: Mapper[U] ): Unit =
  {
    var i = 0; var r=0
    while( r < size )
    {
      var c = 0
      while( c <= r )
      {
        consumer(vals(i),r,c)
        c += 1; i += 1
      }
      r += 1
    }
  }

  @inline override def mapReduce[M <: R, @specialized R]( mapper: Mapper[M] )
                                                        ( reducer: (R, R) => R ): R =
  {
    if( 0 == size )
      throw new IllegalStateException
    var i = nEntries(size) - 1
    var r = size-1
    var c = size-1
    var result: R = mapper(vals(i), r,c)
    while( r >= 0 )
    {
      while( c > 0 )
      {
        c -= 1; i -= 1
        result = reducer( result, mapper(vals(i), r,c) )
      }
      c = r; r -= 1
    }
    result
  }

  override def copy    = new LMatRM(_size,  vals.clone)
  override def unary_- = new LMatRM(_size, -vals)

  override def det =
  {
    var det = 1.0
    var  i = nEntries(size)
    var di = size
    while( di > 0 )
    {
      det *= vals(i)
      i -= di; di -= 1
    }
    det
  }

  override def triSolve( rhs: Vec ): Unit =
  {
    if( rhs.length != size )
      throw new IllegalArgumentException
    // forward substitution
    var iRow = 0; var iCol = 0 // <- indices within this lower triangular matrix
    while( iCol < size )
    {
      var j = 0
      while( j < iCol )
      {
        rhs(iCol) -= rhs(j) * vals(iRow+j)
        j += 1
      }
      rhs(iCol) /= vals(iRow+iCol)
      iCol += 1; iRow += iCol
    }
  }

  override def triSolve( rhs: Mat ): Unit =
  {
    if( rhs.nRows != size )
      throw new IllegalArgumentException
    val len = rhs.nRows * rhs.nCols
    rhs match {
      case MatRM(_,nCols,rVals) =>
        var i = 0; var rEnd = 0
        while( rEnd < len )
        {
          val rStart = rEnd; rEnd += nCols // <- start & end of currently solved row in rhs
          var j = 0
          while( j < rStart )
          {
            var k = rStart
            while( k < rEnd )
            {
              rVals(k) -= rVals(j) * vals(i)
              j += 1; k += 1
            }
            i += 1
          }
          while( j < rEnd )
          {
            rVals(j) /= vals(i)
            j += 1
          }
          i += 1
        }
      case MatCM(_,_,rVals) =>
        val jEnd = nEntries(size)
        var i = 0
        while( i < len )
        {
          val cStart = i; var j = 0 // <- cStart: start of currently solved column in rhs
          while( j < jEnd )
          {
            var k = cStart
            while( k < i )
            {
              rVals(i) -= rVals(k) * vals(j)
              j += 1; k += 1
            }
            rVals(i) /= vals(j)
            i += 1; j += 1
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
    var i = nEntries(size); var c = size
    while( i > 0 )
    {
      i -= 1; c -= 1
      weights(c) /= vals(i)
      vals(i) = - weights(c) * solution(c)
      var j = c
      while( j > 0 )
      {
        i -= 1; j -= 1
        weights(j) -= weights(c) * vals(i)
        vals(i)  =  - weights(c) * solution(j)
      }
    }
  }

  override def triInvert(): Unit =
  {
    // The basic idea is to solve this * this^-1 = I
    // with an efficient forward substitution scheme
    var rEnd = 0
    var rSize = 0
    while( rSize < size )
    {
      val rStart = rEnd
      rEnd += rSize; rSize += 1
      var i = 0
      var j = rStart
      while( j < rEnd )
      {
        var k = rStart
        while( k < j )
        {
          vals(k) -= vals(j)*vals(i)
          i += 1
          k += 1
        }
        vals(j) *= -vals(i)
        i += 1
        j += 1
      }

      vals(rEnd) = 1 / vals(rEnd)
      while( i < rEnd )
      {
        vals(i) *= vals(rEnd)
        i += 1
      }

      rEnd += 1
    }
  }

  override def diag(): Vec  =
  {
    val result = Vec.zeros(size)
    var rSize = size
    var i = -1 + nEntries(size)
    while( i >= 0 )
    {
      result(rSize-1) = vals(i)
      i -= rSize; rSize -= 1
    }
    result
  }

  @inline override def diagForeach[@specialized U]( consumer: Double => U ): Unit =
  {
    var rSize = size
    var i = -1 + nEntries(size)
    while( i >= 0 )
    {
      consumer(vals{i})
      i -= rSize; rSize -= 1
    }
  }

  @inline override def diagForeach[@specialized U]( consumer: (Double,Int) => U ): Unit =
  {
    var rSize = size
    var i = -1 + nEntries(size)
    while( i >= 0 )
    {
      consumer(vals{i}, rSize-1)
      i -= rSize; rSize -= 1
    }
  }

  @inline override def diagModify( modifier: (Double,Int) => Double ): Unit =
  {
    var rSize = size
    var i = -1 + nEntries(size)
    while( i >= 0 )
    {
      vals(i) = modifier(vals{i}, rSize-1)
      i -= rSize; rSize -= 1
    }
  }

  override def choleskyDecomp(): Unit =
  {
    var iRow = 0; var iCol = 0
    while( iCol < size )
    {
      var jRow = 0; var jCol = 0
      while( jCol <= iCol )
      {
        var sum = vals(iRow + jCol)
        var k = 0
        while( k < jCol )
        {
          sum -= vals(iRow + k) * vals(jRow + k)
          k += 1
        }
        if( iCol != jCol ) sum /= vals(jRow + jCol)
        else if( sum > +0.0 ) sum = √(sum)
        else
          throw MatrixNotPositiveDefinite()
        vals(iRow + jCol) = sum
        jCol += 1; jRow += jCol
      }
      iCol += 1; iRow += iCol
    }
  }

  override def choleskySolve( rhs: Vec ) =
  {
    // forward substitution
    triSolve(rhs)
    // backward substitution
    val len = nEntries(size)
    var iRow = len; var iCol = size
    while( iCol > 0 )
    {
      iRow -= iCol; iCol -= 1
      var jRow = len - size
      var jCol = size-1
      while( jCol > iCol )
      {
        rhs(iCol) -= rhs(jCol) * vals(jRow+iCol)
        jRow -= jCol; jCol -= 1
      }
      rhs(iCol) /= vals(iRow+iCol)
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
        var i = nEntries(size)
        var rStart = len
        while( rStart > 0 )
        {
          val rEnd = rStart; rStart -= nCols
          i -= 1
          var j = rEnd
          while( j > rStart )
          {
            j -= 1
            rVals(j) /= vals(i)
          }
          while( j > 0 )
          {
            i -= 1
            var k = rEnd
            while( k > rStart )
            {
              j -= 1; k -= 1
              rVals(j) -= rVals(k) * vals(i)
            }
          }
        }
      case MatCM(_,_,rVals) =>
        val jEnd = nEntries(size)
        var i = len
        while( i > 0 )
        {
          val cStart = i-size
          var j = jEnd
          while( j > 0 )
          {
            i -= 1; j -= 1
            rVals(i) /= vals(j)
            var k = i
            while( k > cStart )
            {
              j -= 1; k -= 1
              rVals(k) -= rVals(i) * vals(j)
            }
          }
        }
    }
  }

  override def choleskySumDerive( weights: LMat ) =
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
    if( weights.size != size )
      throw new IllegalArgumentException

    weights match {
      case LMatRM(_,rVals) =>
        var i = nEntries(size)-1
        var rSize = size
        while( i >= 0 )
        {
          rVals(i) /= 2 * vals(i)
          rSize -= 1
          val rStart = i-rSize
          var j = i
          while( j > rStart )
          {
            j -= 1
            rVals(j) -= 2 * rVals(i) * vals(j)
          }
          i -= 1; j -= 1
          while( j >= 0 )
          {
            rVals(i) /= vals(j)
            rVals(j) -= rVals(i) * vals(i)
            var k = i
            while( k > rStart )
            {
              j -= 1; k -= 1
              rVals(j) -= rVals(i) * vals(k)
              rVals(k) -= rVals(i) * vals(j)
            }
            i -= 1; j -= 1
          }
        }
    }
  }

  override def choleskyInvert(): Unit =
  {
    // STEP 1: calculate this^-T
    triInvert()
    // STEP 2: with the now inverted this, calculate (thisᵀ ⋅ this), i.e. ◥⋅◣
    var rEnd = 0 // <- the end of the current row (exclusive)
    var rSize = 0 // <- size of the current row
    while( rSize < size ) // <-
    {
      val rStart = rEnd // <- the beginning of the current row (is the end of the last)
      rEnd += rSize; rSize += 1
      var i = 0
      var j = rStart
      while( j < rEnd )
      {
        var k = rStart
        while( k <= j )
        {
          vals(i) += vals(j)*vals(k)
          k += 1
          i += 1
        }
        j += 1
      }
      while( i <= rEnd )
      {
        vals(i) *= vals(rEnd)
        i += 1
      }
      rEnd += 1
    }
  }

  def choleskyRowAppend(row: Vec ) =
  {
    if( size+1 != row.length) throw new IllegalArgumentException("New row must be of length (size+1).")
    if( size   == MAX_SIZE  ) throw new IllegalArgumentException("Matrix size exceeds limit.")

    // TODO make this proper triangular matrix sizes (otherwise elements of the array is never used)
    val len = nEntries(size)
    if( len+size+1 > vals.length )
    {
      val Vals = Vec.zeros( nEntries{size*3/2 max size+1 min MAX_SIZE} ) // <- TODO: maybe there's better ways but appending is O(n²) anyways...
      Vec._copy(vals,0, Vals,0, len)
      vals = Vals
    }
    Vec._copy(row,0, vals,len, size+1)
    _size += 1

    val iRow = len; val iCol = size-1
    var jRow = 0;   var jCol = 0
    while( jCol <= iCol )
    {
      var sum = vals(iRow + jCol)
      var k = 0
      while( k < jCol )
      {
        sum -= vals(iRow + k) * vals(jRow + k)
        k += 1
      }
      if( iCol != jCol ) sum /= vals(jRow + jCol)
      else if( sum > +0.0 ) sum = √(sum)
      else
        throw new IllegalArgumentException("Matrix not positive definite.")
      vals(iRow + jCol) = sum
      jCol += 1; jRow += jCol
    }
  }

  def choleskyRowRemoveLast(): Unit =
  {
    _size -= 1
    if( nEntries(size) < vals.length/4 )
      vals = vals.slice(0, nEntries{size*3/2 max size})
  }

  override def symRowForeach[@specialized U]( row: Int )( consumer: (Double,Int) => U ) =
  {
    var i = nEntries(row)
    var c = 0
    while( c < row )
    {
      consumer(vals(i),c)
      i += 1
      c += 1
    }
    while( c < size )
    {
      consumer(vals(i),c)
      c += 1
      i += c
    }
  }

  override def symMatProdTrace( b: (Int,Int) => Double ): Double =
  {
    var tr = 0.0
    var i = nEntries(size)
    var r = size
    while( r > 0 )
    {
      i -= 1
      r -= 1
      var c = r
      tr += vals(i) * b(r,c)
      while( c > 0 )
      {
        i -= 1
        c -= 1
        tr += 2 * vals(i) * b(r,c)
      }
    }
    tr
  }

  override def symSquare: LMat =
  {
    val result = new LMatRM( size, Vec.zeros(nEntries{size}) )

    var cStart = nEntries(size)
    var cSize = size

    while( cSize > 0 ) // <- for each row (last to first)
    {
      val cEnd = cStart-1
      var i = cEnd
      var j = cEnd
      cStart -= cSize; cSize -= 1

      while( j >= cStart )
      {
        var k = j

        while( k >= cStart ) {
                                   result.vals(j) += vals(i)*vals(k)
          if( j != cEnd && k < j ) result.vals(k) += vals(i)*vals(j)
          if( i != cEnd )          result.vals(i) += vals(j)*vals(k)
          i -= 1; k -= 1
        }

        j -= 1
      }
    }
    result
  }

  @inline override def symSquareLike( product: ProductLike ): LMatRM =
  {
    val result = new LMatRM( size, Vec.zeros(nEntries{size}) )

    var cStart = nEntries(size)
    var cSize = size
    @inline def row_jk = cSize

    while( cSize > 0 ) // <- for each column (left -> right)
    {
      val cEnd = cStart-1
      var i = cEnd
      var j = cEnd
      cStart -= cSize; cSize -= 1
      var ri_cj = row_jk

      while( j >= cStart )
      {
        var k = j
        var col_ik = ri_cj

        while( k >= cStart ) {
//          assert( vals(i) == this(ri_cj,  col_ik) )
//          assert( vals(j) == this(row_jk, ri_cj ) )
//          assert( vals(k) == this(row_jk, col_ik) )
                                   result.vals(j) += product( vals(k), vals(i), row_jk, col_ik, ri_cj )
          if( j != cEnd && k < j ) result.vals(k) += product( vals(j), vals(i), row_jk, ri_cj,  col_ik)
          if( i != cEnd )          result.vals(i) += product( vals(j), vals(k), ri_cj,  row_jk, col_ik)
          i -= 1; k -= 1; col_ik -= 1
        }

        j -= 1; ri_cj -= 1
      }
    }
    result
  }
}
object LMatRM extends LMatFactory
{
  def tabulate( size: Int )( idx2val: (Int,Int) => Double ): LMatRM =
  {
    val result = zeros(size)
    result.modify( (_,i,j) => idx2val(i,j) )
    result
  }

  override def zeros( size: Int ): LMatRM =
  {
    if( size < 0        ) throw new IllegalArgumentException("size may not be negative.")
    if( size > MAX_SIZE ) throw new IllegalArgumentException("Matrix size exceeds limit.")
    new LMatRM( size, Vec.zeros(nEntries{size}) )
  }

  /** A row major way of creating a new lower triangular matrix using a curried function
    * to initialize the entry values.
    *
    * @param size
    * @param rowCol2val The tabulation function row -> column -> entry_value.
    * @return
    */
  def tabulateRM( size: Int )( rowCol2val: Int => Int => Double ): LMatRM =
  {
    val result = zeros(size)
    var i = result.vals.length
    var r = size
    while( r > 0 )
    {
      var c = r; r -= 1
      val col2val = rowCol2val(r)
      while( c > 0 )
      {
        c -= 1; i -= 1
        result.vals(i) = col2val(c)
      }
    }
    result
  }

  def unapply( lMat: LMatRM )
    = Some(lMat.size,lMat.vals)
}