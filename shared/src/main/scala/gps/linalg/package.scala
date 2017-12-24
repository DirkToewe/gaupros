package gps

import gps.opt.NoConvergence

import scala.annotation.tailrec
import scala.language.implicitConversions
import scala.math.{pow, sqrt => √}

/** Package containing rudimentary linear algebra utility method.
  *
  * Created by Dirk Toewe on 18.06.17.
  */
package object linalg
{
  private[linalg] val UNDEFINED_ARRAY = Vec.zeros(0)

  private[linalg] def isDefined  ( vec: Vec ) = UNDEFINED_ARRAY.vals ne vec.vals
  private[linalg] def isUndefined( vec: Vec ) = UNDEFINED_ARRAY.vals eq vec.vals

  implicit def vec2seq( vec: Vec ): IndexedSeq[Double] = vec.seq
//  implicit def seq2vec( seq: Seq[Double] ) = new Vec(seq.toArray)

  /** Solves Ax = b using the conjugate gradient method with Jacobi preconditioning.
    * The preconditioning should fare well for diagonally dominant matrices.
    *
    * <h2>See Also:</h2>
    *   <a href="https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method">Wikipedia - Conjugate Gradient Method</a>
    *   <br>
    *   <a href="https://en.wikipedia.org/wiki/Preconditioning#Jacobi_.28or_diagonal.29_preconditioner">Wikipedia - Preconditioning</a>
    *
    * @param A      Function that returns the entries of the symmetric positive definite matrix A.
    * @param b      The right hand side vector of the linear equation system.
    * @param x0     The start value for x for the iteration.
    * @param tolRel
    * @param tolAbs
    * @return
    * @throws NoConvergence If maxIter iterations is exceeded without convergence according to tolRel and tolAbs.
    */
  @throws[NoConvergence]
  def cg_jac(
    A: (Int,Int) => Double, b: Vec, x0: Vec=UNDEFINED_ARRAY,
    tolRel: Double=1e-4, tolAbs: Double=1e-4, maxIter: Long=Long.MaxValue
  ): Vec = {
    @inline def n = b.length

    val rem = Vec.zeros(n) // <- remainder of Kahan summation
    /** Writes the matrix-vector product A⋅u into result.
      */
    def A_*#( u: Vec, result: Vec ): Unit =
    {
      // reset the remainder
      var i = n
      while( i > 0 ) {
        i -= 1
        result(i) = 0
        rem(i) = 0
      }

      i = n
      while( i > 0 )
      {
        var j = i; i -= 1
        while( j > 0 )
        {
          j -= 1
          val Aij = A(i,j)

          var a = Aij*u(j) - rem(i)
          var b = result(i) + a
          rem(i) = {b - result(i)} - a
          result(i) = b
          // utilize symmetry to reduce costly calls of A
          if( j != i ) {
            a = Aij*u(i) - rem(j)
            b = result(j) + a
            rem(j) = {b - result(j)} - a
            result(j) = b
          }
        }
      }
    }

    // the preconditioner is a diagonal matrix
    val M = Vec.tabulate(n){ i => A(i,i) }

    val x = if( isUndefined(x0) ) Vec.zeros(n) else x0.clone
    val Ap = Vec.zeros(n)
    A_*#(x,/*->*/Ap)
    val r = b - Ap
    val z = r / M
    val p = z.clone

    val iMax = maxIter min b.length
    var iter = iMax
    while( iter > 0 )
    {
      iter -= 1

      A_*#(p,/*->*/Ap)

      val rz = r⋅z
      val α = rz / (p⋅Ap)

      var converged = true

      var i = n
      while( i > 0 )
      {
        i -= 1
        x(i) += α *  p(i)
        r(i) -= α * Ap(i)
        z(i)  = r(i) / M(i)
        converged &&= r(i).abs <= tolAbs + tolRel * b(i).abs
      }

      if(converged) return x

      val β = (r⋅z) / rz
      i = n
      while( i > 0 )
      {
        i -= 1
        p(i) = p(i) * β  +  z(i)
      }
    }
    throw NoConvergence(x)
  }

  def allClose( a: Vec, b: Vec, tolRel: Double=1e-5, tolAbs: Double=1e-8, equalNaN: Boolean=false ): Boolean =
  {
    if( a.length != b.length )
      throw new IllegalArgumentException
    @tailrec def loop( i: Int ): Boolean
      = i < 0 || isClose(a(i), b(i), tolRel, tolAbs, equalNaN) && loop(i-1)
    loop(a.length-1)
  }

  @inline def isClose( a: Double, b: Double, tolRel: Double=1e-5, tolAbs: Double=1e-8, equalNaN: Boolean=false ): Boolean
    = if( equalNaN && a.isNaN )
        b.isNaN
      else
        (a-b).abs <= tolAbs + tolRel * (a.abs + b.abs)/2

  implicit class ScalarVectorOps( val scalar: Double ) extends AnyVal
  {
    def * ( vec: Vec )
    = vec*scalar
  }

  /** Returns the euclidean norm of a vector.
    *
    * @param u The vector.
    * @return  The euclidean norm of u.
    */
  def norm( u: Vec, p: Double = 2 ) =
  {
    assert( p >= 1 )

    var i = u.length - 1
    var max = u(i).abs
    while( i > 0 )
    {
      i -= 1
      max = max max u(i).abs
    }

    p match {
      case 2 =>
        val s = sum(u.length){ i => val ui = u(i) / max; ui*ui }
        max * √(s)
      case 1 => sum(u.length){ u(_).abs }
      case _ =>
        val s = sum(u.length){ i => pow( u(i).abs / max, p ) }
        max * pow(s, 1/p)
    }
  }

  def normPow( u: Vec, p: Double = 2 ) =
  {
    assert( 0 < p )
    var result = 0.0
    p match {
      case 2 => sum(u.length){ i => val ui = u(i); ui*ui  }
      case 1 => sum(u.length){ u(_).abs }
      case _ => sum(u.length){ i => pow(u(i).abs, p) }
    }
  }

  @inline def sum( n: Int )( summands: Int => Double ): Double =
  {
    // Kahan Summation. Usually great! Might be bad on underflow and overflow...
    var sum = 0.0
    var rem = 0.0
    var i = n
    while( i > 0 )
    {
      i -= 1
      val a = summands(i) - rem
      val b = sum + a
      rem = (b - sum) - a
      sum = b
    }
    sum
  }

  def distance( u: Vec, v: Vec, p: Double = 2 ): Double =
  {
    norm(u - v, p)
  }

  def distancePow( u: Vec, v: Vec, p: Double = 2 ): Double =
  {
    assert( u.length == v.length )
    normPow(u - v, p)
  }
}
