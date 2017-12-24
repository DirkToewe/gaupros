package gps.regression

import gps.kernel.Kernel
import gps.linalg
import gps.linalg.Vec

import scala.math.{exp, log, Pi => π}

/**
  *
  * Created by Dirk Toewe on 21.07.17.
  */
class LargeGPR[@specialized X]( x: Array[X], y: Vec, yShift: Double, kernel: Kernel[X] )
{
  assert( x.length == y.length )
  assert( ! yShift.isInfinite )
  assert( ! yShift.isNaN )

  private val K = (i: Int, j: Int) => kernel(x{i},x{j}, i,j)
  /** The "weights" vector.
    */
  private val w = linalg.cg_jac(K,y)


  def nSamples = x.length

  def apply( x: X ): Double = mean(x)

  /** Computes the mean for the specified input. <b>O(n)</b> operation.
    *
    * @param x The input value.
    * @return E[y|x]
    */
  def mean( x: X ): Double =
  {
    var i = nSamples; var result = yShift
    while( i > 0 )
    {
      i -= 1
      result += kernel(this.x(i),x, i,-1) * w(i)
    }
    result
  }

  /** Computes variance for the specified input. <b>O(n<sup>2</sup>)</b> operation.
    *
    * @param x The input value.
    * @return V[y|x]
    */
  def variance( x: X ): Double =
  {
    val ks = Vec.tabulate(nSamples){ i => kernel(this.x(i),x, i,-1) }
    kernel(x,x) - {linalg.cg_jac(K,ks) ⋅ ks} max 0
  }

  /** Computes mean and variance for the specified input. <b>O(n<sup>2</sup>)</b> operation.
    *
    * @param x The input value.
    * @return (E[y|x], V[y|x])
    */
  def mean_var( x: X ): (Double,Double) =
  {
    val ks = Vec.tabulate(nSamples){ i => kernel(this.x(i),x, i,-1) }
    val mean = ks ⋅ w
    ( yShift+mean, kernel(x,x) - {linalg.cg_jac(K,ks) ⋅ ks} max 0 )
  }

  /** Returns the probability of y given x, i.e. p(y|x).
    *
    * @param x
    * @param y
    * @return
    */
  def p( x: X, y: Double ): Double =
    exp( logp(x,y) )

  /** Returns the logarithm of the probability of y given x, i.e. log(p(y|x)).
    *
    * @param x
    * @param y
    * @return
    */
  def logp( x: X, y: Double ): Double =
  {
    val (mean,σ2) = mean_var(x)
    val dy = y - mean
    val `2σ²` = 2*σ2
    - dy*dy / `2σ²` - log(π*`2σ²`) / 2
  }
}