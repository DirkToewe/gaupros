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

package gps.regression

import gps.kernel.{Kernel, Sym}
import gps.linalg.{LMat, Vec}
import gps.opt.{BFGS_B_Optimizer, ObjectiveWithHessian}
import gps.regression.gpr.{LikelihoodFunction, LogP_LOO, LogP_LOO_Shifted, LogP_Marginal, LogP_Marginal_Shifted}

import scala.Double.{PositiveInfinity => ∞}
import scala.collection.immutable.Map
import scala.math.{exp, log, Pi => π}

/** A <a href="https://en.wikipedia.org/wiki/Gaussian_process">Gaussian Process Regression</a> implementation.
  *
  * <h2>See also</h2>
  *   <dt>[1]
  *   <dd><a href="http://www.gaussianprocess.org/gpml/chapters/RW2.pdf">
  *     C. E. Rasmussen & C. K. I. Williams, Gaussian Processes for Machine Learning,
  *     Chapter 2, the MIT Press, 2006, ISBN 026218253X
  *   </a>
  *   <dt>[2]
  *   <dd><a href="http://www.ressources-actuarielles.net/EXT/ISFA/1226.nsf/0/e7dc33e4da12c5a9c12576d8002e442b/$FILE/Jones01.pdf">
  *     Donald R. Jones, A, A Taxonomy of Global Optimization Methods Based on Response Surfaces,
  *     Kluwer Academic Publishers, 2001, Journal of Global Optimization 21: 345–383
  *   </a>
  * <dl>
  *
  * <p> Created by Dirk Toewe on 18.06.17.
  *
  * @param x The training inputs, i.e.
  * @param y The (shifted) training outputs. Noisy if `var_noise > 0`
  * @param yShift The value by which y is shifted. The Gaussian Process is fit on (y - y_shift). y_shift is then
  *               again added to Gaussian Process model prediction.
  * @param kernel The covariance (or kernel) function.
  * @param K_chol The cholesky decomposition of the coviance (or Gram) matrix.
  * @tparam X The type of training inputs.
  */
class GPR[@specialized X](x: Array[X], y: Vec, yShift: Double, kernel: Kernel[X], K_chol: LMat, val params: Map[Symbol,Double] )
{
  assert( x.length == y.length )
  assert( ! yShift.isInfinite )
  assert( ! yShift.isNaN )
  /** The "weights" vector.
    */
  private val w = y.clone
  K_chol.choleskySolve(w)

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
    K_chol.triSolve(ks)
    kernel(x,x) - (ks ⋅ ks) max 0
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
    K_chol.triSolve(ks)
    ( yShift+mean, kernel(x,x) - (ks ⋅ ks) max 0 )
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
object GPR
{
  /**
    *
    * @param x The training inputs.
    * @param y The training outputs. Noisy if `var_noise > 0`
    * @param yShift The value by which y is shifted. The Gaussian Process is fit on (y - y_shift). y_shift is then
    *                again added to Gaussian Process model prediction.
    * @param kernel The covariance (or kernel) function.
    * @tparam X The type of training inputs.
    * @return
    */
  def apply[@specialized X]( x: Array[X], y: Vec, yShift: Double, kernel: Kernel[X] ): GPR[X] =
  {
    if( x.length != y.length )
      throw new IllegalArgumentException("x and y must have same length.")

    val xArr = x.clone
    val yArr = y.clone
    var i = yArr.length
    while( i > 0 )
    {
      i -= 1
      yArr(i) -= yShift
    }
    val K = LMat.tabulate(x.length){ (i,j) => kernel(xArr{i},xArr{j}, i,j) }
    K.choleskyDecomp()
    new GPR[X](xArr, yArr, yShift, kernel, K, Map.empty)
  }

  /**
    *
    * @param x
    * @param y
    *
    * @param kernel     The covariance function in dependence of its parameters (parameters => kernel).
    *
    * @param param_min  The lower bound of the individual kernel (hyper)parameters. Parameter i may no be less than param_min[i].
    * @param param_max  The upper bound of the individual kenrel (hyper)parameters. Parameter i may no be greater than param_max[i].
    * @param param_init The initial parameter values for the optimization.
    *
    * @tparam X The regression input type, i.e. the type of the explanatory variable(s).
    */
  def fit[@specialized X](
    x: Array[X],
    y: Vec,
    y_shift: Double,

    kernel: Kernel[X],

    likelihood_fn: (Array[X],Vec,Double,Kernel[X]) => LikelihoodFunction = logp_marginal(_: Array[X],_: Vec,_: Double,_: Kernel[X]),

    param_init: Map[Symbol,Double] = Map.empty.withDefault(_ => 1.0 ), // <- assume one to be a reasonable initial value
    param_min : Map[Symbol,Double] = Map.empty.withDefault(_ => 1e-7), // <- assume all parameters to be positive
    param_max : Map[Symbol,Double] = Map.empty.withDefault(_ => ∞   )
  ): GPR[X] = {
    if( x.length != y.length )
      throw new IllegalArgumentException("x and y must have same length.")

    val likelihood = likelihood_fn(x,y, y_shift, kernel)

    val p0   = Vec( likelihood.params.map(param_init): _* )
    val pMin = Vec( likelihood.params.map(param_min ): _* )
    val pMax = Vec( likelihood.params.map(param_max ): _* )

    val p = BFGS_B_Optimizer() maximize (likelihood, p0,pMin,pMax)
    val params = (likelihood.params zip p).toMap

    val cov = kernel subs (params.toSeq :_*)

    val xArr = x.clone
    val yArr = y.clone
    var i = yArr.length
    while( i >  0 ) {
           i -= 1
      yArr(i) -= y_shift
    }
    val K = LMat.tabulate(x.length){ (i,j) => cov(xArr{i},xArr{j}, i,j) }
    K.choleskyDecomp()

    new GPR(xArr,yArr, y_shift, cov, K, params)
  }

  def fit_shifted[@specialized X](
    x: Array[X],
    y: Vec,
    y_shift: Symbol,

    kernel: Kernel[X],

    likelihood_fn: (Array[X],Vec,Symbol,Kernel[X]) => LikelihoodFunction,

    param_init: Map[Symbol,Double] = Map.empty,
    param_min : Map[Symbol,Double] = Map.empty,
    param_max : Map[Symbol,Double] = Map.empty
  ): GPR[X] = {
    if( x.length != y.length )
      throw new IllegalArgumentException("x and y must have same length.")

    val likelihood = likelihood_fn(x,y, y_shift, kernel)

    val parMin = (sym: Symbol) => param_min .getOrElse( sym, if(sym==y_shift) -(∞) else 1e-7 )
    val parMax =                  param_max .getOrElse(_: Symbol,∞)
    val par0   = (sym: Symbol) => param_init.getOrElse(sym, {
      val result = if( sym != y_shift ) 1.0 else y.sum / y.length
      result max parMin(sym) min parMax(sym)
    })

    val pMin = Vec( likelihood.params.map{parMin}: _* )
    val pMax = Vec( likelihood.params.map{parMax}: _* )
    val p0   = Vec( likelihood.params.map{par0  }: _* )

    val p = BFGS_B_Optimizer() maximize (likelihood, p0,pMin,pMax)
    val params = (likelihood.params zip p).toMap

    val cov = kernel subs (params.toSeq :_*)

    val yShift = params(y_shift)
    val xArr = x.clone
    val yArr = y.clone
    var    i = yArr.length
    while( i >  0 ) {
           i -= 1
      yArr(i) -= yShift
    }
    val K = LMat.tabulate(x.length){ (i,j) => cov(xArr{i},xArr{j}, i,j) }
    K.choleskyDecomp()

    new GPR(xArr,yArr, yShift, cov, K, params)
  }

  def logp_marginal[@specialized X]( x: Array[X], y: Vec, yShift: Double, cov: Kernel[X] ): LikelihoodFunction with ObjectiveWithHessian = new LogP_Marginal        (x,y,yShift,cov)
  def logp_marginal[@specialized X]( x: Array[X], y: Vec, yShift: Symbol, cov: Kernel[X] ): LikelihoodFunction                           = new LogP_Marginal_Shifted(x,y,yShift,cov)
  def logp_loo     [@specialized X]( x: Array[X], y: Vec, yShift: Double, cov: Kernel[X] ): LikelihoodFunction with ObjectiveWithHessian = new LogP_LOO             (x,y,yShift,cov)
  def logp_loo     [@specialized X]( x: Array[X], y: Vec, yShift: Symbol, cov: Kernel[X] ): LikelihoodFunction                           = new LogP_LOO_Shifted     (x,y,yShift,cov)
}
