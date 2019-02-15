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

package gps.regression.gpr

import gps.kernel._
import gps.linalg._
import gps.util.Seqs.binSearch

import scala.math.{log, Pi => π}

/** Leave-one-out log-likelihood function.
  *
  * Created by Dirk Toewe on 10.08.17.
  */
class LogP_LOO_Shifted[X](x: Array[X], y: Vec, yShift: Symbol, kernel: Kernel[X] )
  extends LikelihoodFunctionAbstract[X](
    x,y, kernel, params = binSearch(kernel.params, yShift)( Ordering by (_.name) ) match {
      case i if i >= 0 => kernel.params
      case i =>
        ( kernel.params.take(~i) :+ yShift ) ++ kernel.params.drop(~i)
    }
  )
{
  private val idx  = binSearch(       params, yShift) ensuring (_ >= 0)
  private val skip = binSearch(kernel.params, yShift) match {
    case i if i >= 0 => -1
    case i           => ~i
  }


  private def fval( K_inv: LMat, α: Vec ): Double =
  {
    var logp = -nSamples * log(2*π)
    K_inv diagForeach { (K_ii,i) => logp += log(K_ii) - α(i)*α(i) / K_ii }
    logp / 2
  }


  private def grad( K_inv: LMat, K_diag: Vec, α: Vec, δ: Vec, kGrad: Int => Kernel[X] ): Vec =
  {
    // the gradient can be expressed as weighted sum of the dK/dp
    val dK = K_inv symSquareLike {
      (K_ik,K_kj, i,k,j) =>
        val result = ( α(k) * {α(i)*K_kj + α(j)*K_ik}  -  K_ik*K_kj * { 1 + α(k)*α(k) / K_diag(k) } ) / K_diag(k)
        if( i != j ) result else result/2
    }

    val g = Vec.tabulate(params.length){
      case `skip` => 0.0
      case i =>
        val dkdp = kGrad(i)
        dK.mapReduce{ (w,i,j) => w * dkdp(x{i},x{j}, i,j) }(_ + _)
    }
    K_diag foreach { (K_ii, i) => g(idx) += δ(i)*α(i) / K_ii }
    g
  }


  override def apply( param: Vec ): Double =
  {
    val yShift = param(idx)
    val k  = cov(param)
    val α = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.choleskySolve(α)
    K.choleskyInvert()

    fval(K,α)
  }


  override def gradient( param: Vec ): Vec =
  {
    val yShift = param(idx)
    val (k,kGrad) = cov_grad(param)
    val α = y - yShift
    val δ = Vec.ones(nSamples)

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.choleskySolve(α)
    K.choleskySolve(δ)
    K.choleskyInvert()

    grad(K,K.diag, α,δ, kGrad)
  }

  override def fval_grad( param: Vec ): (Double, Vec) =
  {
    val yShift = param(idx)
    val (k,kGrad) = cov_grad(param)
    val α = y - yShift
    val δ = Vec.ones(nSamples)

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.choleskySolve(α)
    K.choleskySolve(δ)
    K.choleskyInvert()

    (
      fval(K,α),
      grad(K,K.diag, α,δ, kGrad)
    )
  }
}
