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

import gps.kernel.Kernel
import gps.linalg._
import gps.opt.ObjectiveWithHessian

import scala.collection.{IndexedSeq => ISeq}
import scala.math.{log, Pi => π}

/** Marginal log-likelihood function.
  *
  * Created by Dirk Toewe on 10.08.17.
  */
class LogP_Marginal[X]( x: Array[X], y: Vec, yShift: Double, kernel: Kernel[X] ) extends LikelihoodFunctionAbstract[X](x,y,kernel,kernel.params)
                                                                                 with    ObjectiveWithHessian
{
  private def fval( K: LMat, β: Vec ): Double =
  {
    var logp = -0.5 * { (β ⋅ β) + nSamples*log(2*π) }
    K.diagForeach( logp -= log(_) )
    logp
  }


  private def grad( K: LMat, β: Vec, kGrad: Int => Kernel[X] ): Vec =
  {
    val dK = K.copy
    dK.triSolveSumDerive( β, β map (-_) )
    dK.diagModify{ (dKij,i) => dKij - 1 / K(i,i) }
    K.choleskySumDerive(dK)

    Vec.tabulate(params.length)(
      kGrad andThen {
        dkdp => dK.mapReduce{ (w,i,j) => w * dkdp(x{i},x{j}, i,j) }(_ + _)
      }
    )
  }


  private def hess( K: LMat, α: Vec, kGrad: Int => Kernel[X], kHess: (Int,Int) => Kernel[X] ): LMat =
  {
    LMatCM.tabulateCM(params.length){
      dq =>
        val dKdq = {
          val dkdq = kGrad(dq)
          Mat.tabulateSym(nSamples){ (i,j) => dkdq(x{i},x{j}, i,j) }
        }
        K.choleskySolve(dKdq)
        dp =>
          val dKdp = {
            val dkdp = kGrad(dp)
            Mat.tabulateSym(nSamples){ (i,j) => dkdp(x{i},x{j}, i,j) }
          }

          val A = dKdp *# dKdq
          var result = -2 * Mat.inner(α,A,α)

          K.choleskySolve(A)
          result += A.trace

          val B = {
            val ddk_dp_dq = kHess(dp,dq)
            Mat.tabulateSym(nSamples){ (i,j) => ddk_dp_dq(x{i},x{j}, i,j) }
          }

          result += Mat.inner(α,B,α)
          K.choleskySolve(B)
          0.5 * { result - B.trace }
    }
  }


  override def apply( param: Vec ): Double =
  {
    val k = cov(param)
    val β = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.triSolve(β)

    fval(K,β)
  }


  override def gradient( param: Vec ): Vec =
  {
    val (k,kGrad) = cov_grad(param)
    val β = y - yShift

    val K = LMatRM.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.triSolve(β)

    grad(K,β,kGrad)
  }


  override def fval_grad( param: Vec ): (Double,Vec) =
  {
    val (k,kGrad) = cov_grad(param)
    val β = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.triSolve(β)

    (
      fval(K,β),
      grad(K,β,kGrad)
    )
  }


  override def hessian( param: Vec ) =
  {
    // see: Moore CJ, Chua AJK, Berry CPL, Gair JR.
    // Fast methods for training Gaussian processes on large datasets.
    // R.Soc. opensci. 3: 160125.
    // http://dx.doi.org/10.1098/rsos.160125
    // 2016, page 4
    val (k,kGrad,kHess) = cov_grad_hess(param)
    val α = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.choleskySolve(α)

    hess(K,α,kGrad,kHess)
  }


  override def fval_grad_hess( param: Vec ) =
  {
    val (k,kGrad,kHess) = cov_grad_hess(param)
    val α = y - yShift
    val β = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.choleskySolve(α)
    K.triSolve(β)

    (
      fval(K,β),
      grad(K,β,kGrad),
      hess(K,α,kGrad,kHess)
    )
  }



// DO NOT DELETE!!! LEAVE AS EXPLANATORY PSEUDOCODE EXAMPLE
//
//  override def gradient_v2( param: Vec ) =
//  {
//    val (k,kGrad) = cov_grad(param)
//    val β = y - yShift
//
//    val K = LMatRM.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
//    K.choleskyDecomp()
//    K.triSolve(β)
//
//    // LogP = -0.5 * sum[i]( β[i]² + 2 * log(K[i,i]) ) + nSamples*log(2*π) }
//    //
//    // => dLogP/dp = - sum[i]( β[i] * dβ[i]/dp + dK[i,i]/dp / K[i,i] )
//    //
//    val dK = LMat.zeros(nSamples)
//    val dβ = β map (-_)
//
//    // BACKTRACKING TRI_SOLVE
//    //
//    // β[i] = (y[i] - K[i,0]*β[0] - K[i,1]*β[1] - ... - K[i,i-1]*β[i-1]) / K[i,i]
//    //
//    // dβ[i]/dp = - dK[i,i]/dp * β[i] / K[i,i]
//    //            - ( dK[i,0]/dp * β[0] + dK[i,1]/dp * β[1] + ... + dK[i,i-1]/dp * β[i-1] ) / K[i,i]
//    //            - ( K[i,0] * dβ[0]/dp + K[i,1] * dβ[1]/dp + ... + K[i,i-1] * dβ[i-1]/dp ) / K[i,i]
//    //
//    for( i <- nSamples-1 to 0 by -1 )
//    {
//      dβ(i) /= K(i,i)
//      for( j <- i   to 0 by -1 ) dK(i,j) = - dβ(i) * β(j)
//      for( j <- i-1 to 0 by -1 ) dβ(j)  -=   dβ(i) * K(i,j)
//    }
//
//    // - dLog(K[i,i])/dp = - dK[i,i]/dp / K[i,i]
//    //
//    for( i <- nSamples-1 to 0 by -1 )
//      dK(i,i) -= 1 / K(i,i)
//
//    // BACKTRACKING CHOLESKY_DECOMPOSTION
//    //
//    // for a matrix A the cholesky decomposition L with A = L *# L^T can be determined as follows
//    //
//    // L[i,j < i] = ( A[i,j] - L[i,0]*L[j,0] - L[i,1]*L[j,1] - ... - L[i,j-1]*L[j,j-1] ) / L[j,j]
//    //
//    // L[i,i] = √( A[i,i] - L[i,0]² - L[i,1]² - ... - L[i,j-1]² )
//    //
//    // =>
//    //
//    // dL[i,j < i]/dp = - dL[j,j]/dp * L[i,j] / L[j,j] + dA[i,j]/dp / L[j,j]
//    //                  - ( dL[i,0]/dp*L[j,0] + dL[i,1]/dp*L[j,1] + ... + dL[i,j-1]/dp*L[j,j-1] ) / L[j,j]
//    //                  - ( L[i,0]*dL[j,0]/dp + L[i,1]*dL[j,1]/dp + ... + L[i,j-1]*dL[j,j-1]/dp ) / L[j,j]
//    //
//    for( i <- nSamples-1 to 0 by -1 )
//    {
//      dK(i,i) /= 2
//
//      for( j <- i to 0 by -1 )
//      {
//        dK(i,j) /= K(j,j)
//
//        if( j < i )
//          dK(j,j) -= dK(i,j) * K(i,j)
//
//        for( k <- j-1 to 0 by -1 )
//        {
//          dK(i,k) -= dK(i,j) * K(j,k)
//          dK(j,k) -= dK(i,j) * K(i,k)
//        }
//      }
//    }
//
//    Array.tabulate[Double](param.length)(
//      kGrad andThen {
//        dkdp => dK.mapReduce{ (w,i,j) => w * dkdp(x{i},x{j}, i,j) }(_ + _)
//      }
//    )
//  }



// DO NOT DELETE!!! LEAVE FOR HISTORIC REASONS AND COMPARISON TO NEWER METHODS
//
//  override def gradient_v1( param: Vec ) =
//  {
//    val (k,kGrad) = cov_grad(param)
//    val α = y - yShift
//
//    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
//    K.choleskyDecomp()
//    K.choleskySolve(α)
//    K.choleskyInvert()
//    K.modify{ (K_ij,i,j) => α(i)*α(j) - K_ij }
//
//    Array.tabulate[Double](param.length)(
//      kGrad andThen ( dkdp => 0.5 * K.symMatProdTrace{ (i,j) => dkdp(x{i},x{j}, i,j) } )
//    )
//  }
}
