package gps.regression.gpr

import gps.kernel._
import gps.linalg._
import gps.opt.ObjectiveWithHessian

import scala.math.{log, Pi => π}

/** Leave-one-out log-likelihood function.
  *
  * Created by Dirk Toewe on 10.08.17.
  */
class LogP_LOO[X]( x: Array[X], y: Vec, yShift: Double, kernel: Kernel[X] ) extends LikelihoodFunctionAbstract[X](x,y,kernel,kernel.params)
                                                                             with    ObjectiveWithHessian
{
  private def fval( K_inv: LMat, α: Vec ) = {
    var logp = -nSamples * log(2*π)
    K_inv diagForeach { (K_ii,i) => logp += log(K_ii) - α(i)*α(i) / K_ii }
    logp / 2
  }

  private def grad( K_inv: LMat, K_diag: Vec, α: Vec, kGrad: Int => Kernel[X] ) =
  {
    // the gradient can be expressed as weighted sum of the dK/dp
    val dK = K_inv symSquareLike {
      (K_ik,K_kj, i,k,j) =>
        val result = ( α(k) * {α(i)*K_kj + α(j)*K_ik}  -  K_ik*K_kj * { 1 + α(k)*α(k) / K_diag(k) } ) / K_diag(k)
        if( i != j ) result else result/2
    }

    Vec.tabulate(params.length)(
      kGrad andThen {
        dkdp => dK.mapReduce{ (w,i,j) => w * dkdp(x{i},x{j}, i,j) }(_ + _)
      }
    )
  }

  private def hess( K: LMat, K_diag: Vec, α: Vec, kGrad: Int => Kernel[X], kHess: (Int,Int) => Kernel[X] ): LMat =
  {
    LMatRM.tabulateRM(params.length){
      dp =>

        val Z = {
          val dkdp = kGrad(dp)
          Mat.tabulateSym(nSamples){ (i,j) => dkdp(x{i},x{j}, i,j) }
        }
        K.choleskySolve(Z)

        val Zα = Z *# α
        val ZK = Z.copy; K.choleskySolve(ZK.T)

        dq =>

          val dK = {
            val dkdq = kGrad(dq)
            Mat.tabulateSym(nSamples){ (i,j) => -dkdq(x{i},x{j}, i,j) }
          }
          K.choleskySolve(dK)

          val dZ = {
            val ddk_dp_dq = kHess(dp,dq)
            Mat.tabulateSym(nSamples){ (i,j) => ddk_dp_dq(x{i},x{j}, i,j) }
          }
          K.choleskySolve(dZ)
          dZ += (dK *# Z)

          val dα = dK *# α

          K.choleskySolve(dK.T)

          val dZα = dZ *# α  +  Z *# dα
          val dZK = dZ.copy; K.choleskySolve(dZK.T); dZK += (Z *# dK)

          // { α(i)*Zα(i) - 0.5 * { 1 + α(i)*α(i) / K_diag(i) } * ZK(i,i) } / K_diag(i)
          var result = 0.0
          var i = nSamples
          while( i > 0 ) {
            i -= 1
            var sum = 0.0
            sum += dα(i)*Zα(i) + α(i)*dZα(i) - 0.5 * { 1 + α(i)*α(i) / K_diag(i) } * dZK(i,i)
            sum -= {     Zα(i) * α(i)        - 0.5 * { 1 + α(i)*α(i) / K_diag(i) } *  ZK(i,i) } / K_diag(i) * dK(i,i)
            sum -= {            dα(i)        - 0.5 *  dK(i,i) * α(i) / K_diag(i) } *  ZK(i,i)   / K_diag(i) * α(i)
            result += sum / K_diag(i)
          }
          result
    }
  }

  override def apply( param: Vec ): Double =
  {
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
    val (k,kGrad) = cov_grad(param)
    val α = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.choleskySolve(α)
    K.choleskyInvert()

    grad(K,K.diag, α, kGrad)
  }

  override def fval_grad( param: Vec ): (Double,Vec) =
  {
    val (k,kGrad) = cov_grad(param)
    val α = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.choleskySolve(α)
    K.choleskyInvert()

    (
      fval(K,α),
      grad(K,K.diag, α, kGrad)
    )
  }

  override def hessian( param: Vec ): LMat =
  {
    val (k,kGrad,kHess) = cov_grad_hess(param)
    val α = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.choleskySolve(α)

    val K_diag = { // <- diag(K^-1)
      val K_inv = K.copy
      K_inv.choleskyInvert()
      K_inv.diag
    }

    hess(K,K_diag, α, kGrad,kHess)
  }

  override def fval_grad_hess( param: Vec ): (Double, Vec, LMat) =
  {
    val (k,kGrad,kHess) = cov_grad_hess(param)
    val α = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.choleskySolve(α)

    val K_inv = K.copy
    K_inv.choleskyInvert()

    val K_diag = K_inv.diag

    (
      fval(K_inv,α),
      grad(K_inv,K_diag, α, kGrad),
      hess(K    ,K_diag, α, kGrad,kHess)
    )
  }

//  DO NOT DELETE! Gradient implementation according to Carl Rasmussen:
//
//  override def gradient( param: Vec ) = // <- TODO: use method above instead after testing
//  {
//    val (k,kGrad) = cov_grad(param)
//    val α = y - yShift
//
//    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
//    K.choleskyDecomp()
//    K.choleskySolve(α)
//
//    val _K = { // <- diag(K^-1)
//      val K_inv = K.copy
//      K_inv.choleskyInvert()
//      K_inv.diag
//    }
//
//    Array.tabulate[Double](param.length)(
//      kGrad andThen {
//        dkdp =>
//          val Z = Mat.tabulateSym(nSamples){ (i,j) => dkdp(x{i},x{j}, i,j) }
//          K.choleskySolve(Z)
//
//          val  Zα = Z *# α
//          val  ZK = Z.copy
//          K.choleskySolve(ZK.T)
//
//          var result = 0.0
//          for( i <- 0 until nSamples )
//            result += { α(i)*Zα(i) - 0.5 * { 1 + α(i)*α(i) / _K(i) } * ZK(i,i) } / _K(i)
//          result
//      }
//    )
//  }
}
