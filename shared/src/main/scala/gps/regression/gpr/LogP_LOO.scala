package gps.regression.gpr

import gps.kernel._
import gps.linalg._

import scala.math.{log, Pi => π}

/**
  * Created by Dirk Toewe on 10.08.17.
  */
class LogP_LOO[-X]( x: Array[X], y: Vec, yShift: Double, kernel: Kernel[X] ) extends LikelihoodFunction[X](x,y,yShift,kernel)
{
  override def apply( param: Vec ) =
  {
    val k  = cov(param)
    val α = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.choleskySolve(α)
    K.choleskyInvert()

    var logp = -nSamples * log(2*π)
    K diagForeach { (K_ii,i) => logp += log(K_ii) - α(i)*α(i) / K_ii }
    logp / 2
  }


  override def gradient( param: Vec ) =
  {
    val (k,kGrad) = cov_grad(param)
    val α = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.choleskySolve(α)
    K.choleskyInvert()

    // the gradient can be expressed as weighted sum of the dK/dp
    val diag = K.diag
    val dK = K symSquareLike {
      (K_ik,K_kj, i,k,j) =>
        val result = ( α(k) * {α(i)*K_kj + α(j)*K_ik}  -  K_ik*K_kj * { 1 + α(k)*α(k) / diag(k) } ) / diag(k)
        if( i != j ) result else result/2
    }

    Vec.tabulate(param.length)(
      kGrad andThen {
        dkdp => dK.mapReduce{ (w,i,j) => w * dkdp(x{i},x{j}, i,j) }(_ + _)
      }
    )
  }


  override def hessian( param: Vec ) =
  {
    val (k,kGrad,kHess) = cov_grad_hess(param)
    val α = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.choleskySolve(α)

    val _K = { // <- diag(K^-1)
      val K_inv = K.copy
      K_inv.choleskyInvert()
      K_inv.diag
    }

    LMatRM.tabulateRM(param.length){
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

          // { α(i)*Zα(i) - 0.5 * { 1 + α(i)*α(i) / _K(i) } * ZK(i,i) } / _K(i)
          var result = 0.0
          var i = nSamples
          while( i > 0 ) {
            i -= 1
            var sum = 0.0
            sum += dα(i)*Zα(i) + α(i)*dZα(i) - 0.5 * { 1 + α(i)*α(i) / _K(i) } * dZK(i,i)
            sum -= {     Zα(i) * α(i)        - 0.5 * { 1 + α(i)*α(i) / _K(i) } *  ZK(i,i) } / _K(i) * dK(i,i)
            sum -= {            dα(i)        - 0.5 *  dK(i,i) * α(i) / _K(i) } *  ZK(i,i)   / _K(i) * α(i)
            result += sum / _K(i)
          }
          result
    }
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