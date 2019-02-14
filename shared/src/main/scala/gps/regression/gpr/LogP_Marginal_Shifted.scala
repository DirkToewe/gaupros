package gps.regression.gpr

import gps.kernel.Kernel
import gps.linalg._
import gps.opt.ObjectiveWithHessian
import gps.util.Seqs.binSearch

import scala.collection.{IndexedSeq => ISeq}
import scala.math.{log, Pi => π}

class LogP_Marginal_Shifted[X]( x: Array[X], y: Vec, yShift: Symbol, kernel: Kernel[X] )
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

  private def fval( K: LMat, β: Vec ): Double =
  {
    var logp = -0.5 * { (β ⋅ β) + nSamples*log(2*π) }
    K.diagForeach( logp -= log(_) )
    logp
  }


  private def grad( K: LMat, β: Vec, nParams: Int, kGrad: Int => Kernel[X] ): Vec =
  {
    val dK = K.copy
    dK.triSolveSumDerive( β, β map (-_) )
    dK.diagModify{ (dKij,i) => dKij - 1 / K(i,i) }
    K.choleskySumDerive(dK)

    val γ = Vec.ones(nSamples)
    K.triSolve(γ)

    val g = Vec.tabulate(nParams){
      case `skip` => 0.0
      case i =>
        val dkdp = kGrad(i)
        dK.mapReduce{ (w,i,j) => w * dkdp(x{i},x{j}, i,j) }(_ + _)
    }
    g(idx) += β ⋅ γ
    g
  }

  override def apply( param: Vec ): Double =
  {
    val yShift = param(idx)
    val k = cov(param)
    val β = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.triSolve(β)

    fval(K,β)
  }


  override def gradient( param: Vec ): Vec =
  {
    val yShift = param(idx)
    val (k,kGrad) = cov_grad(param)
    val β = y - yShift

    val K = LMatRM.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.triSolve(β)

    grad(K,β,param.length,kGrad)
  }


  override def fval_grad( param: Vec ): (Double,Vec) =
  {
    val yShift = param(idx)
    val (k,kGrad) = cov_grad(param)
    val β = y - yShift

    val K = LMat.tabulate(nSamples){ (i,j) => k(x{i},x{j}, i,j) }
    K.choleskyDecomp()
    K.triSolve(β)

    (
      fval(K,β),
      grad(K,β,param.length,kGrad)
    )
  }
}
