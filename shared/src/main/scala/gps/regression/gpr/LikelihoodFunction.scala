package gps.regression.gpr

import gps.kernel.Kernel
import gps.linalg.Vec
import gps.opt.ObjectiveWithHessian

import scala.collection.{IndexedSeq => ISeq}

/**
  * Created by Dirk Toewe on 20.08.17.
  */
abstract class LikelihoodFunction[@specialized -X]( x: Array[X], y: Vec, yShift: Double, kernel: Kernel[X] ) extends ObjectiveWithHessian
{
  if( x.length != y.length ) throw new IllegalArgumentException

  protected val nSamples = x.length

  def params = kernel.params

  protected val grad_cache:      ISeq[Kernel[X]]  = Array.tabulate(params.length){ kernel pDiff params(_) }
  protected val hess_cache: ISeq[ISeq[Kernel[X]]] = Array.tabulate(params.length){
    i =>
      val dkdp = kernel pDiff params(i)
      Array.tabulate(i+1){ dkdp pDiff params(_) }: IndexedSeq[Kernel[X]]
  }

  protected def cov( param: Vec ) =
  {
    if( params.length != param.length )
      throw new IllegalArgumentException(s"${params.length} parameters expected but ${param.length} given.")
    val sym2val = Map( params zip param map { case (k,v) => k -> (v: Kernel[X]) } :_* )
    kernel subs sym2val
  }

  protected def cov_grad( param: Vec ): (Kernel[X], Int => Kernel[X]) = {
    if( params.length != param.length )
      throw new IllegalArgumentException(s"${params.length} parameters expected but ${param.length} given.")
    val sym2val = Map( params zip param map { case (k,v) => k -> (v: Kernel[X]) } :_* )
    (
      kernel subs sym2val,
      grad_cache map (_ subs sym2val)
    )
  }

  protected def cov_grad_hess( param: Vec ): (Kernel[X], Int => Kernel[X], (Int,Int) => Kernel[X]) =
  {
    if( params.length != param.length )
      throw new IllegalArgumentException(s"${params.length} parameters expected but ${param.length} given.")
    val sym2val = Map( params zip param map { case (k,v) => k -> (v: Kernel[X]) } :_* )
    (
      kernel subs sym2val,
      grad_cache map (_ subs sym2val), // <- cache 1st order derivatives as they may be used multiple
      hess_cache(_)(_) subs sym2val    //    times in the 2nd derivatives of a likelihood function
    )
  }
}