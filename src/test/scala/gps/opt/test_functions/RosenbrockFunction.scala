package gps.opt.test_functions

import gps.linalg.{LMat, Vec}
import gps.opt.ObjectiveWithHessian

/**
  * Created by Dirk Toewe on 26.08.17.
  */
object RosenbrockFunction extends ObjectiveWithHessian
{
  val b = 100

  private def sqr( x: Double ) = x*x

  override def apply( x: Vec ) =
  {
    assert( x.length > 1 )
    (1 until x.length).map{
      i => b * sqr( x(i) - sqr{x(i-1)} ) + sqr( 1 - x(i-1) )
    }.sum
  }

  override def gradient( x: Vec ) =
  {
    assert( x.length > 1 )
    val d = x.length
    Vec.tabulate(d){
      i =>
        var result = 0.0
        if( i > 0  ) result += 2*b * ( x(i)   - sqr{x(i-1)} )
        if( i < d-1) result -= 4*b * ( x(i+1) - sqr{x(i  )} ) * x(i)  +  2*{1 - x(i)}
        result
    }
  }

  override def hessian( x: Vec ) =
  {
    assert( x.length > 1 )
    val d = x.length
    LMat.tabulate(d){
      (i,j) =>
        var result = 0.0
        if( i > 0 ) {
          if (j == i  ) result += 2*b
          if (j == i-1) result -= 4*b * x(i-1)
        }
        if( i < d-1 && j == i ) result -= 4*b* x(i+1) - 12*b * sqr{x(i)}  -  2
        result
    }
  }
}