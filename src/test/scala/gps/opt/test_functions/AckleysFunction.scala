package gps.opt.test_functions

import gps.linalg._
import gps.opt.ObjectiveWithHessian

import scala.math.{E, cos, exp, sin, Pi => π, sqrt => √}

/**
  * Created by Dirk Toewe on 26.08.17.
  */
object AckleysFunction extends ObjectiveWithHessian
{
  val a = 20
  val b = 0.2
  val c = 2*π

  override def apply( x: Vec ) =
  {
    val d = x.length
    val sumSqr = sum(d){ i =>  x(i)*x(i)  / d }
    val sumCos = sum(d){ i => cos(c*x{i}) / d }
    -a*exp( -b * √(sumSqr) ) - exp(sumCos) + a + E
  }

  override def gradient( x: Vec ) =
  {
    val d = x.length
    val sumCos =    sum(d){ i => cos(c*x{i}) / d }
    val sumSqrt= √{ sum(d)( i =>  x(i)*x(i)  / d ) }
    if( sumSqrt == 0 )
      Vec.zeros(d)
    else
      Vec.tabulate(d){
        i =>
          a*b/d * exp( -b * sumSqrt ) / sumSqrt * x(i)  +  c/d * exp(sumCos) * sin(c*x{i})
      }
  }

  override def hessian( x: Vec ) =
  {
    val d = x.length

    val sumSqr = sum(d)( i =>  x(i)*x(i)  / d )
    val sumCos = sum(d){ i => cos(c*x{i}) / d }
    val sumSqrt= √{sumSqr}

    LMat.tabulate(d){
      (i,j) =>
        var result = -a*b * exp( -b * sumSqrt ) / sumSqr * x(i) * x(j) * (b + 1/sumSqrt) +
                     -c*c * exp(sumCos) * sin(c*x{i}) * sin(c*x{j})
        result /= d
        if( i == j )
          result += a*b * exp( -b * sumSqrt ) / sumSqrt  +  c*c * exp(sumCos) * cos(c*x{i})
        result / d
    }
  }
}