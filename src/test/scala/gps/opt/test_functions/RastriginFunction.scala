package gps.opt.test_functions

import gps.linalg.{LMat, Vec}
import gps.opt.ObjectiveWithHessian

import scala.math.{cos, sin, Pi => π}

/** The n-dimensional <a href="https://en.wikipedia.org/wiki/Rastrigin_function">Rastrigin function</a>.
  * May be called with arguments of any dimensions.
  *
  * Created by Dirk Toewe on 26.08.17.
  */
object RastriginFunction extends ObjectiveWithHessian
{
  val A = 10.0

  override def apply( x: Vec ) =
  {
    val n = x.length
    A*n + (0 until n).map( i => x(i)*x(i) - A*cos{2*π*x(i)} ).sum
  }

  override def gradient( x: Vec )
    = Vec.tabulate(x.length){
        i => 2*x(i) + 2*π*A * sin( 2*π * x{i} )
      }

  override def hessian(x: Vec)
    = LMat.tabulate(x.length){
        (i,j) => if( i != j ) 0 else 2 + 4*π*π*A * cos( 2*π * x{i} )
      }
}
