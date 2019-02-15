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

package gps.opt

import gps.linalg._
import gps.opt.test_functions.{RastriginFunction, RosenbrockFunction}
import utest._

import scala.Double.{PositiveInfinity => ∞}
import scala.util.Random

/**
  *
  * Created by Dirk Toewe on 29.08.17.
  */
object UnboxedObjectiveSharp_test extends TestSuite
{
  private def testTrafo(unboxed: UnboxedObjectiveSharp with ObjectiveWithHessian, boxed: ObjectiveWithHessian ): Unit =
  {
    val rng = new Random(1337)
    val n = unboxed.x_min.length

    for( run <- 0 until 1024 )
    {
      val y = Vec.tabulate(n){ _ => rng.nextDouble*4 - 1.5 }

      // TEST RANDOMLY
      def test() =
      {
        val x = unboxed.transform_back(y)

        for( i <- 0 until n )
        {
          assert( x(i) >= unboxed.x_min(i) )
          assert( x(i) <= unboxed.x_max(i) )
        }

        assert{ isClose( boxed(x), unboxed(y) ) }
        assert{ isClose( boxed(x), unboxed.fval_grad(y)._1 ) }
        assert{ isClose( boxed(x), unboxed.fval_grad_hess(y)._1 ) }

        assert{ allClose(
          unboxed.transform_back(y),
          unboxed.transform_back( unboxed.transform(x) )
        ) }
      }

      test()

      // TEST BORDERS
      for( _ <- 1 to 8 ) // <- constrain up to 4 borders
      {
        y( rng.nextInt(n) ) = if(rng.nextBoolean) 0 else 1
        test()
      }
    }
  }

  private def testGrad( unboxed: UnboxedObjectiveSharp with ObjectiveWithHessian ): Unit =
  {
    val rng = new Random(1337)
    val n = unboxed.x_min.length

    val diff = gps.diff.gradient(unboxed)

    for( run <- 0 until 1024 )
    {
      val y = Vec.tabulate(n){ _ => rng.nextDouble*4 - 1.5 }

      // TEST RANDOMLY
      def test() =
      {
        assert{ allClose( unboxed.gradient(y),          diff(y), 1e-2, 1e-2 ) }
        assert{ allClose( unboxed.fval_grad(y)._2,      diff(y), 1e-2, 1e-2 ) }
        assert{ allClose( unboxed.fval_grad_hess(y)._2, diff(y), 1e-2, 1e-2 ) }
      }

      test()

      // TEST BORDERS
      for( _ <- 1 to 8 ) // <- constrain up to 4 borders
      {
        y( rng.nextInt(n) ) = if(rng.nextBoolean) 0 else 1
        test()
      }
    }
  }


  private def testHess( unboxed: UnboxedObjectiveSharp with ObjectiveWithHessian ): Unit =
  {
    val rng = new Random(1337)
    val n = unboxed.x_min.length

    val diff = gps.diff.jacobian(unboxed.gradient)

    for( run <- 0 until 1024 )
    {
      val y = Vec.tabulate(n){ _ => rng.nextDouble*4 - 1.5 }

      // TEST RANDOMLY
      def test() =
      {
        val D = diff(y)
        D.foreach{
          (Dij,i,j) => assert( isClose(Dij, D(j,i), 1e-4, 1e-4 ) )
        }
        def tst( H: LMat )
        = H.foreach{
          (Hij,i,j) => assert{ isClose(Hij, {D(i,j) + D(j,i)}/2, 1e-2, 1e-2) }
        }
        tst( unboxed hessian y )
        tst( unboxed.fval_grad_hess(y)._3 )
      }

      test()

      // TEST BORDERS
      for( _ <- 1 to 8 ) // <- constrain up to 4 borders
      {
        y( rng.nextInt(n) ) = if(rng.nextBoolean) 0 else 1
        test()
      }
    }
  }

  private def test( objective: ObjectiveWithHessian, x_min: Vec, x_max: Vec ): Unit =
  {
    val rng = new Random(1337)

    assert( x_min.length == x_max.length )
    val n = x_min.length

    val xMin = Vec.zeros(n)
    val xMax = Vec.zeros(n)

    for( run <- 0 until 32 )
    {
      for( i <- 0 until n )
      {
        rng.nextInt(4) match {
          case 0 =>
            xMin(i) = - ∞
            xMax(i) = + ∞
          case 1 =>
            xMin(i) = rng.nextDouble*(x_max(i) - x_min(i)) + x_min(i)
            xMax(i) = + ∞
          case 2 =>
            xMin(i) = - ∞
            xMax(i) = rng.nextDouble*(x_max(i) - x_min(i)) + x_min(i)
          case 3 =>
            val dx = x_max(i) - x_min(i)
            val width = rng.nextDouble * dx
            xMin(i) = rng.nextDouble*( dx - width ) + x_min(i)
            xMax(i) = xMin(i) + width
        }
      }

      val marginMax = (xMax-xMin).min min (x_max-x_min).min

      val unboxed = UnboxedObjectiveSharp(objective, xMin,xMax, rng.nextDouble * marginMax/2)

      testTrafo(unboxed, objective)
      testGrad(unboxed)
      testHess(unboxed)

      println( f"[UnboxedObjective/${objective.getClass.getSimpleName}${x_min.length}D] run$run%4d check!")
    }
  }

  override def tests = this{
    'rosenbrock {
      test( RosenbrockFunction, Vec(-2,-2), Vec(+2,+2) )
      test( RosenbrockFunction, Vec(-2,-2,-2), Vec(+2,+2,+2) )
      test( RosenbrockFunction, Vec(-2,-2,-2,-2), Vec(+2,+2,+2,+2) )
    }
    'rastrigin {
      test( RastriginFunction, Vec(-5,-5), Vec(+5,+5) )
      test( RastriginFunction, Vec(-5,-5,-5), Vec(+5,+5,+5) )
      test( RastriginFunction, Vec(-5,-5,-5,-5), Vec(+5,+5,+5,+5) )
    }
  }
}