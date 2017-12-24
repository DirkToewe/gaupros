package gps.opt

import gps.linalg._
import gps.opt.test_functions.{AckleysFunction, RastriginFunction, RosenbrockFunction}
import utest._

import scala.math.{atan2, cos, sin, Pi => π, toRadians => deg2rad}
import scala.util.Random

/**
  * Created by Dirk Toewe on 25.08.17.
  */
object LineSearchStrongWolfe_test extends TestSuite
{
  private def test( objective: ObjectiveWithGradient ) =
  {
    val lineSearch = LineSearchStrongWolfe(objective, c1=0.4, c2=0.8)
    // be more lenient during checking since numeric issues might make trouble
    val c1 = 0.3
    val c2 = 0.9

    val rng = new Random(1337)

    for( run <- 0 to 128 )
    {
      for( _ <- 0 to 1024 )
      {
        val x0 = Vec(
          rng.nextDouble*10 - 5,
          rng.nextDouble*10 - 5
        )
        val (f0,g0) = objective.fval_grad(x0)

        val dirAngle = atan2( g0(1), g0(0) ) + deg2rad( rng.nextDouble*80 - 40 )
        val dir = Vec(
          cos(dirAngle),
          sin(dirAngle)
        ) * -(rng.nextDouble*1.5 + 0.5)

        if( (g0⋅dir) < 1e-8 )
        {
          val (x,f,g) = lineSearch(dir, x0, f0, g0)

          val α = {(x - x0) ⋅ dir} / (dir⋅dir)

          val  φ_0= f0
          val dφ_0= g0⋅dir
          val  φ  = f
          val dφ  = g⋅dir

          assert( φ - φ_0 <=  c1 * dφ_0 * α )
          assert(dφ.abs   <= -c2 * dφ_0     )
        }
      }
      println( f"[LineSearchStrongWolfe/test${objective.getClass.getSimpleName}] run${run}%4d check!")
    }
  }

  override def tests = this{
    'testAckleys   { test(   AckleysFunction) }
    'testRastrigin { test( RastriginFunction) }
    'testRosenbrock{ test(RosenbrockFunction) }
  }
}