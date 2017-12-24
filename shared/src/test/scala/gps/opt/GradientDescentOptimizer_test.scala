package gps.opt

import gps.linalg._
import gps.opt.test_functions.RosenbrockFunction
import utest._

import scala.util.Random

/**
  * Created by Dirk Toewe on 27.08.17.
  */
object GradientDescentOptimizer_test extends TestSuite
{
  override def tests = this{
    'rosenbrock{
      val rng = new Random(1337)
      val gdOpt = GradientDescentOptimizer.withExponentialDecayAbs(2, 1e-6, 8*4096)

      val objective = RosenbrockFunction

      for( run <- 0 until 16 )
      {
        val d = 2+run

        for( _ <- 0 until 512 )
        {
          val x0 = Vec.tabulate(d)( _ => rng.nextDouble*3 - 1.5 )

          val x = gdOpt.minimize(objective,x0)

          if( d <= 3 )
            assert( allClose( Vec.ones(d), x ) )
          else
            assert( allClose( Vec.zeros(d), objective.gradient(x), 1e-2, 1e-2 ) )
        }
        println(f"[GD_OPT/rosenbrock] run$run%4d check!")
      }
    }
  }
}
