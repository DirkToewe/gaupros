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
import gps.opt.test_functions.RosenbrockFunction
import utest._

import scala.util.Random

/**
  * Created by Dirk Toewe on 26.08.17.
  */
object BFGS_B_Optimizer_test extends TestSuite
{
  override def tests = this{
    'rosenbrock{
      val rng = new Random(1337)
      val bfgs = BFGS_B_Optimizer()

      val objective = RosenbrockFunction

      for( run <- 0 until 32 )
      {
        val d = 2+run

        for( _ <- 0 until 512 )
        {
          val x0 = Vec.tabulate(d)( _ => rng.nextDouble*3 - 1.5 )

          val x = bfgs.minimize(objective,x0)

          if( d <= 3 )
            assert( allClose( Vec.ones(d), x ) )
          else
            assert( allClose( Vec.zeros(d), objective.gradient(x), 1e-4, 1e-4 ) )
        }
        println(f"[BFGS/rosenbrock] run$run%4d check!")
      }
    }
  }
}
