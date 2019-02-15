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
import utest._

import scala.util.Random

/** Super class for testing of ... uh well test functions.
  *
  * Created by Dirk Toewe on 26.08.17.
  */
class ObjectiveWithHessian_test(
  objective: ObjectiveWithHessian,
  x_min: Vec, x_max: Vec,
  tolRel: Double, tolAbs: Double
) extends TestSuite
{
  assert( x_min.length == x_max.length )

  val n = x_min.length;
  for( i <- 0 until n )
    assert{ x_min(i) < x_max(i) }

  override def tests = this{
    'gradient{
      val diff = gps.diff.gradient(objective)
      val rng = new Random(1337)

      for( _ <- 0 to 64*1024 )
      {
        val x = Vec.tabulate(n)(
          i => rng.nextDouble*{x_max(i) - x_min(i)} + x_min(i)
        )
        val d = diff(x)
        val g = objective.gradient(x)
//        println("{{")
//        println(d.toSeq)
//        println("=?=")
//        println(g.toSeq)
//        println("}}")
        assert( allClose(d,g, tolRel,tolAbs)  )
      }
    }

    'hessian{
      val hess = gps.diff.jacobian(objective.gradient)
      val rng = new Random(1337)

      for( _ <- 0 to 64*1024 )
      {
        val x = Vec.tabulate(n)(
          i => rng.nextDouble*{x_max(i) - x_min(i)} + x_min(i)
        )
        val D = hess(x)
        val H = objective.hessian(x)
//        println("{{")
//        println(D)
//        println("=?=")
//        println(H)
//        println("}}")
        H foreach {
          (Hij, i,j) => assert{ isClose( Hij, {D(i,j) + D(j,i)} / 2, tolRel,tolAbs) }
        }
      }
    }
  }
}
