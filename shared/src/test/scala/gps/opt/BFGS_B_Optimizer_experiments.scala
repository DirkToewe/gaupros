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

import gps.linalg.Vec
import gps.opt.test_functions.RosenbrockFunction

/**
  * Created by Dirk Toewe on 26.08.17.
  */
object BFGS_B_Optimizer_experiments
{
  private def plotOpt( objective: ObjectiveWithHessian, x0: Vec ) =
  {
    val obj = ObjectiveWithLogger(objective, print=true)
    val x = try {
      BFGS_B_Optimizer()
        .minimize(obj,x0)
    }
    finally {
      println()
      println(s"${obj.nCallsGrad} Function Evaluations")
      obj.plot()
    }
    println(s"Solution: ${x.toSeq}")
  }

  def main( args: Array[String] ): Unit =
  {
    plotOpt( RosenbrockFunction, Vec(-0.5, 1.49) )
  }
}