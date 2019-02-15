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

package gps.opt.test_functions

import gps.linalg.Vec
import gps.util.PlotlyUtil

/**
  * Created by Dirk Toewe on 26.08.17.
  */
object AckleysFunction_experiments
{
  def main( args: Array[String] ) =
  {
    println( AckleysFunction( Vec(0,0) ) )
    println( AckleysFunction.gradient( Vec(0,0) ).toSeq )
    println( AckleysFunction.hessian( Vec(1e-300,0) ) )

    val xRange = -16.0 to 16 by 0.1
    val yRange = -16.0 to 16 by 0.1
    val surf = yRange map { y => xRange map ( x => AckleysFunction{ Vec(x,y) } ) }

    PlotlyUtil.plot(
      layout = """{
      |          title: "Ackley's Function",
      |          scene: {
      |            aspectratio: { x: 1, y: 1, z: 0.5 },
      |            aspectmode: 'manual'
      |          }
      |        }""".stripMargin,
      s"""{
      |            type: 'surface',
      |            colorscale: 'Viridis',
      |            opacity: 1.0,
      |            x: [ ${ xRange mkString ", " } ],
      |            y: [ ${ yRange mkString ", " } ],
      |            z: [
      |              ${ surf map {_ mkString ("[", ", ", "]")} mkString ",\n              " }
      |            ]
      |          }""".stripMargin
    )
  }
}
