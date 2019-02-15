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
object RastriginFunction_experiments
{
  def main( args: Array[String] ) =
  {
    val xRange = -5.0 to 5 by 0.1
    val yRange = -5.0 to 5 by 0.1
    val surf = yRange map { y => xRange map ( x => RastriginFunction{ Vec(x,y) } ) }

    PlotlyUtil.plot(
      layout = """{
      |          title: 'Rastrigin Function',
      |          scene: {
      |            aspectratio: { x: 1, y: 1, z: 0.25 },
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