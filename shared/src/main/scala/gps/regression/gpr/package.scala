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

package gps.regression

import gps.linalg.Vec
import gps.util.PlotlyUtil

import scala.math.{sqrt => √}

/**
  * Created by Dirk Toewe on 29.08.17.
  */
package object gpr
{
  private[gpr] implicit val SYMBOL_ORDER: Ordering[Symbol] = Ordering by (_.name)

  def plot1d( gpr: GPR[Double], x: Vec, y: Vec, marginAbs: Double=0.0, marginRel: Double=0.1 ) =
  {
    assert{ x.length == y.length }

    val (xMin,xMax) = (x.min, x.max)
    val margin = marginAbs + (xMax-xMin)*marginRel
    val xPlot = xMin-margin  to  xMax+margin  by  (xMax-xMin) / 1024
    val Seq(yMean,yLo,yHi) = {
       for(
         (mean,variance) <- xPlot.map(gpr.mean_var);
         σ = 2 * √(variance)
       )
         yield Seq( mean, mean - 2*σ, mean + 2*σ )
    }.transpose

    PlotlyUtil.plot(
      layout="{ title: 'Gaussian Process Regression Example' }",
      s"""
         |{
         |  type: 'scatter',
         |
         |  name: 'Mean - 2σ',
         |
         |  fill: 'none',
         |  mode: 'lines',
         |  line: { width: 0 },
         |
         |  legendgroup: 'sigma_band',
         |  x: [ ${xPlot mkString ", "} ],
         |  y: [ ${yLo mkString ", "} ],
         |}
      """.stripMargin,
      s"""
         |{
         |  type: 'scatter',
         |
         |  name: 'Mean + 2σ',
         |  legendgroup: 'sigma_band',
         |
         |  fill: 'tonexty',
         |  mode: 'none',
         |
         |  x: [ ${xPlot mkString ", "} ],
         |  y: [ ${yHi mkString ", "} ],
         |}
      """.stripMargin,
      s"""
         |{
         |  type: 'scatter',
         |  name: 'Mean',
         |  mode: 'lines',
         |  x: [ ${xPlot mkString ", "} ],
         |  y: [ ${yMean mkString ", "} ],
         |}
      """.stripMargin,
      s"""
         |{
         |  type: 'scatter',
         |  name: 'Samples',
         |  mode: 'markers',
         |  x: [ ${ x mkString ", " } ],
         |  y: [ ${ y mkString ", " } ],
         |}
      """.stripMargin
    )
  }
}
