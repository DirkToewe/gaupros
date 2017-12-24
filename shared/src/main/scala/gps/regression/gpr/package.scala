package gps.regression

import gps.util.PlotlyUtil

import scala.math.{sqrt => √}

/**
  * Created by Dirk Toewe on 29.08.17.
  */
package object gpr
{
  def plot1d( gpr: GPR[Double], x: Vec, y: Vec, marginAbs: Double=0.0, marginRel: Double=0.01 ) =
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
