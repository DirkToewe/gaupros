package gps.opt.test_functions

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