package gps.opt.test_functions

import gps.util.PlotlyUtil

/**
  * Created by Dirk Toewe on 26.08.17.
  */
object RosenbrockFunction_experiments
{
  def main( args: Array[String] ) =
  {
    val xRange = -1.25 to 1.25 by 0.1
    val yRange = -1.0  to 2    by 0.1
    val surf = yRange map { y => xRange map ( x => RosenbrockFunction{ Vec(x,y) } ) }

    PlotlyUtil.plot(
      layout = "{ title: 'Rosenbrock Function' }",
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
