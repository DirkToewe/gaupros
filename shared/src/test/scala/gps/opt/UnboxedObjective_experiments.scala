package gps.opt

import gps.linalg.Vec
import gps.opt.test_functions.RosenbrockFunction
import gps.util.PlotlyUtil

import scala.math.{Pi => π}

/**
  * Created by Dirk Toewe on 28.08.17.
  */
object UnboxedObjective_experiments
{
  def main( args: Array[String] ): Unit =
  {
    val (xMin,xMax) = (-1.5, 0.5) // (-1.5,+0.75)
    val (yMin,yMax) = (-3.0, 2.0) // (-1.0, 2.0 )

    val objective = RosenbrockFunction

    val box = BoxConstraints(
      Vec(xMin,yMin),
      Vec(xMax,yMax)
    )

    val boxed = UnboxedObjective(objective, box.x_min, box.x_max)

    {
      val xRange = (xMin max -3) to (xMax min 3) by 0.01
      val yRange = (yMin max -3) to (yMax min 3) by 0.01
      val surf = yRange map { y => xRange map ( x => objective{ Vec(x,y) } ) }

      PlotlyUtil.plot(
        layout = "{ title: 'Box Constrained Function', showlegend: true }",
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

    {
      val xRange = -π/2 to +π*3/2 by 0.01
      val yRange = -π/2 to +π*3/2 by 0.01
      val surf = yRange map { y => xRange map ( x => boxed{ Vec(x,y) } ) }

      val xR = 0.0 to +π by 0.01
      val yR = 0.0 to +π by 0.01
      val boundaryX = xR               ++ yR.map(_ => π) ++ xR.reverse     ++ yR.map(_ => 0.0)
      val boundaryY = xR.map(_ => 0.0) ++ yR             ++ xR.map(_ => π) ++ yR.reverse
      val boundaryZ = for( (x,y) <- boundaryX zip boundaryY ) yield boxed( Vec(x,y) )

      PlotlyUtil.plot(
        layout = "{ title: 'Box Constrained Function', showlegend: true }",
        s"""{
        |            type: 'surface',
        |            colorscale: 'Viridis',
        |            opacity: 1.0,
        |            x: [ ${ xRange mkString ", " } ],
        |            y: [ ${ yRange mkString ", " } ],
        |            z: [
        |              ${ surf map {_ mkString ("[", ", ", "]")} mkString ",\n              " }
        |            ]
        |          }""".stripMargin,
        s"""{
        |            type: 'scatter3d',
        |            colorscale: 'Viridis',
        |            opacity: 1.0,
        |            mode: 'lines',
        |            x: [ ${ boundaryX mkString ", " } ],
        |            y: [ ${ boundaryY mkString ", " } ],
        |            z: [ ${ boundaryZ mkString ", " } ]
        |          }""".stripMargin
      )
    }

    val logged = ObjectiveWithLogger(objective, print=true)
    val x = try {
      BFGS_B_Optimizer().minimize(
        UnboxedObjective(logged, box.x_min, box.x_max),
        x_init = boxed transform Vec(-0.2,1.9)
      )
    }
    finally {
      println()
      println(s"${logged.nCallsGrad} Function Evaluations")
      logged.plot( marginAbs=0, marginRel=0 )
    }
    println(s"Solution: ${ boxed.transform_back(x).toSeq}")
  }
}