package gps.opt

import gps.opt.test_functions.RosenbrockFunction

/**
  * Created by Dirk Toewe on 27.08.17.
  */
object GradientDescentOptimizer_experiments
{
  private def plotOpt( objective: ObjectiveWithHessian, x0: Vec ) =
  {
    val obj = ObjectiveWithLogger(objective, print=true)
    val x = try {
      GradientDescentOptimizer
//        .withExponentialDecayHybrid(1e-3,1e-3, 1e-2,1e-4, 2 * 1024)
        .withExponentialDecayAbs(4,1e-4, 1024)
        .minimize( obj, x0 )
//        .minimize( obj, x0, x_min = Vec(-2, -2), x_max = Vec(0, +2) )
    }
    finally {
      println()
      println(s"${obj.nCallsGrad} Function Evaluations")
      obj.plot()
    }
    println(s"Solution: ${x.toSeq}")
  }

  def main( args: Array[String] ) =
  {
    plotOpt( RosenbrockFunction, Vec(-0.5, 1.49) )
  }
}