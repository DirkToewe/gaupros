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