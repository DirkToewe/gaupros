package gps.regression

import gps.kernel._
import gps.linalg.Vec
import gps.opt.{BFGS_B_Optimizer, ObjectiveWithLogger}
import gps.util.PlotlyUtil

import scala.math.sin
import scala.util.Random

/**
  * Created by Dirk Toewe on 27.08.17.
  */
object GPR_experiments
{
  def main( args: Array[String] ): Unit =
  {
    plotFit1D()
  }

  def plotFit1D() =
  {
    val rng = new Random(1337)
    val nSamples = 1000
    val x = Array.tabulate(nSamples)( -6 + 12.0 * _ / (nSamples-1) )
    val y = Vec( x map { sin(_) - 0.1 + 0.2 * rng.nextDouble }: _* )
    val yShift = y.sum / y.length

    val cov = Noise('var_noise) + 'var_func * Exp( - AbsDelta.pow(1.9) * 'theta )

    val gpr = GPR.fit(x,y, yShift, cov)
    println(gpr.params)

//    gps.regression.gpr.plot1d(gpr,Vec(x: _*),y)
  }

  def plotTraining() =
  {
    val rng = new Random(1337)
    val nSamples = 128
    val x = (0 until nSamples).map( -6 + 12.0 * _ / (nSamples-1) ).toArray
    val y = Vec( x map { sin(_) - 0.1 + 0.2 * rng.nextDouble }: _* )
    val yShift = 0 // y.sum / y.length

    var cov = Noise( 'var_noise.pow(2) ) + 'var_func * Exp( - AbsDelta.pow(1.9) * 'theta )
    println(cov.params)
    cov = cov subs 'var_func -> 1

//    val likelihood = GPR.logp_marginal(x,y, yShift, cov)
    val likelihood = GPR.logp_loo(x,y, yShift, cov)
    val logged = ObjectiveWithLogger(likelihood, print=true)
    //    val result = BFGS_Optimizer().maximize(likelihood, Array(0.75, 1, 0.1) )
    val result = BFGS_B_Optimizer( tolRel=1e-5, tolAbs=1e-6 )
      .maximize( logged, x_init = Vec(0.75, 0.1), x_min = Vec(1e-2,1e-2) )

    val kernel = cov.subs( cov.params zip result :_* )
    val gpr = GPR(x,y, yShift, kernel)
    val Y = x map ( gpr(_) )

    println()
    println( s"${logged.nCallsGrad} function calls" )
    println(s"Solution: ${result.toSeq}")
    println(s"Likelihood: ${likelihood(result)}")
    println(s"Likelihood Gradient: ${likelihood.gradient(result).toSeq}")

    PlotlyUtil.plot(
      layout="{ title: 'Gaussian Process Regression Example' }",
      s"""
         |{
         |  type: 'scatter',
         |  name: 'training_data',
         |  mode: 'markers',
         |  x: [ ${ x mkString ", " } ],
         |  y: [ ${ y mkString ", " } ],
         |}
      """.stripMargin,
      s"""
         |{
         |  type: 'scatter',
         |  name: 'gpr',
         |  mode: 'lines',
         |  x: [ ${ x mkString ", " } ],
         |  y: [ ${ Y mkString ", " } ],
         |}
      """.stripMargin
    )

    logged.plot( marginAbs=0, marginRel=0 )
  }
}
