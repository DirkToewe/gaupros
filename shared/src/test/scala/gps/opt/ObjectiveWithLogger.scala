package gps.opt

import gps.util.PlotlyUtil

import scala.collection.mutable.ArrayBuffer

/** Wrapper around objective function that keeps track of the number of time it is being called
  * and for which x-values it is being called.
  *
  *
  * Created by Dirk Toewe on 26.08.17.
  */
class ObjectiveWithLogger private[opt]( objective: ObjectiveFunction, print: Boolean ) extends ObjectiveFunction
{
  protected var _nCallsFunc = 0L // <\
  protected var _nCallsGrad = 0L // <== LongAdder would be more efficient but may work in ScalaJS
  protected var _nCallsHess = 0L // </

  def nCallsFunc = _nCallsFunc
  def nCallsGrad = _nCallsGrad
  def nCallsHess = _nCallsHess

  val inputs = ArrayBuffer.empty[Vec]

  protected[opt] def log( x: Vec ) =
  {
    if(print)
      println(x.toSeq + " -> " + objective(x) )
    inputs += x.clone
  }

  override def apply( x: Vec ) =
  {
    synchronized{
      log(x)
      _nCallsFunc += 1
    }
    objective(x)
  }

  /** Assuming the wrapped objective function has a two-dimensional input,
    * this function creates a 3D-HTML-plot of the call history and shows
    * it in the browser.
    *
    * @param marginAbs The absolute margin by which the plot is extended
    *                  around the bounding box of the call history input.
    * @param marginRel The relative margin by which the plot is extended
    *                  around the bounding box of the call history input.
    *                  More specifically this margin is relative to the
    *                  x and y size of said bounding box.
    */
  def plot( marginAbs: Double = 0, marginRel: Double = 0.01 ): Unit =
  {
    assert{ inputs.forall(_.length == 2) }

    val x = inputs.map{_(0)}
    val y = inputs.map{_(1)}
    val z = inputs map objective

    var xMin = x.min
    var xMax = x.max
    val marginX = marginAbs + (xMax-xMin)*marginRel
    xMin -= marginX
    xMax += marginX

    var yMin = y.min
    var yMax = y.max
    val marginY = marginAbs + (yMax-yMin)*marginRel
    yMin -= marginY
    yMax += marginY

    val xRange = xMin to xMax by (xMax-xMin)/1024
    val yRange = yMin to yMax by (yMax-yMin)/1024

    val markerSize = 3

    val surf = yRange.par map { y => xRange map ( x => objective{ Vec(x,y) } ) }

    PlotlyUtil.plot(
      layout = "{ title: 'Optimization Process', showlegend: true }",
      s"""{
      |            type: 'surface',
      |            colorscale: 'Viridis',
      |            name: 'surface',
      |            x: [ ${ xRange mkString ", " } ],
      |            y: [ ${ yRange mkString ", " } ],
      |            z: [
      |              ${ surf map {_ mkString ("[", ", ", "]")} mkString ",\n              " }
      |            ]
      |          }""".stripMargin,
      s"""{
      |            type: 'scatter3d',
      |            name: 'iterations',
      |            mode: 'markers',
      |            marker: {
      |              size: $markerSize,
      |              color: [ ${ 0 until inputs.length mkString ", " } ]
      |            },
      |            x: [ ${ x mkString ", " } ],
      |            y: [ ${ y mkString ", " } ],
      |            z: [ ${ z mkString ", " } ],
      |            text: [ ${ 0 until inputs.length map ("'iteration "+_+"'") mkString ", " } ]
      |          }""".stripMargin
    )
  }
}
private[opt] class ObjectiveLogGrad protected[opt](
  objective: ObjectiveWithGradient,
  print: Boolean=false
) extends ObjectiveWithLogger(objective,print) with ObjectiveWithGradient
{
  override def apply( x: Vec )
    = super[ObjectiveWithLogger].apply(x)

  override def gradient( x: Vec ) =
  {
    synchronized{
      log(x)
      _nCallsGrad += 1
    }
    objective.gradient(x)
  }

  override def fval_grad(x: Vec) =
  {
    synchronized{
      log(x)
      _nCallsFunc += 1
      _nCallsGrad += 1
    }
    objective.fval_grad(x)
  }
}
private[opt] class ObjectiveLogHess(
  objective: ObjectiveWithHessian,
  print: Boolean=false
) extends ObjectiveLogGrad(objective,print) with ObjectiveWithHessian
{
  override def apply    ( x: Vec ) = super[ObjectiveLogGrad].apply(x)
  override def gradient ( x: Vec ) = super[ObjectiveLogGrad].gradient(x)
  override def fval_grad( x: Vec ) = super[ObjectiveLogGrad].fval_grad(x)

  override def hessian(x: Vec) =
  {
    synchronized{
      log(x)
      _nCallsHess += 1
    }
    objective.hessian(x)
  }

  override def fval_grad_hess(x: Vec) =
  {
    synchronized{
      log(x)
      _nCallsFunc += 1
      _nCallsGrad += 1
      _nCallsHess += 1
    }
    objective.fval_grad_hess(x)
  }
}
object ObjectiveWithLogger
{
  def apply( objective: ObjectiveFunction ): ObjectiveWithLogger
    = apply(objective,false)

  def apply( objective: ObjectiveWithGradient ): ObjectiveWithLogger with ObjectiveWithGradient
    = apply(objective,false)

  def apply( objective: ObjectiveWithHessian ): ObjectiveWithLogger with ObjectiveWithHessian
    = apply(objective,false)

  def apply( objective: ObjectiveFunction, print: Boolean ): ObjectiveWithLogger
    = objective match {
        case obj: ObjectiveWithGradient => apply(obj,print)
        case _ => new ObjectiveWithLogger(objective,print)
      }

  def apply( objective: ObjectiveWithGradient, print: Boolean ): ObjectiveWithLogger with ObjectiveWithGradient
    = objective match {
        case obj: ObjectiveWithHessian => apply(obj,print)
        case _ => new ObjectiveLogGrad(objective,print)
      }

  def apply( objective: ObjectiveWithHessian, print: Boolean ): ObjectiveWithLogger with ObjectiveWithHessian
    = new ObjectiveLogHess(objective,print)
}