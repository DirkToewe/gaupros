package gps.opt

import gps.linalg.LMat

import scala.math.{acos, cos, sin, sqrt => √}

/** The UnboxedObjective offers allows to use unconstrained optimizers for constrained optimization
  * problems by wrapping the objective function with its constraints in a particular way. The resulting
  * objective function (called unboxed objective function) is defined for all real numbers. The underlying,
  * wrapped function however is only evaluated at inputs contained by the constraint box. Every minimum
  * and maximum and maximum of the boxed function has a corresponding maximum on the wrapped function.
  * Every maximum or minimum of the box constrained function has at least one corresponding maximum/minimum
  * on the unboxed function, including optima on the box's boundary. Therefore the unboxed objective function
  * allows to find the optima of the boxed function.
  */
object UnboxedObjective
{
  def apply( objective: ObjectiveFunction, x_min: Vec, x_max: Vec ): UnboxedObjective
    = objective match {
        case obj: ObjectiveWithGradient => this(obj,x_min,x_max)
        case _ =>
          val box = BoxConstraints(x_min,x_max)
          new UnboxedObjective(objective, box.x_min, box.x_max)
      }

  def apply( objective: ObjectiveWithGradient, x_min: Vec, x_max: Vec ): UnboxedObjective with ObjectiveWithGradient
    = objective match {
        case obj: ObjectiveWithHessian => this(obj,x_min,x_max)
        case _ =>
          val box = BoxConstraints(x_min,x_max)
          new UnboxedObjectiveWithGradient(objective, box.x_min, box.x_max)
      }

  def apply( objective: ObjectiveWithHessian, x_min: Vec, x_max: Vec ): UnboxedObjective with ObjectiveWithHessian =
  {
    val box = BoxConstraints(x_min,x_max)
    new UnboxedObjectiveWithHessian(objective, box.x_min, box.x_max)
  }

  private def cube( x: Double ) = x*x*x
  private def  sqr( x: Double ) = x*x

  private[opt] trait Trafo
  {
    def apply  ( x: Double ): Double
    def inverse( y: Double ): Double
    def deriv1 ( x: Double ): Double
    def deriv2 ( x: Double ): Double
  }


  private[opt] object TrafoIdentity extends Trafo
  {
    override def apply  ( x: Double ) = x
    override def inverse( y: Double ) = y
    override def deriv1 ( x: Double ) = 1
    override def deriv2 ( x: Double ) = 0
  }

  private[opt] class TrafoMin( min: Double ) extends Trafo
  {
    override def apply  ( x: Double ) = min + x*x / (1+x.abs)
    override def inverse( y: Double ) =
    {
      assert( y >= min )
      val z = y - min
      0.5 * ( z + √{z*(z+4)} )
    }
    override def deriv1 ( x: Double ) = x*(x.abs + 2) / sqr(x.abs+1) // <- deriv[ x² / (1+|x|) ] = x*(|x|+2) / (1+|x|)²
    override def deriv2 ( x: Double ) = 2 / cube(x.abs + 1)
  }

  private[opt] class TrafoMax( max: Double ) extends Trafo
  {
    override def apply  ( x: Double ) = max - x*x / (1+x.abs)
    override def inverse( y: Double ) =
    {
      assert( y <= max )
      val z = max - y
      0.5 * ( z + √{z*(z+4)} )
    }
    override def deriv1 ( x: Double ) = -x*(x.abs + 2) / sqr(x.abs+1) // <- deriv[ x² / (1+|x|) ] = x*(|x|+2) / (1+|x|)²
    override def deriv2 ( x: Double ) = -2 / cube(x.abs + 1)
  }

  private[opt] class Trafo2Bounds( min: Double, max: Double ) extends Trafo
  {
    assert( min <= max )
    override def apply  ( x: Double ) = min + (max-min)/2 * {1-cos(x)}
    override def inverse( y: Double ) =
    {
      assert( y >= min )
      assert( y <= max )
      acos{ 1 - 2 * (y-min) / (max-min) }
    }
    override def deriv1 ( x: Double ) = (max-min)/2 * sin(x)
    override def deriv2 ( x: Double ) = (max-min)/2 * cos(x)
  }
}
/** A wrapper wrapper function f(y ∈ ℝ) around a boxed b(x ∈ [x_min,x_max]).
  * Each minimum and maximum of f(y') has a corresponding extremum on b( x(y') )
  * where `transform_back` is function x(y).
  *
  * This allows to solve box constrained problems using unconstrained optimizers.
  *
  * Created by Dirk Toewe on 29.08.17.
  */
class UnboxedObjective protected[opt](unboxed: ObjectiveFunction, val x_min: Vec, val x_max: Vec ) extends ObjectiveFunction
{
  import UnboxedObjective._

  assert( x_min.length == x_max.length )
  @inline protected[opt] def n = x_min.length

  protected[opt] val trafos = Array.tabulate(n){
    i =>
      val xMin = x_min(i)
      val xMax = x_max(i)

      if( xMin.isInfinite )

        if( xMax.isInfinite ) TrafoIdentity
        else                  new TrafoMax(xMax)

      else if( xMax.isInfinite ) new TrafoMin(xMin)
      else                       new Trafo2Bounds(xMin,xMax)
  }

  /** Transforms an input of this (unboxed) function into the equivalent input for the wrapped (box constrained) function.
    *
    * This is used, for example, to transform the result of the optimization back into boxed coordinates.
    *
    * @param y
    * @return
    */
  def transform_back( y: Vec ): Vec
    = Vec.tabulate(n){ i => trafos(i){y(i)} }

  /** Transforms an input of the wrapped (box constrained) function into the equivalent input for this (unboxed) function.
    *
    * This is used, for example, to transform the stating point of the optimization into unboxed coordinates.
    *
    * @param x
    * @return
    */
  def transform( x: Vec ): Vec
    = Vec.tabulate(n){
        i =>
          assert{ x(i) >= x_min(i) }
          assert{ x(i) <= x_max(i) }
          trafos(i).inverse{x(i)}
      }

  override def apply( y: Vec ): Double
    = unboxed{ transform_back(y) }
}


private[opt] class UnboxedObjectiveWithGradient private[opt](
  unboxed: ObjectiveWithGradient, x_min: Vec, x_max: Vec
) extends UnboxedObjective(unboxed,x_min,x_max) with ObjectiveWithGradient
{
  override def apply( y: Vec ) = super[UnboxedObjective].apply(y)

  override def gradient( y: Vec ) =
  {
    val g = unboxed.gradient{ transform_back(y) }
    Vec.tabulate(n){ i => g(i) * trafos(i).deriv1(y{i}) }
  }

  override def fval_grad( y: Vec ) =
  {
    val (f,g) = unboxed.fval_grad{ transform_back(y) }
    val G = Vec.tabulate(n){ i => g(i) * trafos(i).deriv1(y{i}) }
    (f,G)
  }
}

private[opt] class UnboxedObjectiveWithHessian private[opt](
  unboxed: ObjectiveWithHessian, x_min: Vec, x_max: Vec
) extends UnboxedObjectiveWithGradient(unboxed,x_min,x_max) with ObjectiveWithHessian
{
  override def apply    ( y: Vec ) = super[UnboxedObjectiveWithGradient].apply(y)
  override def gradient ( y: Vec ) = super[UnboxedObjectiveWithGradient].gradient(y)
  override def fval_grad( y: Vec ) = super[UnboxedObjectiveWithGradient].fval_grad(y)

  override def hessian( y: Vec ) = fval_grad_hess(y)._3

  override def fval_grad_hess( y: Vec ): (Double, Vec, LMat) =
  {
    val (f,g,h) = unboxed.fval_grad_hess( transform_back(y) )
    assert( n == g.length )
    assert( n == h.size   )

    val G = Vec.tabulate(n){ i => trafos(i).deriv1(y{i}) }

    val H = h.copy
    H modify {
      (Hij, i,j) =>
        var result = Hij * G(i) * G(j)
        if( i == j )
          result += trafos(i).deriv2{y(i)} * g(i)
        result
    }

    var i = n
    while( i > 0 )
    {
      i -= 1
      G(i) *= g(i)
    }

    (f,G,H)
  }
}