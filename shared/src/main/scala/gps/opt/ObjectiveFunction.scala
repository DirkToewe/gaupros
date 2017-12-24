package gps.opt

/** A ℝ<sup>n</sup> → ℝ function trait that is to be minimized/maximized during optimization.
  *
  * Created by Dirk Toewe on 17.07.17.
  */
trait ObjectiveFunction extends (Vec => Double)
{
  def apply( x: Vec ): Double

  def unary_- : ObjectiveFunction = x => -this(x)
}