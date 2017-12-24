package gps.opt

import gps.linalg.Vec

/** An optimizer allowing box constraints.
  *
  * Created by Dirk Toewe on 29.08.17.
  */
trait Optimizer_B[-F <: ObjectiveFunction] extends Optimizer[F]
{
  override def maximize( objective: F, x_init: Vec ) = maximize(objective, x_init, UNDEFINED_ARRAY, UNDEFINED_ARRAY)
  override def minimize( objective: F, x_init: Vec ) = minimize(objective, x_init, UNDEFINED_ARRAY, UNDEFINED_ARRAY)

  def maximize( objective: F, x_init: Vec=UNDEFINED_ARRAY, x_min: Vec=UNDEFINED_ARRAY, x_max: Vec=UNDEFINED_ARRAY ): Vec
  def minimize( objective: F, x_init: Vec=UNDEFINED_ARRAY, x_min: Vec=UNDEFINED_ARRAY, x_max: Vec=UNDEFINED_ARRAY ): Vec
}
object Optimizer_B {}