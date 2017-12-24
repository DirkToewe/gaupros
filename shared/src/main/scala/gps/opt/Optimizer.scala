package gps.opt

/**
  * Created by Dirk Toewe on 29.08.17.
  */
trait Optimizer[-F <: ObjectiveFunction]
{
  def maximize( objective: F, x_init: Vec ): Vec
  def minimize( objective: F, x_init: Vec ): Vec
}