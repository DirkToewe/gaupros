/* This file is part of GauProS.
 *
 * GauProS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GauProS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GauProS.  If not, see <https://www.gnu.org/licenses/>.
 */

package gps.opt

import gps.linalg._

import scala.Double.{PositiveInfinity => ∞}
import scala.annotation.tailrec
import scala.math.pow

trait GradientDescentOptimizer
{
  /** The method used to calculate the learning rate used to find the next point.
    *
    * @param step The current iteration's index.
    * @param x The current function input.
    * @param gradient The gradient of the function at the current input x.
    * @return The learning rate for the current iteration. If 0 is returned, the
    *         iteration stops.
    */
  def learningRate( step: Long, x: Vec, gradient: Vec ): Double

  def maximize(
    objective: ObjectiveWithGradient,
    x_init  : Vec = UNDEFINED_ARRAY,
    x_min   : Vec = UNDEFINED_ARRAY,
    x_max   : Vec = UNDEFINED_ARRAY
  ) = optimize(+1, objective, x_init, x_min, x_max)

  def minimize(
    objective: ObjectiveWithGradient,
    x_init  : Vec = UNDEFINED_ARRAY,
    x_min   : Vec = UNDEFINED_ARRAY,
    x_max   : Vec = UNDEFINED_ARRAY
  ) = optimize(-1, objective, x_init, x_min, x_max)

  private def optimize(
    dir: Int,
    objective: ObjectiveWithGradient,
    x_init  : Vec = UNDEFINED_ARRAY,
    x_min   : Vec = UNDEFINED_ARRAY,
    x_max   : Vec = UNDEFINED_ARRAY
  ): Vec =
  {
    assert( dir == -1 || dir == +1 )

    if( null == objective) throw new NullPointerException

    val (x,box) = BoxConstraints(x_init,x_min,x_max)

    @tailrec def iterate( step: Long ): Vec =
    {
      val g = objective.gradient(x)
      box.clipDir(x,g, sign=dir)

      val learnRate = dir * learningRate(step,x,g)

      if( 0 == learnRate )
        x
      else {
        if( ! g.forall{! _.isNaN} ) throw new IllegalArgumentException("gradient contains NaN.")
        if(      learnRate.isNaN  ) throw new IllegalArgumentException( "learning_rate is NaN.")

        var i = x.length
        while( i > 0 ) {
          i -= 1
          x(i) += g(i) * learnRate
        }
        box.clip(x)

        iterate(step+1)
      }
    }
    iterate(0)
  }
}
object GradientDescentOptimizer
{
  def withExponentialDecay( learnRate_start: Double, learnRate_end: Double, nSteps: Long ): GradientDescentOptimizer =
  {
    assert( 1 < nSteps )
    assert(learnRate_start > 0)
    assert(learnRate_end   > 0)
    val lr = learnRate_end / learnRate_start
    val div = nSteps - 1.0
    (step, _, _) =>
      if (step >= nSteps)
        0
      else {
        val result = learnRate_start * pow(lr, step / div)
        assert( ! result.isNaN )
        result
      }
  }

  def withExponentialDecayAbs( dx_start: Double, dx_end: Double, nSteps: Long, tol: Double = 1e-7 ): GradientDescentOptimizer =
  {
    assert( 1 < nSteps )
    val dx = dx_end / dx_start
    val div = nSteps - 1.0
    (step, _, g) =>
      if (step >= nSteps)
        0
      else {
        val norm_g = norm(g)
        if( norm_g <= tol )
          0
        else
          dx_start * pow(dx, step / div) / norm_g
      }
  }

  def withExponentialDecayHybrid(
    learnRate_start: Double, learnRate_end: Double,
    dx_start: Double, dx_end: Double,
    nSteps: Long, tol: Double = 1e-7
  ): GradientDescentOptimizer =
  {
    assert( 1 < nSteps )
    val lr = learnRate_end / learnRate_start
    val dx = dx_end / dx_start
    val div = nSteps - 1.0
    (step, _, g) =>
      if (step >= nSteps)
        0
      else {
        val norm_g = norm(g) // <- a "normalized" norm would be better like mean squared error
        if( norm_g <= tol )
          0
        else
          learnRate_start * pow(lr, step / div) + dx_start * pow(dx, step / div) / norm_g
      }
  }

  /** Uses learning rate l(s) = l(0) / (1 + a*s)^p.
    *
    * @param learnRate_start
    * @param learnRate_after
    * @param decay_steps
    * @param max_steps
    * @param p
    * @return
    */
  def withNegativePowerDecay(
    learnRate_start: Double, learnRate_after: Double,
    decay_steps: Long, max_steps: Long,
    p: Double = 1
  ): GradientDescentOptimizer =
  {
    assert( 0 < learnRate_start)
    assert( 0 < learnRate_after)
    assert( 0 < p )
    val a = { pow(learnRate_start/learnRate_after, 1/p ) - 1} / {decay_steps-1} // <- FIXME
    (step,_,_) =>
      if (step >= max_steps)
        0
      else
        learnRate_start / pow(1 + a * step, p)
  }
}