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

import gps.linalg.Vec

import scala.Double.{PositiveInfinity => ∞}

/** Super-trait for exact and inexact <a href="https://en.wikipedia.org/wiki/Line_search">line search</a> methods.
  *
  * Created by Dirk Toewe on 27.08.17.
  */
trait LineSearchMethod
{
  // TODO: make it work for maximization as well ???

  /** Performs a line search to minimize/reduce an object function a long the given direction. Some line search
    * methods ensure certain conditions like the strong Wolfe conditions.
    *
    * @param dir The search direction along which the objective is to be minimized/reduced.
    * @param x_0 The starting point for the search.
    * @param φ_0 The objective function value at the starting point.
    * @param g_0 The objective gradient at the starting point.
    * @param α_init The initial value for the line search. Should be 1 for Newton and Quasi-Newton methods.
    * @param α_max The largest alpha that is to be considered for the line search.
    * @return x-Value, function value and gradient of the new point found by this line search method.
    */
  @throws[LineSearchFailure]
  def apply( dir: Vec, x_0: Vec, φ_0: Double, g_0: Vec, α_init: Double = 1, α_max: Double = ∞ ): (Vec,Double,Vec)
}
