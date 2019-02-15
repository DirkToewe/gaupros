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

/** A ℝ<sup>n</sup> → ℝ function trait that is to be minimized/maximized during optimization.
  *
  * Created by Dirk Toewe on 17.07.17.
  */
trait ObjectiveFunction extends (Vec => Double)
{
  def apply( x: Vec ): Double

  def unary_- : ObjectiveFunction = x => -this(x)
}