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