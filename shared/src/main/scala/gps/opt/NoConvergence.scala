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

/** Exception thrown by iterative optimizers and solvers if the maximum number of iterations
  * is exceeded without achieving the demanded tolerance in the solution.
  *
  * @param x The iteration result which does not meet the tolerance requirements.
  *
  * Created by Dirk Toewe on 23.08.17.
  */
case class NoConvergence( x: Vec ) extends RuntimeException {}