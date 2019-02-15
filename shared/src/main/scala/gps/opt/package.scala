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

package gps

import gps.linalg.Vec

/** Rudimentary optimization methods and utilities.
  *
  * Created by Dirk Toewe on 20.06.17.
  */
package object opt
{
  private[opt] val UNDEFINED_ARRAY = Vec.zeros(0)

  private[opt] def isDefined  ( vec: Vec ) = UNDEFINED_ARRAY.vals ne vec.vals
  private[opt] def isUndefined( vec: Vec ) = UNDEFINED_ARRAY.vals eq vec.vals

  def regulaFalsi( f: Double => Double, lo: Double, hi: Double )
    = ???
}