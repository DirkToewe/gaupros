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

package gps.kernel

/**
  * Created by Dirk Toewe on 05.08.17.
  */
protected[kernel] abstract class Func[-X](arg: Kernel[X] ) extends Kernel[X]
{
  override def apply( xi: X, xj: X )
    = this(xi,xj, -1,-1)

  override def apply( xi: X, xj: X, i: Int, j: Int ): Double

  override protected[kernel] final def toString(xi: String, xj: String )
    = s"${getClass.getSimpleName}(${arg.toString(xi,xj)})"

  override protected[kernel] final val argDependent = arg.argDependent
  override final val params = arg.params
}
