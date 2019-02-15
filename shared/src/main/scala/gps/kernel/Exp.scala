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
  * Created by Dirk Toewe on 30.07.17.
  */
private[kernel] class Exp[-X] private(val e: Kernel[X] ) extends Func[X](e)
{
  override def apply( xi: X, xj: X, i: Int, j: Int )
    = math.exp{ e(xi,xj, i,j) }

  override def pDiff( sym: Symbol ): Kernel[X]
    = this * (e pDiff sym)

  override def subs[Y <: X](map: PartialFunction[Symbol,Kernel[Y]] ): Kernel[Y]
    = Exp(e subs map)
}
object Exp
{
  def apply[X]( exponent: Kernel[X] ): Kernel[X]
    = exponent match {
        case Log(kernel) => kernel
        case Const(e) => math.exp(e)
        case kernel => new Exp(kernel)
      }

  def unapply[X]( exp: Exp[X] ) = Some(exp.e)
}