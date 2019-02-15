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

import scala.math.{Pi => π}

/**
  * Created by dtitx on 05.08.17.
  */
class Sin[@specialized -X]( k: Kernel[X] ) extends Func[X](k)
{
  override def apply( xi: X, xj: X, i: Int, j: Int )
    = math.sin{ k(xi,xj, i,j) }

  override def pDiff(sym: Symbol ): Kernel[X]
    = Sin(π/2 + k) * (k pDiff sym)

  override def subs[Y <: X]( map: PartialFunction[Symbol,Kernel[Y]] ): Kernel[Y]
    = Sin(k subs map)
}
object Sin
{
  def apply[X]( kernel: Kernel[X] ): Kernel[X]
    = kernel match {
        case Const(a) => math.sin(a)
        case Const(a) + b =>
          var angle = a % (2*π)
          if( angle >   π ) angle -= 2*π
          if( angle <= -π ) angle += 2*π
          new Sin( angle + b )
        case _ => new Sin(kernel)
      }
}