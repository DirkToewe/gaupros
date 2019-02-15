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
  * Created by Dirk Toewe on 31.07.17.
  */
private[kernel] class Log[@specialized -X] private(val k: Kernel[X] ) extends Func[X](k)
{
  override def apply( xi: X, xj: X, i: Int, j: Int )
    = math.log{ k(xi,xj, i,j) }

  override def pDiff(sym: Symbol ): Kernel[X]
    = (k pDiff sym) / k

  override def subs[Y <: X](map: PartialFunction[Symbol,Kernel[Y]] ): Kernel[Y]
    = Log(k subs map)
}
object Log
{
  def apply[X]( kernel: Kernel[X] ): Kernel[X]
    = kernel match {
        case Exp(a) => a
        case Const(e) => math.log(e)
        case Pow(b,Const(1)) => new Log(b)
        case Pow(b,e) => e*Log(b)
        case _ => new Log(kernel)
      }

  def unapply[X]( log: Log[X] ) = Some(log.k)
}