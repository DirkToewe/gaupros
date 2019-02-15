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
private[kernel] class Noise[-X] private(val m: Kernel[X] ) extends Kernel[X]
{
  override def apply( xi: X, xj: X )
    = 0

  override def apply( xi: X, xj: X, i: Int, j: Int )
    = if( 0 > i || i != j ) 0 else m(xi,xj, i,j)

  override def pDiff( sym: Symbol ): Kernel[X]
    = Noise(m pDiff sym)

  override def subs[Y <: X](map: PartialFunction[Symbol,Kernel[Y]] ): Kernel[Y]
    = Noise(m subs map)

  override protected[kernel] def toString(xi: String, xj: String ) = s"Noise(${ m.toString(xi,xj) })"
  override protected[kernel] val argDependent = m.argDependent
  override val params = m.params
}
object Noise
{
  def apply[X]( magnitude: Kernel[X] ): Kernel[X]
    = magnitude match {
        case Const(0) => 0
        case Noise(Noise(k))     => Noise(k)
        case Noise(a) *       b  => Noise(a*b)
        case Noise(a) + Noise(b) => Noise(a+b)
        case kernel => new Noise(kernel)
      }

  def unapply[X]( noise: Noise[X] )
    = Some(noise.m)
}