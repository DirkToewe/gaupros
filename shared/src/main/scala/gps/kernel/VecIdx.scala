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

import gps.linalg.Vec

/** A Kernel[Vec] that applies a Kernel[Double] on a specific index entry of a the argument Vecs.
  *
  * Created by Dirk Toewe on 10.08.17.
  */
protected[kernel] class VecIdx private[kernel]( val idx: Int, val op: Kernel[Double] ) extends Kernel[Vec]
{
  assert( idx >= 0 )

  override def apply( xi: Vec, xj: Vec )
    = this(xi,xj, -1,-1)

  override def apply( xi: Vec, xj: Vec, i: Int, j: Int ) =
  {
    assert(xi.length > idx)
    assert(xi.length == xj.length)
    op(xi(idx),xj(idx), i,j)
  }

  override def subs[Y <: Vec]( map: PartialFunction[Symbol,Kernel[Y]] ): Kernel[Y] =
  {
    val _map: PartialFunction[Symbol,Kernel[Double]] = {
      case sym =>
        val result = map applyOrElse (sym, (_: Any) => sym: Kernel[Double] )
        if( result.argDependent )
          throw new IllegalArgumentException( s"Cannot substitute '$sym' by '$result' since it is part of a nested kernel." )
        result.asInstanceOf[Kernel[Double]]
    }
    VecIdx(idx, op subs _map)
  }

  override def pDiff(sym: Symbol )
    = VecIdx(idx, op pDiff sym)

  override val params = op.params

  override protected[kernel] def toString(xi: String, xj: String ) =
  {
    val Xi = s"${xi.substring(0, xi.length-1)},$idx]" ensuring  (xi endsWith "]")
    val Xj = s"${xj.substring(0, xj.length-1)},$idx]" ensuring  (xj endsWith "]")
    op.toString(Xi,Xj)
  }
}
object VecIdx
{
  def apply( idx: Int, op: Kernel[Double] ): Kernel[Vec]
    = op match {
        case Independent(a)     => a
        case Independent(a) + b => a + VecIdx(idx,b)
        case Independent(a) * b => a * VecIdx(idx,b)
        case Pow(a,Independent(b)) =>  VecIdx(idx,a) pow b
        case _ => new VecIdx(idx,op)
      }

  def unapply( vecIndex: VecIdx )
    = Some(vecIndex.idx, vecIndex.op)
}