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

/** A sub-kernel of an kernel with array input. The sub-kernel is only applied to one index of
  * the array input pair.
  *
  * Created by Dirk Toewe on 10.08.17.
  */
protected[kernel] class ArrIdx[@specialized X] private[kernel](val idx: Int, val mergeOp: Kernel[X] ) extends Kernel[Array[X]]
{
  assert( idx >= 0 )

  override def apply( xi: Array[X], xj: Array[X] )
    = this(xi,xj, -1,-1)

  override def apply( xi: Array[X], xj: Array[X], i: Int, j: Int ) =
  {
    assert(xi.length > idx)
    assert(xi.length == xj.length)
    mergeOp(xi(idx),xj(idx), i,j)
  }

  override def subs[Y <: Array[X]]( map: PartialFunction[Symbol,Kernel[Y]] ): Kernel[Y] =
  {
    val _map: PartialFunction[Symbol,Kernel[X]] = {
      case sym =>
        val result = map applyOrElse (sym, (_: Any) => sym: Kernel[X] )
        if( result.argDependent )
          throw new IllegalArgumentException( s"Cannot substitute '$sym' by '$result' since it is part of a nested kernel." )
        result.asInstanceOf[Kernel[X]]
    }
    ArrIdx(idx, mergeOp subs _map)
  }

  override def pDiff(sym: Symbol )
  = ArrIdx(idx, mergeOp pDiff sym)

  override val params = mergeOp.params

  override protected[kernel] def toString(xi: String, xj: String ) =
  {
    val Xi = s"${xi.substring(0, xi.length-1)},$idx]" ensuring  (xi endsWith "]")
    val Xj = s"${xj.substring(0, xj.length-1)},$idx]" ensuring  (xj endsWith "]")
    mergeOp.toString(Xi,Xj)
  }
}
object ArrIdx
{
  def apply[X]( idx: Int, mergeOp: Kernel[X] ): Kernel[Array[X]]
    = mergeOp match {
        case Independent(a)     => a
        case Independent(a) + b => a + ArrIdx(idx,b)
        case Independent(a) * b => a * ArrIdx(idx,b)
        case Pow(a,Independent(b)) =>  ArrIdx(idx,a) pow b
        case _ => new ArrIdx(idx,mergeOp)
      }

  def unapply[X]( ofArrElem: ArrIdx[X] )
    = Some(ofArrElem.idx, ofArrElem.mergeOp)
}