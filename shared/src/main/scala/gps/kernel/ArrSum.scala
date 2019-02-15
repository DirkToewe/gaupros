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

/** A Kernel[Array[X]] that uses a nested Kernel[X] to zip the entries of the array
  * arguments pairwise and sums them up afterwards.
  *
  * Created by Dirk Toewe on 30.07.17.
  */
private[kernel] class ArrSum[X](val mapOp: Kernel[X] ) extends Kernel[Array[X]]
{
  override def apply( xi: Array[X], xj: Array[X] )
    = this(xi,xj, -1,-1)

  override def apply( xi: Array[X], xj: Array[X], i: Int, j: Int ) =
  {
    assert(xi.length > 0)
    assert(xi.length == xj.length)

    def distPow( k: Int, len: Int ): Double
      = if( len == 1 )
          mapOp(xi{k},xj{k}, i,j)
        else {
          val l = len / 2
          distPow(k,l) +
          distPow(k+l, len-l)
        }

    distPow(0, xi.length)
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
    ArrSum(mapOp subs _map)
  }

  override def pDiff(sym: Symbol )
    = ArrSum(mapOp pDiff sym)

  override val params = mapOp.params

  override protected[kernel] def toString(xi: String, xj: String ) =
  {
    val Xi = xi.substring(0, xi.length-1) + ",k]"  ensuring  (xi endsWith "]")
    val Xj = xj.substring(0, xj.length-1) + ",k]"  ensuring  (xj endsWith "]")
    s"Σₖ{ ${mapOp.toString(Xi,Xj)} }"
  }
}
object ArrSum
{
  def apply[X]( mergeOp: Kernel[X] ): Kernel[Array[X]]
    = mergeOp match {
        case Independent(a) * b => a * ArrSum(b)
        case _ => new ArrSum(mergeOp)
      }

  def unapply[X]( argSum: ArrSum[X] )
    = Some(argSum.mapOp)
}