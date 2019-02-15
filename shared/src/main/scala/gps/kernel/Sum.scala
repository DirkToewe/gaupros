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
  * Created by dtitx on 30.07.17.
  */
private[kernel] class Sum[X]( val a: Kernel[X], val b: Kernel[X] ) extends BinOp[X](a,b)
{
  override def apply( xi: X, xj: X, i: Int, j: Int ): Double
    = a(xi,xj, i,j) +
      b(xi,xj, i,j)

  protected override def eval[Y]( a: Kernel[Y], b: Kernel[Y] ): Kernel[Y]
    = a + b

  protected override def diff[Y <: X]( da: Kernel[Y], db: Kernel[Y] ): Kernel[Y]
    = da + db

  override protected[kernel] def toString(xi: String, xj: String ) =
  {
    val A = a.toString(xi, xj)
    val B = b.toString(xi, xj)
    s"$A + $B"
  }
}
private[kernel] object Sum
{
  def apply[X]( summand1: Kernel[X], summand2: Kernel[X] ): Kernel[X]
    = (summand1,summand2) match
      {
        case (      a ,Const(0)) => a
        case (Const(a),Const(b) +> c) => a+b + c
        case (Const(0),      b ) => b
        case ( a+b, c ) => a + (b+c)
        // move constants to the very left
        case ( a, Const(b) +> c ) => b + (a+c)
        case ( a1 *< Sym(a2), b1 *< Sym(b2) +> c ) =>
          if( a2 == b2 ) (a1+b1) * a2 + c
          else if( b2.name < a2.name ) b1*b2 + a1*a2 + c
          else new Sum(summand1,summand2)
        case ( Noise(a), Noise(b) ) => Noise(a+b)
        case (a,b) => new Sum(a,b)
      }
}
object +
{
  def unapply[X]( sum: Sum[X] )
    = Some(sum.a, sum.b)
}
object +>
{
  /** Always matches and addition adding a zero to the right if need be.
    *
    * @param kernel
    * @tparam X
    * @return
    */
  def unapply[X]( kernel: Kernel[X] )
    = kernel match {
        case a+b => Some(a,b)
        case _ => Some(kernel,Const(0))
      }
}
object +<
{
  /** Always matches and addition adding a zero to the left if need be.
    *
    * @param kernel
    * @tparam X
    * @return
    */
  def unapply[X]( kernel: Kernel[X] )
    = kernel match {
        case a+b => Some(a,b)
        case _ => Some(Const(0),kernel)
      }
}
