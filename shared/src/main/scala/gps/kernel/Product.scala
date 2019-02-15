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
private[kernel] class Product[-X] private(val a: Kernel[X], val b: Kernel[X] ) extends BinOp[X](a,b)
{
  override def apply( xi: X, xj: X, i: Int, j: Int ): Double
    = a(xi,xj, i,j) *
      b(xi,xj, i,j)

  protected override def eval[Y]( a: Kernel[Y], b: Kernel[Y] ): Kernel[Y]
    = a * b

  protected override def diff[Y <: X]( da: Kernel[Y], db: Kernel[Y] ): Kernel[Y]
    = da*b + a*db

  override protected[kernel] def toString(xi: String, xj: String ) =
  {
    def parenthesize( k: Kernel[X] ) =
    {
      val str = k.toString(xi,xj)
      k match {
        case (_) + (_) => s"($str)"
        case _ => str
      }
    }
    val A = parenthesize(a)
    val B = parenthesize(b)
    s"$A * $B"
  }
}
private[kernel] object Product
{
  def apply[X]( factor1: Kernel[X], factor2: Kernel[X] ): Kernel[X]
    = (factor1,factor2) match
      {
        case (      a, Const(1)) => a
        case (Const(a),Const(b) *> c ) => a*b * c
        // bring all constants to the very left
        case (      a, Const(b) *> c ) => b * (a*c)
        case (Const(1), kernel ) => kernel
        case (Const(0),      _ ) => 0
        case (Noise(a),      b ) => Noise(a*b)
        case (      a, Noise(b)) => Noise(a*b)
        case (a*b, c) => a * (b*c)
        // sort symbols alphabetically
        case (a: Const,b) => new Product(a,b)
        case ( Pow1(Sym(a1),a2), Pow1(Sym(b1),b2) *> c ) =>
               if( a1   ==   b1      )  a1.pow(a2+b2) * c
          else if( a1.name > b1.name ) (b1 pow b2) * (a1 pow a2) * c
          else new Product(factor1,factor2)
        case (a, ( b @ Pow1(Sym(_),_) ) *> c ) => new Product(b,a*c)
        case (a,b) => new Product(a,b)
      }
}
object *
{
  def unapply[X]( product: Product[X] )
    = Some(product.a, product.b)
}
protected[kernel] object *>
{
  /** Always matches an addition, inserting the identity (1) on the right if need be.
    *
    * @param kernel
    * @tparam X
    * @return
    */
  def unapply[X]( kernel: Kernel[X] )
    = kernel match {
        case a*b => Some(a,b)
        case _ => Some(kernel, Const(1))
      }
}
protected[kernel] object *<
{
  /** Always matches an addition, inserting the identity (1) on the left if need be.
    *
    * @param kernel
    * @tparam X
    * @return
    */
  def unapply[X]( kernel: Kernel[X] )
    = kernel match {
        case a*b => Some(a,b)
        case _ => Some(Const(1), kernel)
      }
}