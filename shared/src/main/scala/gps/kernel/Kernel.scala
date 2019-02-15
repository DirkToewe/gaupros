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

import scala.collection.immutable.IndexedSeq

/** A base trait for parametrized and unparametrized machine learning Kernels.
  * The kernels are derivable w.r.t. their parameters using the pDiff function.
  *
  * Created by Dirk Toewe on 17.07.17.
  */
trait Kernel[@specialized -X]
{
  def apply( xi: X, xj: X ): Double

  def apply( xi: X, xj: X, i: Int, j: Int ): Double
    = apply(xi,xj)

  /** A sorted, read-only list of the names of all hyper-parameters of the kernel.
    */
  def params: IndexedSeq[Symbol] = NO_PARAMS

  /** Returns false if it is guaranteed to be constant wrt. the inputs xi and xj.
    *
    * @return false if it is guaranteed to be constant wrt. the inputs xi and xj.
    */
  protected[kernel] def argDependent: Boolean = true

  /** Returns the partial derivative of this Kernel in the given hyper-parameter.
    *
    * @param sym The name of the hyper-parameter.
    * @return d(this)/d(sym)
    */
  def pDiff( sym: Symbol ): Kernel[X] = 0

  /** Substitutes hyper-parameters with the kernels given by the partial function.
    *
    * @param map The partial function defining the replaced symbols and their replacements.
    * @tparam Y The new input type of the resulting kernel.
    * @return
    */
  def subs[Y <: X]( map: PartialFunction[Symbol,Kernel[Y]] ): Kernel[Y] = this

  final def subs[Y <: X, K <% Kernel[Y]]( entries: (Symbol,K)* ): Kernel[Y] =
  {
    val kv = for( (k,v) <- entries ) yield k -> (v: Kernel[Y])
    subs{ Map(kv:_*) }
  }

  override final    def toString = "(X[i],X[j]) => " + toString("X[i]", "X[j]")
  protected[kernel] def toString(xi: String, xj: String ): String = "?"

  def pow[Y <: X]( exponent: Kernel[Y]): Kernel[Y] = Pow(this,exponent)
  def pow        ( exponent: Symbol   ): Kernel[X] = pow[X](exponent)
  def pow        ( exponent: Double   ): Kernel[X] = pow[X](exponent)

  def unary_+ : Kernel[X] =          this
  def unary_- : Kernel[X] = Negative(this)

  def + [Y <: X]( summand   : Kernel[Y] ): Kernel[Y] =        Sum(this,summand)
  def - [Y <: X]( subtrahend: Kernel[Y] ): Kernel[Y] = Difference(this,subtrahend)
  def * [Y <: X]( factor    : Kernel[Y] ): Kernel[Y] =    Product(this,factor)
  def / [Y <: X]( divisor   : Kernel[Y] ): Kernel[Y] =   Quotient(this,divisor)

  def + ( summand   : Symbol ): Kernel[X] =        Sum(this,summand)
  def - ( subtrahend: Symbol ): Kernel[X] = Difference(this,subtrahend)
  def * ( factor    : Symbol ): Kernel[X] =    Product(this,factor)
  def / ( divisor   : Symbol ): Kernel[X] =   Quotient(this,divisor)

  def + ( summand   : Double ): Kernel[X] =        Sum(this,summand)    
  def - ( subtrahend: Double ): Kernel[X] = Difference(this,subtrahend)
  def * ( factor    : Double ): Kernel[X] =    Product(this,factor)     
  def / ( divisor   : Double ): Kernel[X] =   Quotient(this,divisor)    
}