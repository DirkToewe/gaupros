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

package gps

import java.util.Arrays.copyOf

import scala.collection.immutable.{IndexedSeq => ISeq}
import scala.language.implicitConversions

/**
  * Created by dtub on 27.07.17.
  */
package object kernel
{
  // implicit def double2kernel[D <% Double]( const: D ): Kernel[Any] = Const(const)
  implicit def double2kernel( const: Double ): Kernel[Any] = Const(const)
  implicit def int2kernel( const: Int ): Kernel[Any] = Const(const)
  implicit def symbol2sym( sym: Symbol ): Kernel[Any] = Sym(sym)

  private[kernel] val NO_PARAMS = ISeq.empty[Symbol]
  private[kernel] implicit val SYMBOL_ORDER: Ordering[Symbol] = Ordering by (_.name)
  private[kernel] def sortedUnion( a: IndexedSeq[Symbol], b: IndexedSeq[Symbol] ): IndexedSeq[Symbol] =
  {
    // use merge sort like merging.
    var ia = 0
    var ib = 0
    val result = new Array[Symbol](a.size + b.size)
    var len = 0; var i = 0
    while( ia < a.size && ib < b.size )
    {
      val c = implicitly[Ordering[Symbol]].compare(a{ia}, b{ib})
      result(i) = if( c < 0 ) a(ia) else b(ib)
      if( c <= 0 ) ia += 1
      if( c >= 0 ) ib += 1
      len += 1; i += 1
    }
    while( ia < a.size ) { result(i) = a(ia); len += 1; i += 1; ia += 1 }
    while( ib < b.size ) { result(i) = b(ib); len += 1; i += 1; ib += 1 }
    assert( ia == a.size )
    assert( ib == b.size )
    assert( i == len )
    copyOf(result,len)
  }
}