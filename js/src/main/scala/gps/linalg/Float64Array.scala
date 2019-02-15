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

package gps.linalg

import scala.scalajs.js
import scala.scalajs.js.annotation.{JSBracketAccess, JSGlobal}

@js.native
@JSGlobal
class Float64Array private[this]() extends js.Any
{
  def this( length: Int ) = this()
  def this( length: Float64Array ) = this()

  val length: Int = js.native

  @JSBracketAccess def  apply( index: Int ): Double = js.native
  @JSBracketAccess def update( index: Int, value: Double ): Unit = js.native

  def map( mapFn: js.Function1[Double,    Double] ): Float64Array = js.native
  def map( mapFn: js.Function2[Double,Int,Double] ): Float64Array = js.native

  def every( predicate: js.Function1[Double,    Boolean] ): Boolean = js.native
  def every( predicate: js.Function2[Double,Int,Boolean] ): Boolean = js.native

  def fill( value: Double, from: Int=0, until: Int=length ): Float64Array = js.native

  def reduce[R]( reduceFn: js.Function2[R,Double,R], initialValue: Int = ??? ): R = js.native

  def slice( from: Int, until: Int=length ): Float64Array = js.native
}
@js.native
@JSGlobal
object Float64Array extends js.Any
{
  def of( values: Double* ): Float64Array = js.native

  def from( source: Float64Array, mapFn: js.Function1[Double,    Double] ): Float64Array = js.native
  def from( source: Float64Array, mapFn: js.Function2[Double,Int,Double] ): Float64Array = js.native
}