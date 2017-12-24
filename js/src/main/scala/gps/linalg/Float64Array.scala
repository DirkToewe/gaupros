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