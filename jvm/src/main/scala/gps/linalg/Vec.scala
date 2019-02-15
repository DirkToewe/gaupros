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

import gps.{linalg => la}

import scala.annotation.tailrec

/** A low-level 64-bit float vector implementation.
  *
  * The JVM implementation is backed by an Array[Double].
  *
  * Created by Dirk Toewe on 06.08.17.
  */
class Vec private[linalg]( val vals: Array[Double] ) extends AnyVal
{
  @inline def length = vals.length

  @inline def apply( idx: Int ) = vals(idx)

  @inline def update( idx: Int, value: Double ) = vals(idx) = value

  def isDefinedAt( idx: Int ): Boolean
    = 0 <= idx && idx < length

  @inline def foreach[U]( method: Double => U ) = {
    var    i = vals.length
    while( i >  0 ) {
           i -= 1
      method( vals(i) )
    }
  }

  @inline def foreach[U]( method: (Double,Int) => U ) = {
    var    i = vals.length
    while( i >  0 ) {
      i -= 1
      method( vals(i), i )
    }
  }

  @inline def forall( predicate: Double => Boolean ) = {
    @tailrec def loop( i: Int ): Boolean
      = i < 0 || predicate{this(i)} && loop(i-1)
    loop(length-1)
  }

  def clone = new Vec(vals.clone)

  def iterator = new Iterator[Double]{
    private var i=0
    def hasNext = i < length
    def    next = {
      if( ! hasNext ) throw new NoSuchElementException
      val result = apply(i); i += 1
      result
    }
  }

  def seq: IndexedSeq[Double] = vals

  def slice( from: Int, until: Int=length ) = {
    val len = until-from
    val result = new Array[Double](len)
    System.arraycopy(vals,from, result,0, len)
    new Vec(result)
  }

  def ⋅ ( v: Vec ) =
  {
    if( this.length != v.length )
      throw new IllegalArgumentException
    la.sum(length){ i => this(i) * v(i) }
  }

  def map( f: Double => Double ) = {
    val result = new Array[Double](length)
    var i = result.length
    while( 0 < i ) {
      i -= 1
      result(i) = f( this(i) )
    }
    new Vec(result)
  }

  def += ( subtrahend: Double ) = {
    var i = length
    while( i > 0 ) {
      i -= 1
      this(i) += subtrahend
    }
  }

  def -= ( summand: Double ) = {
    var i = length
    while( i > 0 ) {
      i -= 1
      this(i) -= summand
    }
  }

  def *= ( factor: Double ) = {
    var i = length
    while( i > 0 ) {
      i -= 1
      this(i) *= factor
    }
  }

  def /= ( divisor: Double ) = {
    var i = length
    while( i > 0 ) {
      i -= 1
      this(i) /= divisor
    }
  }

  def += ( summands: Vec ) =
  {
    if( summands.length != length )
      throw new IllegalArgumentException
    var i = length
    while( i > 0 ) {
      i -= 1
      this(i) += summands(i)
    }
  }

  def -= ( subtrahends: Vec ) =
  {
    if( subtrahends.length != length )
      throw new IllegalArgumentException
    var i = length
    while( i > 0 ) {
      i -= 1
      this(i) -= subtrahends(i)
    }
  }

  def *= ( factors: Vec ) =
  {
    if( factors.length != length )
      throw new IllegalArgumentException
    var i = length
    while( i > 0 ) {
      i -= 1
      this(i) *= factors(i)
    }
  }

  def /= ( divisors: Vec ) =
  {
    if( divisors.length != length )
      throw new IllegalArgumentException
    var i = length
    while( i > 0 ) {
      i -= 1
      this(i) /= divisors(i)
    }
  }

  def unary_-
    = map(-_)

  def + ( summand: Vec )
    = if( summand.length != length )
        throw new IllegalArgumentException
      else
        Vec.tabulate(length){ i => this(i) + summand(i) }

  def - ( subtrahend: Vec )
    = if( subtrahend.length != length )
        throw new IllegalArgumentException
      else
        Vec.tabulate(length){ i => this(i) - subtrahend(i) }

  def * ( factor: Vec )
    = if( factor.length != length )
        throw new IllegalArgumentException
      else
        Vec.tabulate(length){ i => this(i) * factor(i) }

  def / ( divisor: Vec )
    = if( divisor.length != length )
        throw new IllegalArgumentException
      else
        Vec.tabulate(length){ i => this(i) / divisor(i) }

  def + ( summand   : Double ) = map (_ + summand)
  def - ( subtrahend: Double ) = map (_ - subtrahend)
  def * ( factor    : Double ) = map (_ * factor)
  def / ( divisor   : Double ) = map (_ / divisor)

  def sum = la.sum(length)( this(_) )
  def min = {
    assert(length > 0)
    var result = Double.PositiveInfinity
    var i = length
    while( i > 0 ) {
      i -= 1
      if( ! {result <= this(i)} ) result = this(i) // <- handles NaN
    }
    result
  }
  def max = {
    assert(length > 0)
    var result = Double.NegativeInfinity
    var i = length
    while( i > 0 ) {
      i -= 1
      if( ! {result >= this(i)} ) result = this(i) // <- handles NaN
    }
    result
  }

  override def toString: String = mkString("[",", ","]")

  def mkString( infix: CharSequence ): String
    = mkString("",infix,"")

  def mkString( infix: CharSequence, format: Double => String ): String
    = mkString("",infix,"",format)

  def mkString( prefix: CharSequence, infix: CharSequence, suffix: CharSequence, format: Double => String = _.toString ): String =
  {
    val _infix = infix.toString
    val result = new StringBuilder(prefix.toString)
    var i = 0
    while( i < length-1 ) {
      result ++= format( this(i) )
      result ++= _infix
      i += 1
    }
    if( 0 < length )
      result ++= format( this(length-1) )
    result ++= suffix.toString
    result.toString
  }
}
object Vec
{
  @inline def apply( vals: Double* ) = new Vec( Array(vals: _*) )

  def tabulate( size: Int )( tabulator: Int => Double ) =
  {
    if( size < 0 )
      throw new IllegalArgumentException
    val result = new Array[Double](size)
    var i = result.length
    while( 0 < i ) {
      i -= 1
      result(i) = tabulator(i)
    }
    new Vec(result)
  }

  def zeros( size: Int )
    = new Vec( new Array[Double](size) )

  def ones( size: Int )
    = Vec.full(size)(1)

  def full( size: Int )( value: Double )
    = tabulate(size)( _ => value )

  class VecExtractor private[Vec]( val vec: Array[Double] ) extends AnyVal
  {
    def isEmpty = false
    def get = vec
  }

  def unapplySeq( vec: Vec ) = new VecExtractor(vec.vals)

  @inline def _copy( src: Vec, srcPos: Int, dest: Vec, destPos: Int, len: Int )
    = System.arraycopy(src.vals,srcPos, dest.vals,destPos, len)
}