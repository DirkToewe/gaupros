package gps.linalg

import gps.{linalg => la}

import scala.annotation.tailrec
import scala.scalajs.js.typedarray.Float64Array

/** A low-level 64-bit float vector implementation.
  *
  * The JavaScript implementation is backed by a Float64Array.
  *
  * Created by Dirk Toewe on 06.08.17.
  */
class Vec private[linalg]( val vals: Float64Array ) extends AnyVal
{
  def length = vals.length

  def apply( idx: Int ) = {
    if( 0 > idx || idx >= length ) throw new IndexOutOfBoundsException
    vals(idx)
  }

  def update( idx: Int, value: Double ) = {
    if( 0 > idx || idx >= length ) throw new IndexOutOfBoundsException
    vals(idx) = value
  }

  def isDefinedAt( idx: Int ): Boolean
    = 0 <= idx && idx < length

  def forall( predicate: Double => Boolean ) = {
    @tailrec def loop( i: Int ): Boolean
      = i < 0 || predicate{vals(i)} && loop(i-1)
    loop(length-1)
  }

  def clone = Vec.tabulate(length){apply}

  def iterator = new Iterator[Double]{
    private var i=0
    def hasNext = i < length
    def next = {
      if( ! hasNext )
        throw new NoSuchElementException
      val result = vals(i); i += 1
      result
    }
  }

  def seq: IndexedSeq[Double] = vals.toIndexedSeq

  def slice( from: Int, until: Int=length )
    = Vec.tabulate(until-from){ i => vals(i + from) }

  def ⋅ ( v: Vec ) =
  {
    if( v.length != length )
      throw new IllegalArgumentException
    la.sum(length){ i => vals(i) * v(i) }
  }

  def map( f: Double => Double )
  = Vec.tabulate(length){ f apply vals(_) }

  def += ( subtrahend: Double ) = {
    var i = length
    while( i > 0 )
    {
      i -= 1
      vals(i) += subtrahend
    }
  }

  def -= ( summand: Double ) = {
    var i = length
    while( i > 0 )
    {
      i -= 1
      vals(i) -= summand
    }
  }

  def *= ( factor: Double ) = {
    var i = length
    while( i > 0 )
    {
      i -= 1
      vals(i) *= factor
    }
  }

  def /= ( divisor: Double ) = {
    var i = length
    while( i > 0 )
    {
      i -= 1
      vals(i) /= divisor
    }
  }

  def += ( summands: Vec ) =
  {
    if( summands.length != length )
      throw new IllegalArgumentException
    var i = length
    while( i > 0 )
    {
      i -= 1
      vals(i) += summands(i)
    }
  }

  def -= ( subtrahends: Vec ) =
  {
    if( subtrahends.length != length )
      throw new IllegalArgumentException
    var i = length
    while( i > 0 )
    {
      i -= 1
      vals(i) -= subtrahends(i)
    }
  }

  def *= ( factors: Vec ) =
  {
    if( factors.length != length )
      throw new IllegalArgumentException
    var i = length
    while( i > 0 )
    {
      i -= 1
      vals(i) *= factors(i)
    }
  }

  def /= ( divisors: Vec ) =
  {
    if( divisors.length != length )
      throw new IllegalArgumentException
    var i = length
    while( i > 0 )
    {
      i -= 1
      vals(i) /= divisors(i)
    }
  }

  def unary_-
    = map(-_)

  def + ( summand: Vec )
    = if( summand.length != length )
        throw new IllegalArgumentException
      else
        Vec.tabulate(length){ i => vals(i) + summand.vals(i) }

  def - ( subtrahend: Vec )
    = if( subtrahend.length != length )
        throw new IllegalArgumentException
      else
        Vec.tabulate(length){ i => vals(i) - subtrahend.vals(i) }

  def * ( factor: Vec )
    = if( factor.length != length )
        throw new IllegalArgumentException
      else
        Vec.tabulate(length){ i => vals(i) * factor.vals(i) }

  def / ( divisor: Vec )
    = if( divisor.length != length )
        throw new IllegalArgumentException
      else
        Vec.tabulate(length){ i => vals(i) / divisor.vals(i) }

  def + ( summand   : Double ) = map (_ + summand)
  def - ( subtrahend: Double ) = map (_ - subtrahend)
  def * ( factor    : Double ) = map (_ * factor)
  def / ( divisor   : Double ) = map (_ / divisor)

  def sum = la.sum(length)( vals(_) )

  def min = {
    assert(length > 0)
    var result = Double.PositiveInfinity
    var i = length
    while( i > 0 ) {
      i -= 1
      if( ! {result <= vals(i)} ) result = vals(i) // <- handles NaN
    }
    result
  }

  def max = {
    assert(length > 0)
    var result = Double.NegativeInfinity
    var i = length
    while( i > 0 ) {
      i -= 1
      if( ! {result >= vals(i)} ) result = vals(i) // <- handles NaN
    }
    result
  }

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
      result ++= format( vals(i) )
      result ++= _infix
      i += 1
    }
    if( 0 < length )
      result ++= format( vals(length-1) )
    result ++= suffix.toString
    result.toString
  }
}
object Vec
{
  def apply( x: Double ) = {
    val result = Vec.zeros(1)
    result(0) = x
    result
  }

  def apply( x: Double, y: Double ) = {
    val result = Vec.zeros(2)
    result(0) = x
    result(1) = y
    result
  }

  def apply( x: Double, y: Double, z: Double ) = {
    val result = Vec.zeros(3)
    result(0) = x
    result(1) = y
    result(2) = z
    result
  }

  def apply( x: Double, y: Double, z: Double, u: Double ) = {
    val result = Vec.zeros(4)
    result(0) = x
    result(1) = y
    result(2) = z
    result(3) = u
    result
  }

  def apply( x: Double, y: Double, z: Double, u: Double, v: Double ) = {
    val result = Vec.zeros(5)
    result(0) = x
    result(1) = y
    result(2) = z
    result(3) = u
    result(4) = v
    result
  }

  def apply( x: Double, y: Double, z: Double, u: Double, v: Double, w: Double ) = {
    val result = Vec.zeros(6)
    result(0) = x
    result(1) = y
    result(2) = z
    result(3) = u
    result(4) = v
    result(5) = w
    result
  }

  def apply( vals: Double* ) = {
    val _vals = if( vals.isInstanceOf[IndexedSeq[_]] ) vals else vals.toIndexedSeq
    tabulate(_vals.length)(_vals)
  }

  def tabulate( size: Int )( tabulator: Int => Double ) =
  {
    if( size < 0 )
      throw new IllegalArgumentException
    val result = new Float64Array(size)
    var i = result.length
    while( 0 < i ) {
      i -= 1
      result(i) = tabulator(i)
    }
    new Vec(result)
  }

  def zeros( size: Int )
    = new Vec( new Float64Array(size) )

  def ones( size: Int )
    = Vec.full(size)(1)

  def full( size: Int )( value: Double )
    = tabulate(size)( _ => value )

  class VecExtractor private[Vec]( val vec: Float64Array ) extends AnyVal
  {
    def isEmpty = false
    def get = vec
  }

  def unapplySeq( vec: Vec ) = new VecExtractor(vec.vals)

  @inline def _copy( src: Vec, srcPos: Int, dst: Vec, dstPos: Int, len: Int ): Unit = {
    if( srcPos < 0 || dstPos < 0 || len < 0 )
      throw new IndexOutOfBoundsException

    if(    srcPos > src.length - len
        || dstPos > dst.length - len )
      throw new IndexOutOfBoundsException

    var i=len
    while( i > 0 ) {
      i -= 1
      dst(i+dstPos) = src(i+srcPos)
    }
  }
}