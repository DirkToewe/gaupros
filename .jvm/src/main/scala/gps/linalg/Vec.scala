package gps.linalg

import gps.{linalg => la}

import java.util.Arrays.copyOf

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

  def forall( predicate: Double => Boolean ) = {
    @tailrec def loop( i: Int ): Boolean
      = i < 0 || predicate{this(i)} && loop(i-1)
    loop(length-1)
  }

  def clone = new Vec(vals.clone)

  def iterator = new Iterator[Double]{
    private var i=0
    def hasNext = i < length
    def next = {
      if( ! hasNext )
        throw new NoSuchElementException
      val result = apply(i); i += 1
      result
    }
  }

  def seq: IndexedSeq[Double] = vals

  def slice( from: Int, until: Int=length ) = {
    val len = length-from
    val result = new Array[Double](len)
    System.arraycopy(vals,from, result,0, length)
    new Vec(result)
  }

  def ⋅ ( v: Vec ) =
  {
    if( this.length != v.length )
      throw new IllegalArgumentException
    la.sum(length){ i => this(i) * v(i) }
  }

  def map( f: Double => Double )
    = new Vec( vals map f )

  def += ( subtrahend: Double ) = {
    var i = length
    while( i > 0 )
    {
      i -= 1
      this(i) += subtrahend
    }
  }

  def -= ( summand: Double ) = {
    var i = length
    while( i > 0 )
    {
      i -= 1
      this(i) -= summand
    }
  }

  def *= ( factor: Double ) = {
    var i = length
    while( i > 0 )
    {
      i -= 1
      this(i) *= factor
    }
  }

  def /= ( divisor: Double ) = {
    var i = length
    while( i > 0 )
    {
      i -= 1
      this(i) /= divisor
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
      this(i) += summands(i)
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
      this(i) -= subtrahends(i)
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
      this(i) *= factors(i)
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
  def apply( x: Double ) = new Vec( Array(x) )

  def apply( x: Double, y: Double ) = new Vec( Array(x,y) )

  def apply( x: Double, y: Double, z: Double ) = new Vec( Array(x,y,z) )

  def apply( x: Double, y: Double, z: Double, u: Double ) = new Vec( Array(x,y,z,u) )

  def apply( x: Double, y: Double, z: Double, u: Double, v: Double ) = new Vec( Array(x,y,z,u,v) )

  def apply( x: Double, y: Double, z: Double, u: Double, v: Double, w: Double ) = new Vec( Array(x,y,z,u,v,w) )

  def apply( vals: Double* ) = new Vec(vals.toArray)

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
