package gps.linalg

import gps.{linalg => la}

import scala.collection.mutable

/** A low-level 64-bit float vector implementation.
  *
  * The JavaScript implementation is backed by a Float64Array.
  *
  * Created by Dirk Toewe on 06.08.17.
  */
class Vec private[linalg]( val vals: Float64Array ) extends AnyVal
{
  @inline def length: Int = vals.length

  // apply *should* throw an error on out of bounds indicies since undefined is no Double...
  @inline def apply( idx: Int ): Double = {
    if( 0 > idx || idx >= length ) throw new IndexOutOfBoundsException
    vals(idx)
  }

  @inline def update( idx: Int, value: Double ): Unit = {
    if( 0 > idx || idx >= length ) throw new IndexOutOfBoundsException
    vals(idx) = value
  }

  def isDefinedAt( idx: Int ): Boolean
    = 0 <= idx && idx < length

  def forall( predicate: Double => Boolean ): Boolean = vals.every(predicate)

  def clone = new Vec( new Float64Array(vals) )

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

  def seq: IndexedSeq[Double] = new mutable.IndexedSeq[Double]{
    override def length = vals.length
    override def apply( idx: Int ) = vals(idx)
    override def update( idx: Int, value: Double ) = vals(idx) = value
  }

  def slice( from: Int, until: Int=length ) = new Vec( vals.slice(from,until) )

  def ⋅ ( v: Vec ) = {
    if( v.length != length )
      throw new IllegalArgumentException
    la.sum(length){ i => vals(i) * v(i) }
  }

  @inline def map( mapFn:  Double      => Double ) = new Vec(vals map mapFn)
  @inline def map( mapFn: (Double,Int) => Double ) = new Vec(vals map mapFn)

  def += (   summand: Double ) = { var i = length; while( i > 0 ) { i -= 1; vals(i) += summand    } }
  def -= (subtrahend: Double ) = { var i = length; while( i > 0 ) { i -= 1; vals(i) -= subtrahend } }
  def *= (    factor: Double ) = { var i = length; while( i > 0 ) { i -= 1; vals(i) *= factor     } }
  def /= (   divisor: Double ) = { var i = length; while( i > 0 ) { i -= 1; vals(i) /= divisor    } }

  def += ( summands: Vec ) =
  {
    if( summands.length != length )
      throw new IllegalArgumentException
    var i = length
    while( i > 0 ) {
      i -= 1
      vals(i) += summands(i)
    }
  }

  def -= ( subtrahends: Vec ) =
  {
    if( subtrahends.length != length )
      throw new IllegalArgumentException
    var i = length
    while( i > 0 ) {
      i -= 1
      vals(i) -= subtrahends(i)
    }
  }

  def *= ( factors: Vec ) =
  {
    if( factors.length != length )
      throw new IllegalArgumentException
    var i = length
    while( i > 0 ) {
      i -= 1
      vals(i) *= factors(i)
    }
  }

  def /= ( divisors: Vec ) =
  {
    if( divisors.length != length )
      throw new IllegalArgumentException
    var i = length
    while( i > 0 ) {
      i -= 1
      vals(i) /= divisors(i)
    }
  }

  @inline def unary_- = map(-_)

  @inline def + (   summands: Vec ) = if(   summands.length == length ) map{_ +    summands(_)} else throw new IllegalArgumentException
  @inline def - (subtrahends: Vec ) = if(subtrahends.length == length ) map{_ - subtrahends(_)} else throw new IllegalArgumentException
  @inline def * (    factors: Vec ) = if(    factors.length == length ) map{_ *     factors(_)} else throw new IllegalArgumentException
  @inline def / (   divisors: Vec ) = if(   divisors.length == length ) map{_ /    divisors(_)} else throw new IllegalArgumentException

  @inline def + (    summand: Double ) = map (_ + summand   )
  @inline def - ( subtrahend: Double ) = map (_ - subtrahend)
  @inline def * (     factor: Double ) = map (_ * factor    )
  @inline def / (    divisor: Double ) = map (_ / divisor   )

  def sum = la.sum(length)( vals(_) )

  @inline def min = vals.reduce[Double]( math.min(_,_) )
  @inline def max = vals.reduce[Double]( math.max(_,_) )

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
  @inline def apply( vals: Double* ) = new Vec( Float64Array.of(vals: _*) )

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

  def zeros( size: Int ) = new Vec( new Float64Array(size) )
  def  ones( size: Int ) =     Vec.full(size)(1)

  def full( size: Int )( value: Double )
    = new Vec( new Float64Array(size).fill(value) )

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