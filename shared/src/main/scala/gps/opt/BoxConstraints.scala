package gps.opt

import gps.linalg.Vec

import scala.Double.{PositiveInfinity => ∞}
import scala.annotation.tailrec

/** A utility class for handling box constraints on optimization.
  *
  * @param x_min
  * @param x_max
  */
class BoxConstraints private( val x_min: Vec, val x_max: Vec )
{
  assert(x_min.length == x_max.length)
  def nVars = x_min.length

  /** Is the given point is outside the box constraints, this method changes its
    * entries to the coordinates of the closest point inside the box.
    *
    * @param x
    */
  def clip( x: Vec ): Unit =
  {
    assert( nVars == x.length )
    var i = nVars
    while( i > 0 )
    {
      i -= 1
      x(i) = x(i)  max  x_min(i)  min  x_max(i)
    }
  }

  /** For a search direction vector, this method set all entries to zero
    * that would lead out of the bounding box.
    *
    * @param x
    * @param dir
    * @param sign
    */
  def clipDir( x: Vec, dir: Vec, sign: Int = +1 ): Unit =
  {
    assert( nVars == x.length )
    assert( nVars == dir.length )
    var i = nVars
    while( i > 0 )
    {
      i -= 1
      if( sign*dir(i) > 0 && x(i) == x_max(i) ) dir(i) = 0
      if( sign*dir(i) < 0 && x(i) == x_min(i) ) dir(i) = 0
    }
  }

  def contains( x: Vec ): Boolean =
  {
    @tailrec def contains( i: Int ): Boolean
      = i < 0 ||
        x(i) >= x_min(i) &&
        x(i) <= x_max(i) &&
        contains(i-1)
    contains(nVars-1)
  }
}
object BoxConstraints
{
  /** Creates a new BoxConstraints object from the given lower and upper bounds as well as the iteration
    * starting point. Tries to cleverly fill the gaps in case x_min, x_max and/or x_init are not specified.
    * At least one of the three however has to be specified.
    *
    * @param x_init
    * @param x_min
    * @param x_max
    * @return
    */
  def apply( x_init: Vec, x_min: Vec, x_max: Vec): (Vec,BoxConstraints) =
  {
    if( null == x_init.vals ) throw new NullPointerException
    if( null == x_min .vals ) throw new NullPointerException
    if( null == x_max .vals ) throw new NullPointerException

    val nVars = if( isDefined(x_init) ) x_init.length
           else if( isDefined(x_min ) ) x_min .length
           else if( isDefined(x_max ) ) x_max .length
           else throw new IllegalArgumentException("All 'param_'-arguments undefined.")

    val xMin = if( isUndefined(x_min) ) Vec.full(nVars)(- ∞)
          else if(x_min.length == nVars) x_min
          else throw new IllegalArgumentException("Not all 'param_'-arguments of same length.")

    val xMax = if( isUndefined(x_max) ) Vec.full(nVars)(+ ∞)
          else if(x_max.length == nVars) x_max
          else throw new IllegalArgumentException("Not all 'param_'-arguments of same length.")

    val x = if( isUndefined(x_init) )
        Vec.tabulate(nVars){
          i =>
            if( ! xMin(i).isInfinite )
              if( ! xMax(i).isInfinite )
                0.5 * xMin(i) + 0.5 * xMax(i)
              else
                xMin(i)
            else if( ! xMax(i).isInfinite )
              xMax(i)
            else
              throw new IllegalArgumentException("param_init unspecified and both param_min and param_max infinite.")
        }
      else if(x_init.length == nVars) x_init.clone
      else throw new IllegalArgumentException("Not all 'param_'-arguments of same length.")

    var i = nVars
    while( i > 0 ) {
      i -= 1
      if( ! {x(i) <= xMax(i)} ) throw new IllegalArgumentException("param_init > param_max or NaN.")
      else
      if( ! {x(i) >= xMin(i)} ) throw new IllegalArgumentException("param_init < param_min or NaN.")
    }

    ( x, new BoxConstraints(xMin,xMax) )
  }

  def apply( x_min: Vec, x_max: Vec ): BoxConstraints =
  {
    if( null == x_min.vals ) throw new NullPointerException
    if( null == x_max.vals ) throw new NullPointerException

    val nVars = if( isDefined(x_min ) ) x_min .length
           else if( isDefined(x_max ) ) x_max .length
           else throw new IllegalArgumentException("All 'x_'-arguments undefined.")

    val xMin = if( isUndefined(x_min) ) Vec.full(nVars)(- ∞)
          else if(x_min.length == nVars) x_min
          else throw new IllegalArgumentException("Not all 'x_'-arguments of same length.")

    val xMax = if( isUndefined(x_max) ) Vec.full(nVars)(+ ∞)
          else if(x_max.length == nVars) x_max
          else throw new IllegalArgumentException("Not all 'x_'-arguments of same length.")

    var i = nVars
    while( i > 0 ) {
      i -= 1
      if( ! {xMin(i) <= xMax(i)} )
        throw new IllegalArgumentException("x_min > x_max or NaN.")
    }

    new BoxConstraints(xMin,xMax)
  }
}