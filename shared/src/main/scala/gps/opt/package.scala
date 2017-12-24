package gps

import gps.linalg.Vec

/** Rudimentary optimization methods and utilities.
  *
  * Created by Dirk Toewe on 20.06.17.
  */
package object opt
{
  private[opt] val UNDEFINED_ARRAY = Vec.zeros(0)

  private[opt] def isDefined  ( vec: Vec ) = UNDEFINED_ARRAY.vals ne vec.vals
  private[opt] def isUndefined( vec: Vec ) = UNDEFINED_ARRAY.vals eq vec.vals

  def regulaFalsi( f: Double => Double, lo: Double, hi: Double )
    = ???
}