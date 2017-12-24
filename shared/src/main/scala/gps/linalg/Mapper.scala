package gps.linalg

/** A functional interface intended for mapping matrix entries.
  *
  * Created by Dirk Toewe on 15/08/17.
  */
trait Mapper[@specialized +M]
{
  @inline def apply( A_ij: Double, i: Int, j: Int ): M
}