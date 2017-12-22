package gps.kernel

/**
  * Created by Dirk Toewe on 01.08.17.
  */
object AbsDelta extends Kernel[Double]
{
  override def apply( xi: Double, xj: Double )
    = (xi - xj).abs

  override protected[kernel] def toString(xi: String, xj: String ) = s"|$xi - $xj|"
}
