package gps.kernel

/**
  * Created by dtitx on 30.07.17.
  */
protected[kernel] class Const private(val value: Double ) extends Kernel[Any]
{
  override def apply( xi: Any, xj: Any )
    = value

  override protected[kernel] def toString(xi: String, xj: String ) = value.toString
  override protected[kernel] def argDependent = false
}
object Const
{
  def apply( value: Double ): Kernel[Any]
    = new Const(value)

  def unapply( const: Const ): Option[Double]
    = Some(const.value)
}