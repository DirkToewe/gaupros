package gps.kernel

/**
  * Created by Dirk Toewe on 05.08.17.
  */
protected[kernel] abstract class Func[-X](arg: Kernel[X] ) extends Kernel[X]
{
  override def apply( xi: X, xj: X )
    = this(xi,xj, -1,-1)

  override def apply( xi: X, xj: X, i: Int, j: Int ): Double

  override protected[kernel] final def toString(xi: String, xj: String )
    = s"${getClass.getSimpleName}(${arg.toString(xi,xj)})"

  override protected[kernel] final val argDependent = arg.argDependent
  override final val params = arg.params
}
