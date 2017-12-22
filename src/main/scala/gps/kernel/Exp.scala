package gps.kernel

/**
  * Created by Dirk Toewe on 30.07.17.
  */
private[kernel] class Exp[-X] private(val e: Kernel[X] ) extends Func[X](e)
{
  override def apply( xi: X, xj: X, i: Int, j: Int )
    = math.exp{ e(xi,xj, i,j) }

  override def pDiff( sym: Symbol ): Kernel[X]
    = this * (e pDiff sym)

  override def subs[Y <: X](map: PartialFunction[Symbol,Kernel[Y]] ): Kernel[Y]
    = Exp(e subs map)
}
object Exp
{
  def apply[X]( exponent: Kernel[X] ): Kernel[X]
    = exponent match {
        case Log(kernel) => kernel
        case Const(e) => math.exp(e)
        case kernel => new Exp(kernel)
      }

  def unapply[X]( exp: Exp[X] ) = Some(exp.e)
}