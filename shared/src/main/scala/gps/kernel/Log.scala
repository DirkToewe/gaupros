package gps.kernel

/**
  * Created by Dirk Toewe on 31.07.17.
  */
private[kernel] class Log[@specialized -X] private(val k: Kernel[X] ) extends Func[X](k)
{
  override def apply( xi: X, xj: X, i: Int, j: Int )
    = math.log{ k(xi,xj, i,j) }

  override def pDiff(sym: Symbol ): Kernel[X]
    = (k pDiff sym) / k

  override def subs[Y <: X](map: PartialFunction[Symbol,Kernel[Y]] ): Kernel[Y]
    = Log(k subs map)
}
object Log
{
  def apply[X]( kernel: Kernel[X] ): Kernel[X]
    = kernel match {
        case Exp(a) => a
        case Const(e) => math.log(e)
        case Pow(b,Const(1)) => new Log(b)
        case Pow(b,e) => e*Log(b)
        case _ => new Log(kernel)
      }

  def unapply[X]( log: Log[X] ) = Some(log.k)
}