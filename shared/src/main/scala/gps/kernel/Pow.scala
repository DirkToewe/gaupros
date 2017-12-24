package gps.kernel

private[kernel] class Pow[-X] private(val b: Kernel[X], val e: Kernel[X] ) extends BinOp[X](b,e)
{
  override def apply( xi: X, xj: X, i: Int, j: Int ): Double
    = math.pow(
        b(xi,xj, i,j),
        e(xi,xj, i,j)
      )

  protected override def eval[Y]( b: Kernel[Y], e: Kernel[Y] ): Kernel[Y]
    = Pow(b,e)

  protected override def diff[Y <: X]( db: Kernel[Y], de: Kernel[Y] ): Kernel[Y]
    = db * Pow(b,e-1)*e  +  Pow(b,e)*Log(b) * de

  override protected[kernel] def toString(xi: String, xj: String ) =
  {
    var bStr = b.toString(xi,xj)
    bStr = b match {
      case Const(_) | Sym(_) | Exp(_) | Log(_) | _: AbsDelta.type => bStr
      case _ => s"($bStr)"
    }
    var eStr = e.toString(xi,xj)
    eStr = e match {
      case Const(2) => "²"
      case Const(3) => "³"
      case Const(_) | Sym(_) | Exp(_) | Log(_) | _: AbsDelta.type => s"^$eStr"
      case _ => s"^($eStr)"
    }
    bStr + eStr
  }
}
private[kernel] object Pow
{
  def apply[X]( base: Kernel[X], exponent: Kernel[X] ): Kernel[X]
    = (base,exponent) match {
        case (      _,  Const(0)) => Const(1)
        case (      _,  Const(1)) => base
        case (Const(b), Const(e)) => math.pow(b,e)
        case _ => new Pow(base,exponent)
      }

  def unapply[X]( pow: Pow[X] )
    = Some(pow.b, pow.e)
}
private[kernel] object Pow1
{
  /** Always matches power, inserting a 1 as exponent if need be.
    *
    * @param kernel
    * @tparam X
    * @return
    */
  def unapply[X]( kernel: Kernel[X] )
    = kernel match {
        case Pow(b,e) => Some(b,e)
        case _ => Some(kernel,Const(1))
      }
}