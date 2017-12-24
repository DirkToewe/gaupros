package gps.kernel

/** A Kernel[Vec] that uses a nested Kernel[Double] to zip the entries of the vectors
  * arguments pairwise and sums them up afterwards.
  *
  * Created by Dirk Toewe on 30.07.17.
  */
private[kernel] class VecSum( val mapOp: Kernel[Double] ) extends Kernel[Vec]
{
  override def apply( xi: Vec, xj: Vec )
  = this(xi,xj, -1,-1)

  override def apply( xi: Vec, xj: Vec, i: Int, j: Int ) =
  {
    assert(xi.length > 0)
    assert(xi.length == xj.length)

    def distPow( k: Int, len: Int ): Double
      = if( len == 1 )
          mapOp(xi{k},xj{k}, i,j)
        else {
          val l = len / 2
          distPow(k,l) +
            distPow(k+l, len-l)
        }

    distPow(0, xi.length)
  }

  override def subs[Y <: Vec]( map: PartialFunction[Symbol,Kernel[Y]] ): Kernel[Vec] =
  {
    val _map: PartialFunction[Symbol,Kernel[Double]] = {
      case sym =>
        val result = map applyOrElse (sym, (_: Any) => sym: Kernel[Double] )
        if( result.argDependent )
          throw new IllegalArgumentException( s"Cannot substitute '$sym' by '$result' since it is part of a nested kernel." )
        result.asInstanceOf[Kernel[Double]]
    }
    VecSum(mapOp subs _map)
  }

  override def pDiff(sym: Symbol )
    = VecSum(mapOp pDiff sym)

  override val params = mapOp.params

  override protected[kernel] def toString(xi: String, xj: String ) =
  {
    val Xi = xi.substring(0, xi.length-1) + ",k]"  ensuring  (xi endsWith "]")
    val Xj = xj.substring(0, xj.length-1) + ",k]"  ensuring  (xj endsWith "]")
    s"Σₖ{ ${mapOp.toString(Xi,Xj)} }"
  }
}
object VecSum
{
  def apply( mapOp: Kernel[Double] ): Kernel[Vec]
    = mapOp match {
        case Independent(a) * b => a * VecSum(b)
        case _ => new VecSum(mapOp)
      }

  def unapply( argSum: VecSum )
    = Some(argSum.mapOp)
}