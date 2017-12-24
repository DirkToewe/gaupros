package gps.kernel

import scala.collection.immutable.IndexedSeq

/**
  * Created by dtitx on 30.07.17.
  */
private[kernel] class Sym private(val sym: Symbol ) extends Kernel[Any]
{
  override def apply( xi: Any, xj: Any )
    = throw new IllegalStateException(s"Parameter '$sym' not specified.")

  override def pDiff( symbol: Symbol )
    = if( symbol == sym ) 1 else 0

  override def subs[Y <: Any](map: PartialFunction[Symbol,Kernel[Y]] ): Kernel[Y]
    = map applyOrElse (sym, (_: Symbol) => this)

  override protected[kernel] def toString(xi: String, xj: String ) = sym.name
  override protected[kernel] def argDependent = false
  override val params = IndexedSeq(sym)
}
object Sym
{
  def apply( name: Symbol )
    = new Sym(name)

  def unapply( sym: Sym )
    = Some(sym.sym)
}
