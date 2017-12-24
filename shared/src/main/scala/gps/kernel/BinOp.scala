package gps.kernel

import gps.util.Seqs.binSearch

import scala.collection.immutable.IndexedSeq

/** Super trait of Kernels that are created from two Kernels using a binary
  * operator like add, multiply, power, ...
  *
  * Created by Dirk Toewe on 31.07.17.
  */
private[kernel] abstract class BinOp[-X] protected(a: Kernel[X], b: Kernel[X] ) extends Kernel[X]
{
  /** Creates a new kernel applying the same operator that was used in this
    * Kernel.
    *
    * @param a
    * @param b
    * @tparam Y
    * @return
    */
  protected def eval[Y]( a: Kernel[Y], b: Kernel[Y] ): Kernel[Y]

  /** Calculates the derivative of this kernel given the derivatives of the two
    * operator arguments/kernels.
    *
    * @param da
    * @param db
    * @tparam Y
    * @return
    */
  protected def diff[Y <: X]( da: Kernel[Y], db: Kernel[Y] ): Kernel[Y]

  override def pDiff(sym: Symbol )
    = if( 0 > binSearch(params,sym) )
        0
      else diff(
        a pDiff sym,
        b pDiff sym
      )

  override def apply( xi: X, xj: X )
    = this(xi,xj, -1,-1)

  override def subs[Y <: X](map: PartialFunction[Symbol,Kernel[Y]] ): Kernel[Y]
    = eval(
        a subs map,
        b subs map
      )

  override protected[kernel] val argDependent
    = a.argDependent || b.argDependent

  override val params = {
    val arr = sortedUnion(a.params, b.params).toArray
    new IndexedSeq[Symbol]{
      override def apply( idx: Int ) = arr(idx)
      override def length = arr.length
    }
  }
}