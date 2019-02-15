/* This file is part of GauProS.
 *
 * GauProS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GauProS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GauProS.  If not, see <https://www.gnu.org/licenses/>.
 */

package gps.regression.gpr

import gps.kernel.Kernel
import gps.linalg.Vec

import scala.collection.{IndexedSeq => ISeq}

private[gpr] abstract class LikelihoodFunctionAbstract[@specialized -X](
  x: Array[X],
  y: Vec,
  kernel: Kernel[X],
  override val params: ISeq[Symbol]
) extends LikelihoodFunction
{
  if( x.length != y.length ) throw new IllegalArgumentException

  protected val nSamples = x.length

  protected val grad_cache:      ISeq[Kernel[X]]  = Array.tabulate(params.length){ kernel pDiff params(_) }
  protected val hess_cache: ISeq[ISeq[Kernel[X]]] = Array.tabulate(params.length){
    i =>
      val dkdp = kernel pDiff params(i)
      Array.tabulate(i+1){ dkdp pDiff params(_) }: ISeq[Kernel[X]]
  }

  protected def cov( param: Vec ): Kernel[X] =
  {
    if( params.length != param.length )
      throw new IllegalArgumentException(s"${params.length} parameters expected but ${param.length} given.")
    val sym2val = Map( params zip param map { case (k,v) => k -> (v: Kernel[X]) } :_* )
    kernel subs sym2val
  }

  protected def cov_grad( param: Vec ): (Kernel[X], Int => Kernel[X]) = {
    if( params.length != param.length )
      throw new IllegalArgumentException(s"${params.length} parameters expected but ${param.length} given.")
    val sym2val = Map( params zip param map { case (k,v) => k -> (v: Kernel[X]) } :_* )

    (
      kernel subs sym2val,
      grad_cache map (_ subs sym2val)
    )
  }

  protected def cov_grad_hess( param: Vec ): (Kernel[X], Int => Kernel[X], (Int,Int) => Kernel[X]) =
  {
    if( params.length != param.length )
      throw new IllegalArgumentException(s"${params.length} parameters expected but ${param.length} given.")
    val sym2val = Map( params zip param map { case (k,v) => k -> (v: Kernel[X]) } :_* )

    (
      kernel subs sym2val,
      grad_cache map (_ subs sym2val), // <- cache 1st order derivatives as they may be used multiple
      hess_cache(_)(_) subs sym2val    //    times in the 2nd derivatives of a likelihood function
    )
  }
}
