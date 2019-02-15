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

package gps.opt

import gps.diff
import gps.linalg.Vec

/**
  * Created by Dirk Toewe on 09.08.17.
  */
trait ObjectiveWithGradient extends ObjectiveFunction
{
  def gradient( x: Vec ): Vec
    = diff.gradient(this)(x)

  /** Returns the function value and the gradient at in one call. This method can
    * be overwritten by implementations to improve performance.
    *
    * @param x
    * @return
    */
  def fval_grad( x: Vec )
    = ( this(x), gradient(x) )

  override def unary_- = new ObjectiveWithGradient
  {
    override def apply   ( x: Vec ) = - ObjectiveWithGradient.this(x)
    override def gradient( x: Vec ) = - ObjectiveWithGradient.this.gradient(x)

    override def fval_grad( x: Vec )
      = ObjectiveWithGradient.this.fval_grad(x) match { case (f,g) => (-f,-g) }
  }
}