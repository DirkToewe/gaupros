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

package gps.opt.test_functions

import gps.linalg.{ LMat, Vec }
import gps.opt.ObjectiveWithHessian

import scala.math.{cos, sin, Pi => π}

/** The n-dimensional <a href="https://en.wikipedia.org/wiki/Rastrigin_function">Rastrigin function</a>.
  * May be called with arguments of any dimensions.
  *
  * Created by Dirk Toewe on 26.08.17.
  */
object RastriginFunction extends ObjectiveWithHessian
{
  val A = 10.0

  override def apply( x: Vec ) =
  {
    val n = x.length
    A*n + (0 until n).map( i => x(i)*x(i) - A*cos{2*π*x(i)} ).sum
  }

  override def gradient( x: Vec )
    = Vec.tabulate(x.length){
        i => 2*x(i) + 2*π*A * sin( 2*π * x{i} )
      }

  override def hessian(x: Vec)
    = LMat.tabulate(x.length){
        (i,j) => if( i != j ) 0 else 2 + 4*π*π*A * cos( 2*π * x{i} )
      }
}
