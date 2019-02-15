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

package gps

import gps.linalg.{LMat, Mat, _}

/**
  * Created by dtult on 16/07/17.
  */
package object diff
{
  /** Takes a ℝ<sup>n</sup> → ℝ function <b>f</b> and returns a ℝ<sup>n</sup> → ℝ<sup>n</sup>
    * function which numerically computes the gradient <b>∇f</b> using the central difference
    * method.
    *
    * The finite difference step is calculated as follows: dx<sub>i</sub> = dxAbs + dxRel*|X<sub>i</sub>|.
    *
    * @param f The function which is to be differentiated numerically.
    * @return <b>∇f</b>, The gradient of f.
    */
  def gradient( f: Vec => Double, dxAbs: Double = 1e-7, dxRel: Double = 1e-7 ): Vec => Vec
    = x => {
        val X = x.clone
        Vec.tabulate(X.length){
          i =>
            val X_i = X(i)
            val dx = dxAbs + dxRel * X_i.abs
            val dX = (X_i + dx) - (X_i - dx)
            X(i) = X_i + dx; val fHi = f(X)
            X(i) = X_i - dx; val fLo = f(X)
            X(i) = X_i
            (fHi - fLo) / dX
        }
      }

  def jacobian( g: Vec => Vec, dxAbs: Double = 1e-7, dxRel: Double = 1e-7 ): Vec => Mat
    = x => {
        val X = x.clone
        MatCM.tabulateCM(g(X).length,X.length){
          i =>
            val X_i = X(i)
            val dx = dxAbs + dxRel * X_i.abs
            X(i) = X_i + dx; val dGdX  = g(X)
            X(i) = X_i - dx;     dGdX -= g(X)
            X(i) = X_i
            dGdX /= (X_i + dx) - (X_i - dx)
            dGdX
        }
      }

  def hessian( f: Vec => Double ): Vec => LMat
    = x => {
        val jac = jacobian( gradient(f) )(x)
        jac foreach {
          (jac_ij,i,j) => assert( isClose(jac(j,i), jac_ij) )
        }
        LMat.tabulate(x.length)( (i,j) => {jac(i,j) + jac(j,i)} / 2 )
      }
}
