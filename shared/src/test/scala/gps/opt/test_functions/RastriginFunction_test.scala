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

import gps.linalg.Vec
import gps.opt.ObjectiveWithHessian_test

/**
  * Created by Dirk Toewe on 26.08.17.
  */
object RastriginFunction1D_test extends ObjectiveWithHessian_test(RastriginFunction, Vec(-5), Vec(+5), 1e-3, 1e-3 )
object RastriginFunction2D_test extends ObjectiveWithHessian_test(RastriginFunction, Vec(-5,-5), Vec(+5,+5), 1e-3, 1e-3 )
object RastriginFunction3D_test extends ObjectiveWithHessian_test(RastriginFunction, Vec(-5,-5,-5), Vec(+5,+5,+5), 1e-3, 1e-3 )
object RastriginFunction4D_test extends ObjectiveWithHessian_test(RastriginFunction, Vec(-5,-5,-5,-5), Vec(+5,+5,+5,+5), 1e-3, 1e-3 )