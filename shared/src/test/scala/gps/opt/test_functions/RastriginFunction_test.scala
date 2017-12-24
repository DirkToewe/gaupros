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