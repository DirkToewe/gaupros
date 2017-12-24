package gps.opt.test_functions

import gps.opt.ObjectiveWithHessian_test

/**
  * Created by Dirk Toewe on 26.08.17.
  */
object RosenbrockFunction2D_test extends ObjectiveWithHessian_test(RosenbrockFunction, Vec(-2,-2), Vec(+2,+2), 1e-4, 1e-4 )
object RosenbrockFunction3D_test extends ObjectiveWithHessian_test(RosenbrockFunction, Vec(-2,-2,-2), Vec(+2,+2,+2), 1e-4, 1e-4 )
object RosenbrockFunction4D_test extends ObjectiveWithHessian_test(RosenbrockFunction, Vec(-2,-2,-2,-2), Vec(+2,+2,+2,+2), 1e-4, 1e-4 )
