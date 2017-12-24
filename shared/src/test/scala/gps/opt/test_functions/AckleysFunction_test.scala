package gps.opt.test_functions

import gps.opt.ObjectiveWithHessian_test

/**
  * Created by Dirk Toewe on 26.08.17.
  */
object AckleysFunction1D_test extends ObjectiveWithHessian_test(AckleysFunction, Vec(-32), Vec(+32), 1e-4, 1e-4 )
object AckleysFunction2D_test extends ObjectiveWithHessian_test(AckleysFunction, Vec(-32,-32), Vec(+32,+32), 1e-3, 1e-3 )
object AckleysFunction3D_test extends ObjectiveWithHessian_test(AckleysFunction, Vec(-32,-32,-32), Vec(+32,+32,+32), 1e-3, 1e-3 )
object AckleysFunction4D_test extends ObjectiveWithHessian_test(AckleysFunction, Vec(-32,-32,-32,-32), Vec(+32,+32,+32,+32), 1e-4, 1e-4 )
