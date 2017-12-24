package gps.opt

import gps.linalg.Vec

/** Exception thrown by iterative optimizers and solvers if the maximum number of iterations
  * is exceeded without achieving the demanded tolerance in the solution.
  *
  * @param x The iteration result which does not meet the tolerance requirements.
  *
  * Created by Dirk Toewe on 23.08.17.
  */
case class NoConvergence( x: Vec ) extends RuntimeException {}