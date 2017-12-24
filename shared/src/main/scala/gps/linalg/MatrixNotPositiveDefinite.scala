package gps.linalg

/** Thrown by methods when a matrix argument is not positive definite (numerically speaking).
  *
  * Created by Dirk Toewe on 28.08.17.
  */
case class MatrixNotPositiveDefinite() extends IllegalArgumentException("Matrix not positive definite.") {}
