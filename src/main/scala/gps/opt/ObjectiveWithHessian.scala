package gps.opt

import gps.linalg._

/**
  * Created by Dirk Toewe on 09.08.17.
  */
trait ObjectiveWithHessian extends ObjectiveWithGradient
{
  /** Returns the hessian matrix at input x. The hessian matrix is assumed to be symmetric.
    * Which means the second (partial) derivatives must be symmetric around x, which is
    * usually the case but not always (see <a href="https://en.wikipedia.org/wiki/Symmetry_of_second_derivatives">Wikipedia</a>).
    *
    * @param x
    * @return
    */
  def hessian( x: Vec ): LMat

  override def gradient( x: Vec ): Vec

  /** Returns the function value, the gradient and the hessian matrix in one call.
    * This method can be overwritten by implementations to improve performance.
    *
    * @param x
    * @return
    */
  def fval_grad_hess( x: Vec )
    = ( this(x), gradient(x), hessian(x) )

  override def unary_- = new ObjectiveWithHessian
  {
    override def apply   ( x: Vec ) = - ObjectiveWithHessian.this(x)
    override def gradient( x: Vec ) = - ObjectiveWithHessian.this.gradient(x)
    override def hessian ( x: Vec ) = - ObjectiveWithHessian.this.hessian(x)

    override def fval_grad( x: Vec )
      = ObjectiveWithHessian.this.fval_grad(x) match { case (f,g) => (-f,-g) }

    override def fval_grad_hess( x: Vec )
      = ObjectiveWithHessian.this.fval_grad_hess(x) match { case (f,g,h) => (-f,-g,-h) }
  }
}