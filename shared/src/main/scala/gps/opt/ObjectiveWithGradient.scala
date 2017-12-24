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