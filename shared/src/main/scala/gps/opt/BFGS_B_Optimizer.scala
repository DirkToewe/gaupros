package gps.opt

import gps.linalg._

import scala.Double.{PositiveInfinity => ∞}
import scala.annotation.tailrec
import scala.math.{sqrt => √}

/** A simple implementation of the <a href="https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm">BFGS algorithm</a>.
  *
  * Created by Dirk Toewe on 22.08.17.
  *
  * @param tolRel Relative tolerance for the gradient to be considered zero.
  *               Stopping condition is: <b>‖∇f‖₂ / √n ≤ tolAbs + tolRel * |f|</b>
  * @param tolAbs Absolute tolerance for the gradient to be considered zero.
  * @param α_init0 The first iteration of BFGS is a gradient descent. α_init0 determines the initial guess for
  *                the line search along the gradient. Reduce if the optimizer overshoots in the first iteration.
  * @param maxIter Maximum allowed number of iterations. If exceeded, `NoConvergence` is thrown.
  * @param lineSearchMethod Line search method to be used. Should satisfy Wolfe or strong Wolfe conditions
  *                         in order to keep the Hessian approximation positive definite after each update.
  */
case class BFGS_B_Optimizer(
  tolRel: Double = 1e-6, tolAbs: Double = 1e-8, α_init0: Double = 1e-3, maxIter: Long = 16*1024,
  lineSearchMethod: ObjectiveWithGradient => LineSearchMethod = LineSearchStrongWolfe(_)
) extends Optimizer_B[ObjectiveWithGradient]
{
  assert( tolRel >= 0 )
  assert( tolAbs >= 0 )
  assert( maxIter >= 0 )
  assert( null != lineSearchMethod )

  override def maximize( objective: ObjectiveWithGradient, x_init: Vec ): Vec
    = minimize(-objective, x_init )

  override def minimize( objective: ObjectiveWithGradient, x_init: Vec ): Vec =
  {
    if( null == objective  ) throw new NullPointerException
    if( null == x_init.vals) throw new NullPointerException

    val nVars = x_init.length

    if( 0 == nVars )
      return Vec.zeros(0)

    // inverse of current hessian approx. Initializing it to I means Gradient Descent on 1st iteration
    val H = LMat.tabulate(nVars){ (i, j) => if (i == j) α_init0 else 0 }

    val lineSearch = lineSearchMethod(objective)

    @tailrec def iterate( iter: Long, x: Vec, f: Double, g: Vec ): Vec =
    {
      if( iter >= maxIter )
        throw new NoConvergence(x)

      if( √(normPow(g) / nVars) <= tolAbs + tolRel * f.abs )
        return x

      // determine search direction from inverse Hessian
      val dir = H symMatMul g
      dir *= -1
      val (_x,_f,_g) = lineSearch(dir, x,f,g, 1,∞)

      // update inverse of Hessian approximation
      val dx = _x - x // <- often referred to as 's'
      val dg = _g - g // <- often referred to as 'y'
      var dxdg = dx⋅dg; assert( dxdg > +0.0, s"dxdg = ${dxdg} <= +0.0" )
      val Hdg = H symMatMul dg

      dx /= dxdg; dxdg += dg⋅Hdg

      H modify {
        (Bij,i,j) =>
          Bij - Hdg(i)*dx(j) - Hdg(j)*dx(i) + dxdg * dx(i)*dx(j)
      }

      iterate(iter+1,_x,_f,_g)
    }

    // start iteration
    val (f,g) = objective.fval_grad(x_init)
    iterate(0, x_init,f,g)
  }

  override def maximize( objective: ObjectiveWithGradient, x_init: Vec=UNDEFINED_ARRAY, x_min: Vec=UNDEFINED_ARRAY, x_max: Vec=UNDEFINED_ARRAY ) =
  {
    val (x0,box) = BoxConstraints(x_init,x_min,x_max)
    val obj = UnboxedObjective(objective, box.x_min, box.x_max)
    val y0 = obj transform x0
    val y = maximize(obj,y0)
    obj transform_back y
  }

  override def minimize( objective: ObjectiveWithGradient, x_init: Vec=UNDEFINED_ARRAY, x_min: Vec=UNDEFINED_ARRAY, x_max: Vec=UNDEFINED_ARRAY ) =
  {
    val (x0,box) = BoxConstraints(x_init,x_min,x_max)
    val obj = UnboxedObjective(objective, box.x_min, box.x_max)
    val y0 = obj transform x0
    val y = minimize(obj,y0)
    obj transform_back y
  }
}