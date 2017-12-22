package gps.opt

import Double.{ PositiveInfinity => ∞ }
import gps.linalg._

import scala.annotation.tailrec

/**
  *
  * Created by Dirk Toewe on 25.08.17.
  *
  * @see Numerical Optimization, (Jorge Nocedal, Stephen J. Wright), Springer, 1999, Algorithm 3.2 & Algorithm 3.3.
  * @param c1 Constant of the <i>sufficient decrease condition</i> (1st strong Wolfe condition).
  *           Has to be greater than 0 and less than c2.
  * @param c2 Constant of the <i>curvature condition</i> (2nd strong Wolfe condition).
  *           Has to be greater than c1 and less than 1.
  * @param c3 Exponential growth factor. Whenever the line search range <b>(α<sub>i-1</sub>, α<sub>i</sub>]</b> is not guaranteed
  *           to contain strong Wolfe points, <b>(α<sub>i</sub>,α<sub>i</sub>*c3]</b> is searched instead.
  */
case class LineSearchStrongWolfe( objective: ObjectiveWithGradient, c1: Double=0.3, c2: Double=0.8, c3: Double=1.5 ) extends LineSearchMethod
{
  assert(c1 >  0)
  assert(c1 < c2)
  assert( 1 > c2)
  assert( 1 < c3)

  override def apply( dir: Vec, x_0: Vec, φ_0: Double, g_0: Vec,
                      α_init: Double = 1, α_max: Double = ∞ ): (Vec,Double,Vec) =
  {
    val dφ_0 = g_0⋅dir
    assert( dφ_0 < -0.0 )

    @tailrec def zoom( α_a: Double, φ_a: Double, α_b: Double, φ_b: Double ): (Vec,Double,Vec) =
    {
      assert( α_a >= +0.0 )
      assert( α_b >= +0.0 )
      val α = (α_a + α_b)/2
      val X = x_0 + α*dir
      val (φ,g) = objective.fval_grad(X)

      if( ! ( α != α_a && α != α_b ) ) // <- handles numeric issues, where the result between
        (X,φ,g)                        //    two values that are too close numerically

      else if( φ - φ_0 > c1 * dφ_0 * α || φ > φ_a )
        zoom(α_a,φ_a, α,φ)

      else {
        val dφ = g⋅dir
        if( dφ.abs <= -c2 * dφ_0 )
          (X,φ,g)
        else if( dφ * (α_b - α_a) >= 0 )
          zoom(α,φ, α_a,φ_a)
        else
          zoom(α,φ, α_b,φ_b)
      }
    }

    /** This method finds a region that is guaranteed to contain valid values according to strong Wolfe.
      */
    @tailrec def bracket( α_min: Double, φ_min: Double, α: Double ): (Vec,Double,Vec) =
    {
      val X = x_0 + α*dir
      val (φ,g) = objective.fval_grad(X)
      if( φ - φ_0 > c1 * dφ_0 * α  ||  α_min > +0.0 && φ > φ_min )
        zoom(α_min,φ_min, α,φ)
      else {
        val dφ = g⋅dir
        if( dφ.abs <= -c2 * dφ_0 )
          (X,φ,g)
        else if( dφ >= 0 )
          zoom(α,φ, α_min,φ_min)
        else if( α < α_max )
          bracket(α,φ, α*c3 min α_max)
        else
          throw new LineSearchFailure(α_max, "Strong wolfe conditions not satisfiable in α-range.")
      }
    }

    bracket(+0.0, φ_0, α_init)
  }
}