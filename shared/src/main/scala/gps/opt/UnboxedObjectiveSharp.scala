package gps.opt

import gps.linalg._

import scala.annotation.tailrec

/**
  *
  */
object UnboxedObjectiveSharp
{
  def apply( objective: ObjectiveFunction, x_min: Vec, x_max: Vec, margin: Double ): UnboxedObjectiveSharp
    = objective match {
        case obj: ObjectiveWithGradient => this(obj,x_min,x_max,margin)
        case _ =>
          val box = BoxConstraints(x_min,x_max)
          new UnboxedObjectiveSharp(objective, box.x_min, box.x_max,margin)
      }

  def apply( objective: ObjectiveWithGradient, x_min: Vec, x_max: Vec, margin: Double ): UnboxedObjectiveSharp with ObjectiveWithGradient
    = objective match {
        case obj: ObjectiveWithHessian => this(obj,x_min,x_max,margin)
        case _ =>
          val box = BoxConstraints(x_min,x_max)
          new UnboxedObjectiveSharpWithGradient(objective, box.x_min, box.x_max,margin)
      }

  def apply( objective: ObjectiveWithHessian, x_min: Vec, x_max: Vec, margin: Double ): UnboxedObjectiveSharp with ObjectiveWithHessian =
  {
    val box = BoxConstraints(x_min,x_max)
    new UnboxedObjectiveSharpWithHessian(objective, box.x_min, box.x_max,margin)
  }

  private[opt] trait Trafo
  {
    def apply  ( x: Double ): Double
    def inverse( y: Double ): Double
    def deriv1 ( x: Double ): Double
    def deriv2 ( x: Double ): Double
  }


  private[opt] object TrafoIdentity extends Trafo
  {
    override def apply  ( x: Double ) = x
    override def inverse( y: Double ) = y
    override def deriv1 ( x: Double ) = 1
    override def deriv2 ( x: Double ) = 0
  }


  private[opt] class Trafo1Bound( min: Double, private var fm: Double ) extends Trafo
  {
    // Looking for Polynomial f(x) where:
    //
    // f(0) = 0
    // f'(0) = 0
    // f''(m) = 0
    // f'(m) = s
    // f(m) = fm
    //
    // => f(x) = a*x² + b*x³
    // with:
    //   2a + 6*b*m = 0  =>  a = -3*b*m
    //   2a*m + 3*b*m² = s  =>  a*m = s => a = s/m
    //
    //   fm = s/m * m² - s/(3m²)*m³ = s*m - s*m/3 => fm = s*m/3 => m = 3*fm/s

    private val s = fm.signum
    private val m = 3*fm/s
    private val a = s / m
    private val b = -a / (3*m)
    fm = m*m * (a + b*m)

    override def apply( X: Double ) =
    {
      val x = X.abs
      min + { if( x < m ) x*x * (a + b*x) else fm + (x-m)*s }
    }

    override def deriv1( X: Double ): Double =
    {
      val x = X.abs
      if( x < m ) X * (2*a + 3*b*x) else s*X.signum
    }

    override def deriv2( X: Double ): Double =
    {
      val x = X.abs
      if( x < m ) 2*a + 6*b*x else 0
    }

    override def inverse( Y: Double ) =
    {
      val y = Y - min
      assert(y*s >= 0)

      // TODO implement and use opt.regulaFalsi instead
      @tailrec def secantSearch( α: Double, β: Double, f_α: Double, f_β: Double ): Double =
      {
        val mid = (α+β) / 2
        if( α != mid && mid != β ) // handles NaN and values in between float gaps
        {
          val x = (f_β*α - f_α*β) / (f_β - f_α)
          val f_x = x*x * (a + b*x)  -  y
          if( f_x == f_β )
            (β+x) / 2
          else
            secantSearch(β,x, f_β,f_x)
        }
        else if( f_α.abs <= f_β.abs ) α else β
      }

      val t = (y-fm) * s

      if( t < 0 )
        secantSearch(0,m, -y,fm-y)
      else
        m + t
    }
  }


  private[opt] class Trafo2Bounds( min: Double, max: Double, private var fm: Double ) extends Trafo
  {
    assert( fm >= 0 )
    assert( fm*2 <= max-min )
    assert( min  <= max )

    // looking for a polynomial f(x) where:
    // f(0) = 0
    // f'(0) = 0
    // f'(m) = ( max-min - 2*f(m) ) / (1 - 2m)
    // f''(m) = 0
    //
    //
    // => f(x) = a*x² + b*x³
    // with:
    //   2a + 6*b*m = 0  =>  a = -3*b*m  =>  b = -a / (3*m)
    //   2*a*m + 3*b*m² = a*m  =  (max-min - 2*a*m² - 2*b*m³) / (1 - 2m) = (max-min - 2*a*m² + 2/3*a*m²) / (1-2m)
    //
    //   => a*m - 2*a*m²  =  max-min - 2*a*m² + 2/3*a*m²
    //   => a*m  =  max-min + 2/3*a*m²
    //   => a*(m - 2/3*m²)  =  max-min
    //   => a = (max - min) * 3 / (3*m - 2*m²)
    //
    // m is given by:
    //   f(m) = fm = a*m² + b*m³ = a*m² - a/3 * m² = 2/3 * a * m²
    //
    //   => fm = 2 * (max-min) * m / (3 - 2*m)
    //   => fm * (3 - 2*m) = 2 * (max-min) * m
    //   => 2*m * (max-min + fm) = 3*fm
    //   => m = fm / (max-min + fm) * 3/2

    def h = max - min
    val m = fm / (h+fm) * 3/2
    val a = -3*h / (m * {2*m - 3})
    val b = -a / (3*m)
    fm = m*m * (a + b*m)

    override def apply( X: Double ) =
    {
      var x = X % 2.0
      if( x > +1 ) x -= 2
      if( x < -1 ) x += 2
      x = x.abs
      if( x <    m ) min + x*x * (a + b*x)
      else if( x <= 1-m ) min + fm + (x-m)/(1-2*m) * (h-2*fm)
      else{    x  = 1-x;  max - x*x * (a + b*x) }
    }

    override def deriv1( X: Double ) =
    {
      var x = X % 2.0
      if( x > +1 ) x -= 2
      if( x < -1 ) x += 2
      val s = x.signum
      x = x.abs
      s * {  if( x <=   m ) x * (2*a + 3*b*x)
      else if( x <  1-m ) (h-2*fm) / (1-2*m)
      else{    x  = 1-x;  x * (2*a + 3*b*x) }
      }
    }

    override def deriv2( X: Double ) =
    {
      var x = X % 2.0
      if( x > +1 ) x -= 2
      if( x < -1 ) x += 2
      x = x.abs
      if( x <=   m ) 2*a + 6*b*x
      else if( x <  1-m ) 0.0
      else{    x  = 1-x; -2*a - 6*b*x }
    }

    override def inverse( y: Double ) =
    {
      assert( y >= min )
      assert( y <= max )

      def inv( y: Double ) =
      {
        @tailrec def secantSearch( α: Double, β: Double, f_α: Double, f_β: Double ): Double =
        {
          val mid = (α+β) / 2
          if( α != mid && mid != β ) // handles NaN and values in between float gaps
          {
            val x = (f_β*α - f_α*β) / (f_β - f_α)
            val f_x = x*x * (a + b*x)  -  y
            if( f_x == f_β )
              (β+x) / 2
            else
              secantSearch(β,x, f_β,f_x)
          }
          else if( f_α.abs <= f_β.abs ) α else β
        }
        assert( y >= 0 )
        assert( y < fm )
        secantSearch(0,m, 0-y,fm-y)
      }

      val lo = min + fm
      val hi = max - fm
      if     ( y <  lo ) inv(y-min)
      else if( y <= hi ) m + (y-lo) / (h-2*fm) * (1-2*m)
      else           1 - inv(max-y)
    }
  }
}
/** A wrapper wrapper function f(y ∈ ℝ) around a boxed b(x ∈ [x_min,x_max]).
  * Each minimum and maximum of f(y') has a corresponding extremum on b( x(y') )
  * where `transform_back` is function x(y).
  *
  * This allows to solve box constrained problems using unconstrained optimizers.
  *
  * Created by Dirk Toewe on 29.08.17.
  */
class UnboxedObjectiveSharp protected[opt](unboxed: ObjectiveFunction, val x_min: Vec, val x_max: Vec, val margin: Double ) extends ObjectiveFunction
{
  import UnboxedObjectiveSharp._

  assert( x_min.length == x_max.length )
  @inline protected[opt] def n = x_min.length

  protected[opt] val trafos = Array.tabulate(n){
    i =>
      val xMin = x_min(i)
      val xMax = x_max(i)

      if( xMin.isInfinite )

        if( xMax.isInfinite ) TrafoIdentity
        else              new Trafo1Bound(xMax, -margin)

      else if( xMax.isInfinite ) new Trafo1Bound(xMin, +margin)
      else                       new Trafo2Bounds(xMin,xMax, margin)
  }

  /** Transforms an input of this (unboxed) function into the equivalent input for the wrapped (box constrained) function.
    *
    * This is used, for example, to transform the result of the optimization back into boxed coordinates.
    *
    * @param y
    * @return
    */
  def transform_back( y: Vec ): Vec
    = Vec.tabulate(n){ i => trafos(i){y(i)} }

  /** Transforms an input of the wrapped (box constrained) function into the equivalent input for this (unboxed) function.
    *
    * This is used, for example, to transform the stating point of the optimization into unboxed coordinates.
    *
    * @param x
    * @return
    */
  def transform( x: Vec ): Vec
    = Vec.tabulate(n){
        i =>
          assert{ x(i) >= x_min(i) }
          assert{ x(i) <= x_max(i) }
          trafos(i).inverse{x(i)}
      }

  override def apply( y: Vec ): Double
    = unboxed{ transform_back(y) }
}


private[opt] class UnboxedObjectiveSharpWithGradient private[opt](
  unboxed: ObjectiveWithGradient, x_min: Vec, x_max: Vec, margin: Double
) extends UnboxedObjectiveSharp(unboxed,x_min,x_max,margin) with ObjectiveWithGradient
{
  override def apply( y: Vec ) = super[UnboxedObjectiveSharp].apply(y)

  override def gradient( y: Vec ) =
  {
    val g = unboxed.gradient{ transform_back(y) }
    Vec.tabulate(n){ i => g(i) * trafos(i).deriv1(y{i}) }
  }

  override def fval_grad( y: Vec ) =
  {
    val (f,g) = unboxed.fval_grad{ transform_back(y) }
    val G = Vec.tabulate(n){ i => g(i) * trafos(i).deriv1(y{i}) }
    (f,G)
  }
}


private[opt] class UnboxedObjectiveSharpWithHessian private[opt](
  unboxed: ObjectiveWithHessian, x_min: Vec, x_max: Vec, margin: Double
) extends UnboxedObjectiveSharpWithGradient(unboxed,x_min,x_max,margin) with ObjectiveWithHessian
{
  override def apply    ( y: Vec ) = super[UnboxedObjectiveSharpWithGradient].apply(y)
  override def gradient ( y: Vec ) = super[UnboxedObjectiveSharpWithGradient].gradient(y)
  override def fval_grad( y: Vec ) = super[UnboxedObjectiveSharpWithGradient].fval_grad(y)

  override def hessian( y: Vec ) = fval_grad_hess(y)._3

  override def fval_grad_hess( y: Vec ): (Double, Vec, LMat) =
  {
    val (f,g,h) = unboxed.fval_grad_hess( transform_back(y) )
    assert( n == g.length )
    assert( n == h.size   )

    val G = Vec.tabulate(n){ i => trafos(i).deriv1(y{i}) }

    val H = h.copy
    H modify {
      (Hij, i,j) =>
        var result = Hij * G(i) * G(j)
        if( i == j )
          result += trafos(i).deriv2{y(i)} * g(i)
        result
    }

    G *= g

    (f,G,H)
  }
}