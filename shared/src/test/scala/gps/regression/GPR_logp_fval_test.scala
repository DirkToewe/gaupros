package gps.regression

import gps.doe.lhs
import gps.kernel._
import gps.linalg._
import gps.opt.ObjectiveWithHessian
import utest._

import scala.util.Random

object GPR_logp_fval_test extends TestSuite
{
  private def testLikelihoodFVal(
    name: String,
    likelihood: (Array[Vec],Vec,Double,Kernel[Vec]) => ObjectiveWithHessian
  ) = {
    for( run <- (0 to 32)/*.par*/ )
    {
      val rng = new Random(run)
      val nSamples = 2 + rng.nextInt(96)
      val nFeatures= 1 + rng.nextInt(8)

      val cov = Noise('var_noise)  +  'var_func * Exp( - VecSum( AbsDelta.pow(1.9+0.11*rng.nextDouble) / 'var_cov ) )
      assert( 3 == cov.params.length )

      val x = lhs( nSamples, nFeatures, rng )
      x.foreach(_ *= 4)
      x.foreach(_ -= 2)

      assert(           x.length == nSamples   )
      assert{ x forall (_.length == nFeatures) }

      val y = Vec.tabulate(nSamples)( _ => rng.nextDouble )

      val logp = likelihood(x,y, /*yShift=*/-1+2*rng.nextDouble, cov)

      var nFailures = 0

      def test( fval: Vec => Double )
        = for(
            _ <- 1 to 32;
            p = Vec(
              1.0 + rng.nextDouble *128.0, // <- var_cov
              0.5 + rng.nextDouble *128.0, // <- var_func
              0.0 + rng.nextDouble *128.0  // <- var_noise
            )
          ) try {
            assert{ isClose(logp(p), fval(p)) }
          }
          catch {
            case MatrixNotPositiveDefinite() => nFailures += 1
          }

      test(logp.fval_grad(_)._1)
      test(logp.fval_grad_hess(_)._1)

      assert( nFailures <= 2 )

      println( f"[GPR/logp${name}FVal] run${run}%4d check!")
    }
  }

  override def tests = this {
    TestableSymbol('logpLOOFVal     ) { testLikelihoodFVal("LOO",      GPR.logp_loo     (_,_,_,_)) }
    TestableSymbol('logpMarginalFVal) { testLikelihoodFVal("Marginal", GPR.logp_marginal(_,_,_,_)) }
  }
}