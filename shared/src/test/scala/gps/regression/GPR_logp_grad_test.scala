package gps.regression

import gps.diff.gradient
import gps.doe.lhs
import gps.kernel._
import gps.linalg._
import gps.regression.gpr.LikelihoodFunction
import utest._

import scala.util.Random

class GPR_logp_grad_test( name: String, likelihood: (Array[Vec],Vec,Double,Kernel[Vec]) => LikelihoodFunction[Vec] ) extends TestSuite
{
  override def tests = this{
    TestableSymbol('testRandomly) {
      for( run <- (0 to 16)/*.par*/ )
      {
        val rng = new Random(run)
        val nSamples = 2 + rng.nextInt(96)
        val nFeatures= 1 + rng.nextInt(8)

        val cov = Noise('var_noise)  +  'var_func * Exp( - VecSum( AbsDelta.pow(1.9+0.11*rng.nextDouble) / 'var_cov ) )
        assert( 3 == cov.params.length )

        val x = lhs( nSamples, nFeatures, rng )
        x foreach (_ *= 4)
        x foreach (_ -= 2)

        assert(          x.length == nSamples )
        assert( x.forall(_.length == nFeatures) )

        val y = Vec.tabulate(nSamples)( _ => rng.nextDouble )

        val logp = likelihood(x,y, /*yShift=*/-1+2*rng.nextDouble, cov)

        var nFailures = 0

        def test( logpGrad: Vec => Vec )
        = for(
          _ <- 1 to 32;
          logpDiff = gradient(logp);
          _ <- 1 to 32;
          p = Vec(
            1.0 + rng.nextDouble *128.0, // <- var_cov
            0.5 + rng.nextDouble *128.0, // <- var_func
            0.0 + rng.nextDouble *128.0  // <- var_noise
          )
        ) try {
          val d = logpDiff(p)
          val g = logpGrad(p)
          assert{ allClose(d, g, tolAbs = 1e-3, tolRel = 1e-3) }
        }
        catch {
          case MatrixNotPositiveDefinite() => nFailures += 1
        }

        test(logp.gradient)
        test(logp.fval_grad(_)._2)
        test(logp.fval_grad_hess(_)._2)

        assert( nFailures <= 3 )

        println( f"[GPR/logp${name}Gradient] run${run}%4d check!")
      }
    }
  }
}
object GPR_logp_grad_marginal_test extends GPR_logp_grad_test("Marginal", GPR.logp_marginal)
object GPR_logp_grad_loo_test      extends GPR_logp_grad_test("LOO",      GPR.logp_loo     )