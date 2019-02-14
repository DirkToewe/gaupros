package gps.regression

import gps.diff.jacobian
import gps.doe.lhs
import gps.kernel._
import gps.linalg._
import gps.opt.ObjectiveWithHessian
import gps.regression.gpr.LikelihoodFunction
import utest._

import scala.util.Random

class GPR_logp_hess_test( name: String, likelihood: (Array[Vec],Vec,Double,Kernel[Vec]) => LikelihoodFunction with ObjectiveWithHessian ) extends TestSuite
{
  override def tests = this{
    TestableSymbol('testRandomly){
      for( run <- (0 to 16)/*.par*/ )
      {
        val rng = new Random(run)
        val nSamples = 2 + rng.nextInt(96)
        val nFeatures = 2 + rng.nextInt(8)

        val cov = Noise('var_noise) + 'var_func * Exp(
          - VecIdx(0, AbsDelta.pow(1.9 + 0.11 * rng.nextDouble) * 't0 )
            - VecIdx(1, AbsDelta.pow(1.9 + 0.11 * rng.nextDouble) / 't1 )
        )

        val x = lhs(nSamples, nFeatures, rng)
        x.foreach(_ *= 4)
        x.foreach(_ -= 2)

        assert(x.length == nSamples)
        assert(x.forall(_.length == nFeatures))

        val y = Vec.tabulate(nSamples)( _ => rng.nextDouble )

        val logp = likelihood(x, y, /*yShift=*/ -1 + 2 * rng.nextDouble, cov)

        var nFailures = 0

        def test( logpHess: Vec => LMat )
        = for (
          _ <- 1 to 32;
          logpDiff = jacobian(logp.gradient);
          _ <- 1 to 32;
          p = Vec(
            0.0 + rng.nextDouble * 2.0, // <- t0
            0.5 + rng.nextDouble * 1024,// <- t1
            0.5 + rng.nextDouble * 128, // <- var_func
            0.0 + rng.nextDouble * 128  // <- var_noise
          )
        ) try {
          val d = logpDiff(p)
          val h = logpHess(p)
          h foreach {
            (h_ij, i, j) => assert(
              isClose(h_ij, {d(i,j) + d(j,i)} / 2, 1e-3, 1e-3)
            )
          }
        }
        catch {
          case MatrixNotPositiveDefinite() => nFailures += 1
        }

        test(logp.hessian(_))
        test(logp.fval_grad_hess(_)._3)

        assert( nFailures <= 2 )

        println( f"[GPR/logp${name}Hessian] run${run}%4d check!")
      }
    }
  }
}
object GPR_logp_hess_marginal_test extends GPR_logp_hess_test("Marginal", GPR.logp_marginal(_,_,_,_) )
object GPR_logp_hess_loo_test      extends GPR_logp_hess_test("LOO",      GPR.logp_loo     (_,_,_,_) )