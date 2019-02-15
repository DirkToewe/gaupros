/* This file is part of GauProS.
 *
 * GauProS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GauProS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GauProS.  If not, see <https://www.gnu.org/licenses/>.
 */

package gps.regression

import gps.diff.gradient
import gps.doe.lhs
import gps.kernel._
import gps.linalg._
import gps.opt.ObjectiveWithHessian
import gps.regression.gpr.LikelihoodFunction
import utest._

import scala.util.Random

class GPR_logp_shifted_test1(
  name: String,
  likelihood         : (Array[Vec],Vec,Symbol,Kernel[Vec]) => LikelihoodFunction,
  likelihoodUnshifted: (Array[Vec],Vec,Double,Kernel[Vec]) => LikelihoodFunction
) extends TestSuite
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

        val logp = likelihood(x,y, 'y_shift, cov)

        var nFailures = 0

        def test( logpGrad: Vec => Vec )
          = for(
              _ <- 1 to 32;
              logpDiff = gradient(logp);
              _ <- 1 to 32;
              p = Vec(
                  1.0 + rng.nextDouble *128.0, // <- var_cov
                  0.5 + rng.nextDouble *128.0, // <- var_func
                  0.0 + rng.nextDouble *128.0, // <- var_noise
                -64.0 + rng.nextDouble *128.0  // <- y_shift
              )
            ) try {
              val d = logpDiff(p)
              val g = logpGrad(p)
              val (f,h) = likelihoodUnshifted(x,y, p(3), cov) fval_grad p.slice(0,3)
              assert{ allClose(d, g, tolAbs = 1e-3, tolRel = 1e-3) }
              assert{ allClose(h, g.slice(0,3) ) }

              assert( isClose(f, logp(p)              ) )
              assert( isClose(f, logp.fval_grad(p)._1 ) )
            }
            catch {
              case MatrixNotPositiveDefinite() => nFailures += 1
            }

        test(logp.gradient)
        test(logp.fval_grad(_)._2)
        logp match {
          case logp: ObjectiveWithHessian =>
            test(logp.fval_grad_hess(_)._2)
          case _ =>
        }

        assert( nFailures <= 3 )

        println( f"[GPR/logp${name}ShiftedGradient#1] run${run}%4d check!")
      }
    }
  }
}

class GPR_logp_shifted_test2(
  name: String,
  likelihood: (Array[Vec],Vec,Symbol,Kernel[Vec]) => LikelihoodFunction
) extends TestSuite
{
  override def tests = this{
    TestableSymbol('testRandomly) {
      for( run <- (0 to 16)/*.par*/ )
      {
        val rng = new Random(run)
        val nSamples = 2 + rng.nextInt(96)
        val nFeatures= 1 + rng.nextInt(8)

        val cov = Noise('var_noise)  +  ('var_func + 64.5) * Exp( - VecSum( AbsDelta.pow(1.9+0.11*rng.nextDouble) / 'var_cov ) )
        assert( 3 == cov.params.length )

        val x = lhs( nSamples, nFeatures, rng )
        x foreach (_ *= 4)
        x foreach (_ -= 2)

        assert(          x.length == nSamples )
        assert( x.forall(_.length == nFeatures) )

        val y = Vec.tabulate(nSamples)( _ => rng.nextDouble )

        val logp = likelihood(x,y, 'var_func, cov)

        var nFailures = 0

        def test( logpGrad: Vec => Vec )
          = for(
              _ <- 1 to 32;
              logpDiff = gradient(logp);
              _ <- 1 to 32;
              p = Vec(
                  1.0 + rng.nextDouble *128.0, // <- var_cov
                -64.0 + rng.nextDouble *128.0, // <- var_func
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
        logp match {
          case logp: ObjectiveWithHessian =>
            test(logp.fval_grad_hess(_)._2)
          case _ =>
        }

        assert( nFailures <= 3 )

        println( f"[GPR/logp${name}ShiftedGradient#2] run${run}%4d check!")
      }
    }
  }
}

object GPR_logp_marginal_shifted_test1 extends GPR_logp_shifted_test1("Marginal", GPR.logp_marginal(_,_,_,_), GPR.logp_marginal(_,_,_,_) )
object GPR_logp_marginal_shifted_test2 extends GPR_logp_shifted_test2("Marginal", GPR.logp_marginal(_,_,_,_) )
object GPR_logp_loo_shifted_test1      extends GPR_logp_shifted_test1("LOO",      GPR.logp_loo     (_,_,_,_), GPR.logp_loo     (_,_,_,_) )
object GPR_logp_loo_shifted_test2      extends GPR_logp_shifted_test2("LOO",      GPR.logp_loo     (_,_,_,_) )
