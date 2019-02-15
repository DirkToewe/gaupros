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

package gps.kernel

import gps.diff.gradient
import gps.linalg._
import utest._

import scala.math.{exp, log, pow}
import scala.util.Random

/**
  * Created by dtitx on 05.08.17.
  */
object Kernel_test extends TestSuite
{
  override def tests = this {

    TestableSymbol('krigingKernel) {

      val cov: Vec => (Vec, Vec, Int, Int) => Double = {
        case Vec(var_noise, var_func, norm_p, t0, t1, t2) =>
          (xi, xj, i, j) =>
            val dx = xi - xj map (_.abs) map { pow(_, norm_p) }
            var result = var_func * exp(-dx(0) / t0 - dx(1) / t1 - dx(2) / t2) ensuring (!_.isNaN)
            if (0 <= i && i == j)
              result += var_noise
            result
      }

      def covDiff(p: Vec)(xi: Vec, xj: Vec, i: Int, j: Int)
      = gradient(cov(_)(xi, xj, i, j))(p)

      val covGrad: Vec => Int => (Vec, Vec, Int, Int) => Double = {
        case Vec(var_noise, var_func, norm_p, t0, t1, t2) =>
          if (var_noise < 0) throw new IllegalArgumentException
          else if (var_func <= 0) throw new IllegalArgumentException
          else if (t0 <= 0) throw new IllegalArgumentException
          else if (t1 <= 0) throw new IllegalArgumentException
          else if (t2 <= 0) throw new IllegalArgumentException
          else {
            case 0 /*var_noise*/ => (_, _, i, j) => if (0 <= i && i == j) 1 else 0
            case 1 /*var_func*/ => (xi, xj, _, _) =>
              val dx = xi - xj map (_.abs) map { pow(_, norm_p) }
              exp(-dx(0) / t0 - dx(1) / t1 - dx(2) / t2)
            case 2 /*norm_p*/ => (xi, xj, _, _) =>
              val dx = xi - xj map (_.abs) map { pow(_, norm_p) }
              val dxdp = xi - xj map (_.abs) map { x => pow(x, norm_p) * log(x) }
              exp(-dx(0) / t0 - dx(1) / t1 - dx(2) / t2) * var_func *
                (-dxdp(0) / t0 - dxdp(1) / t1 - dxdp(2) / t2)
            case 3 /*t0*/ => (xi, xj, _, _) =>
              val dx = xi - xj map (_.abs) map { pow(_, norm_p) }
              exp(-dx(0) / t0 - dx(1) / t1 - dx(2) / t2) * var_func * (dx(0) / t0 / t0)
            case 4 /*t1*/ => (xi, xj, _, _) =>
              val dx = xi - xj map (_.abs) map { pow(_, norm_p) }
              exp(-dx(0) / t0 - dx(1) / t1 - dx(2) / t2) * var_func * (dx(1) / t1 / t1)
            case 5 /*t1*/ => (xi, xj, _, _) =>
              val dx = xi - xj map (_.abs) map { pow(_, norm_p) }
              exp(-dx(0) / t0 - dx(1) / t1 - dx(2) / t2) * var_func * (dx(2) / t2 / t2)
          }
      }

      val rng = new Random(1337)

      for (
        _ <- 1 to 512;
        t0 = 1e-6 + 10 * rng.nextDouble; var_noise= 1e-6 + 10 * rng.nextDouble;
        t1 = 1e-6 + 10 * rng.nextDouble; var_func = 1e-6 + 10 * rng.nextDouble;
        t2 = 1e-6 + 10 * rng.nextDouble; norm_p   = 1.9  + 0.2* rng.nextDouble;
        p = Vec(var_noise, var_func, norm_p, t0, t1, t2);
        _ <- 1 to 256;
        i = -1 + rng.nextInt(32); xi = Vec.tabulate(3)( _ => -10 + 20 * rng.nextDouble );
        j = -1 + rng.nextInt(32); xj = Vec.tabulate(3)( _ => -10 + 20 * rng.nextDouble );
        cd = covDiff(p)(xi, xj, i, j);
        cg = Vec.tabulate(cd.length)( covGrad(p)(_)(xi, xj, i, j) )
      ) {
        assert{ allClose(cd,cg,1e-4,1e-6) }
      }

      val normPow = AbsDelta pow 'norm_p
      val kernel = Noise('var_noise) + 'var_func * Exp(
        - VecIdx(0, normPow / 't0)
        - VecIdx(1, normPow / 't1)
        - VecIdx(2, normPow / 't2)
      )
      assert(kernel.params.length == 6)

      for (
        _ <- 1 to 512;
        t0 = 1e-6 + 10 * rng.nextDouble; var_noise= 1e-6 + 10 * rng.nextDouble;
        t1 = 1e-6 + 10 * rng.nextDouble; var_func = 1e-6 + 10 * rng.nextDouble;
        t2 = 1e-6 + 10 * rng.nextDouble; norm_p   = 1.9  + 0.2* rng.nextDouble;
        p = Vec(var_noise, var_func, norm_p, t0, t1, t2);
        _ <- 1 to 256;
        i = -1 + rng.nextInt(32); xi = Vec.tabulate(3)( _ => -10 + 20 * rng.nextDouble );
        j = -1 + rng.nextInt(32); xj = Vec.tabulate(3)( _ => -10 + 20 * rng.nextDouble );
        fc = cov(p)(xi, xj, i, j);
        fk = kernel subs(
          'var_noise -> var_noise, 'var_func -> var_func, 'norm_p -> norm_p,
          't0 -> t0, 't1 -> t1, 't2 -> t2
        ) apply (xi, xj, i, j)
      ) {
        assert{ isClose(fc,fk) }
      }

      for (
        _ <- 1 to 512;
        t0 = 1e-6 + rng.nextDouble() * 10; var_noise= 1e-6 + 10 * rng.nextDouble;
        t1 = 1e-6 + rng.nextDouble() * 10; var_func = 1e-6 + 10 * rng.nextDouble;
        t2 = 1e-6 + rng.nextDouble() * 10; norm_p   = 1.9  + 0.2* rng.nextDouble;
        p = Vec(var_noise, var_func, norm_p, t0, t1, t2);
        _ <- 1 to 256;
        i = -1 + rng.nextInt(32); xi = Vec.tabulate(3)( _ => -10 + 20 * rng.nextDouble );
        j = -1 + rng.nextInt(32); xj = Vec.tabulate(3)( _ => -10 + 20 * rng.nextDouble );
        (sym, k) <- Seq('var_noise, 'var_func, 'norm_p, 't0, 't1, 't2).zipWithIndex;
        cg = covGrad(p)(k)(xi, xj, i, j);
        kg = kernel pDiff sym subs(
          'var_noise -> var_noise, 'var_func -> var_func, 'norm_p -> norm_p,
          't0 -> t0, 't1 -> t1, 't2 -> t2
        ) apply (xi, xj, i, j)
      ) {
        assert{ isClose(cg,kg) }
      }
    }



    TestableSymbol('powExpKernel) {

      val cov: Vec => (Vec, Vec, Int, Int) => Double = {
        case Vec(var_noise, var_func, norm_p, theta) =>
          (xi, xj, i, j) =>
            val dx = xi - xj map (_.abs) map { pow(_, norm_p) }
            var result = var_func * exp(-dx.sum * theta) ensuring (!_.isNaN)
            if (0 <= i && i == j)
              result += var_noise
            result
      }

      def covDiff(p: Vec)(xi: Vec, xj: Vec, i: Int, j: Int)
        = gradient(cov(_)(xi, xj, i, j))(p)

      val covGrad: Vec => Int => (Vec, Vec, Int, Int) => Double = {
        case Vec(var_noise, var_func, norm_p, theta) =>
          if (var_noise < 0) throw new IllegalArgumentException
          else if (var_func <= 0) throw new IllegalArgumentException
          else if (theta < 0) throw new IllegalArgumentException
          else {
            case 0 /*var_noise*/ => (_, _, i, j) => if (0 <= i && i == j) 1 else 0
            case 1 /*var_func*/ => (xi, xj, _, _) =>
              val dx = xi - xj map (_.abs) map { pow(_, norm_p) }
              exp(-dx.sum * theta)
            case 2 /*norm_p*/ => (xi, xj, _, _) =>
              val dx = xi - xj map (_.abs) map { pow(_, norm_p) }
              val dxdp = xi - xj map (_.abs) map { x => pow(x, norm_p) * log(x) }
              exp(-dx.sum * theta) * var_func * (-dxdp.sum * theta)
            case 3 /*theta*/ => (xi, xj, _, _) =>
              val dx = xi - xj map (_.abs) map { pow(_, norm_p) }
              exp(-dx.sum * theta) * var_func * -dx.sum
          }
      }

      val rng = new Random(1337)

      for (
        _ <- 1 to 512;
        theta    = 1e-6 + 10 * rng.nextDouble;
        var_noise= 1e-6 + 10 * rng.nextDouble;
        var_func = 1e-6 + 10 * rng.nextDouble;
        norm_p   = 1.9  + 0.2* rng.nextDouble;
        p = Vec(var_noise, var_func, norm_p, theta);
        _ <- 1 to 256;
        i = -1 + rng.nextInt(32); xi = Vec.tabulate(3)( _ => -10 + 20 * rng.nextDouble );
        j = -1 + rng.nextInt(32); xj = Vec.tabulate(3)( _ => -10 + 20 * rng.nextDouble );
        k <- 0 until 4;
        cd = covDiff(p)(xi, xj, i, j)(k);
        cg = covGrad(p)(k)(xi, xj, i, j)
      ) assert{ isClose(cd,cg, tolAbs=1e-4, tolRel=1e-6) }

      val normPow = AbsDelta pow 'norm_p
      val kernel = Noise('var_noise) + 'var_func * Exp( -VecSum('theta * normPow) )
      assert(kernel.params.length == 4)

      for (
        _ <- 1 to 512;
        theta    = 1e-6 + 10 * rng.nextDouble;
        var_noise= 1e-6 + 10 * rng.nextDouble;
        var_func = 1e-6 + 10 * rng.nextDouble;
        norm_p   = 1.9  + 0.2* rng.nextDouble;
        p = Vec(var_noise, var_func, norm_p, theta);
        _ <- 1 to 256;
        i = -1 + rng.nextInt(32); xi = Vec.tabulate(3)( _ => -10 + 20 * rng.nextDouble );
        j = -1 + rng.nextInt(32); xj = Vec.tabulate(3)( _ => -10 + 20 * rng.nextDouble );
        fc = cov(p)(xi, xj, i, j);
        fk = kernel subs(
          'var_noise -> var_noise, 'var_func -> var_func, 'norm_p -> norm_p, 'theta -> theta
        ) apply (xi, xj, i, j)
      ) assert{ isClose(fc,fk) }

      for (
        _ <- 1 to 512;
        theta    = 1e-6 + 10 * rng.nextDouble;
        var_noise= 1e-6 + 10 * rng.nextDouble;
        var_func = 1e-6 + 10 * rng.nextDouble;
        norm_p   = 1.9  + 0.2* rng.nextDouble;
        p = Vec(var_noise, var_func, norm_p, theta);
        _ <- 1 to 256;
        i = -1 + rng.nextInt(32); xi = Vec.tabulate(3)( _ => -10 + 20 * rng.nextDouble );
        j = -1 + rng.nextInt(32); xj = Vec.tabulate(3)( _ => -10 + 20 * rng.nextDouble );
        (sym, k) <- Seq('var_noise, 'var_func, 'norm_p, 'theta).zipWithIndex;
        cg = covGrad(p)(k)(xi, xj, i, j);
        kg = kernel pDiff sym subs(
          'var_noise -> var_noise, 'var_func -> var_func, 'norm_p -> norm_p, 'theta -> theta
        ) apply (xi, xj, i, j)
      ) assert{ isClose(kg,cg) }
    }

  }
}