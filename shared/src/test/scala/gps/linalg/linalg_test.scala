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

package gps.linalg

import gps.doe.lhs
import gps.kernel._
import gps.linalg
import utest.{TestSuite, TestableSymbol}

import scala.math.{cos, sin, toRadians => deg2rad}
import scala.util.Random

/**
  * Created by Dirk Toewe on 23.08.17.
  */
object linalg_test extends TestSuite
{
  override def tests = this{

    TestableSymbol('norm) {

      val rng = new Random(1337)

      for( _ <- 0 to 1024*1024 )
      {
        val r = rng.nextDouble*4
        val θ = deg2rad(rng.nextDouble*180)
        val ϕ = deg2rad(rng.nextDouble*360)
        val u = Vec(
          r * sin(θ) * cos(ϕ),
          r * sin(θ) * sin(ϕ),
          r * cos(θ)
        )
        val nrm= norm(u)
        val np = normPow(u)
        assert{ isClose(nrm,r) }
        assert{ isClose(np,r*r) }
      }
    }


    TestableSymbol('sum) {
      assert{
        isClose( 1, sum(1000*1000){Array.fill(1000*1000)(1e-6)} )
      }

      val rng = new Random(1337)

      for( _ <- 0 to 1024 )
      {
        val n = rng.nextInt(1024)
        val summands = Array.fill(n){rng.nextInt(32*1024) - 16L}
        assert{
          isClose(summands.sum, sum(n){summands(_)} )
        }
      }
    }


    TestableSymbol('cg_jac) {

      val rng = new Random(1337)

      for( run <- (0 until 32)/*.par*/ )
      {
        val nSamples = rng.nextInt(1024)+2
        val nFeatures= rng.nextInt(   8)+1

        val x = lhs(nSamples,nFeatures,rng)
        x.foreach(_ *= 16)
        x.foreach(_ -=  8)

        val noise   = rng.nextDouble*4  + 0.1
        val variance= rng.nextDouble*31 + 1
        val cov = Noise(noise) + variance * Exp(
          - (0 until nFeatures)
              .map{ VecIdx(_, AbsDelta) }
              .map{ _.pow{rng.nextDouble*0.11 + 1.9} }
              .map( _ * rng.nextDouble*4 )
              .reduce(_ + _)
        )

//        val A = (i: Int, j: Int) => cov(x{i},x{j}, i,j)
        val A = Mat.tabulateSym(nSamples){ (i,j) => cov(x{i},x{j}, i,j) }
        for( _ <- 0 until 128 )
        {
          val b = Vec.tabulate(nSamples)( _ => rng.nextDouble()*16 - 8 )
          val x = linalg.cg_jac(A,b, tolRel=1e-4, tolAbs=1e-5 )

          for( i <- 0 until nSamples )
          {
            val β = (0 until nSamples).map{ j => A(i,j)*x(j) }.sum
            val bi = b(i)
            assert( isClose(β,bi, tolRel=1e-3, tolAbs=1e-5) )
          }
        }
        println( f"[cg_jac] run$run%4d check!")
      }
    }
  }
}
