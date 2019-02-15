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

import utest._

import scala.util.Random

/**
  * Created by dtitx on 06.08.17.
  */
class Mat_test( mat: MatFactory ) extends TestSuite
{
  override def tests = this{
    'tabulate {
      for( nRows <- (1 until 64)/* ++ (MAX_SIZE-16 to MAX_SIZE) */ )
      {
        for( nCols <- (1 until 64)/* ++ (MAX_SIZE-16 to MAX_SIZE) */ )
        {
          val M = mat.tabulate(nRows,nCols)( (i,j) => nCols*(i+1) + j+1 )
          for(
            r <- 0 until nRows;
            c <- 0 until nCols
          ) assert( nCols*(r+1) + (c+1)  ==  M(r,c) )
        }
        println( f"[${mat.getClass.getSimpleName}/tabulate] run${nRows}%4d check!")
      }
    }


    'foreach {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val nRows = rng.nextInt(128) + 1
        val nCols = rng.nextInt(128) + 1
        val L = IndexedSeq.tabulate(nRows, nCols) {
          (i,j) => rng.nextInt(1000)
        }
        val M = mat.tabulate(nRows,nCols){L(_)(_)}

        var sum = 0L
        M foreach {
          (mat_ij,i,j) =>
            sum += mat_ij.toInt
            assert( mat_ij == L{i}{j} )
        }
        assert( sum == L.map(_.sum).sum )

        println( f"[${mat.getClass.getSimpleName}/foreach] run$run%4d check!")
      }
    }


    'diagReduce {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128) + 1
        val A = mat.tabulate(size,size){ (_,_) => rng.nextInt(200) - 100 }

        val sum = A diagReduce (_ + _)
        val check = (0 until size).map{ i => A(i,i) }.sum
        assert{ isClose(sum,check) }

        println( f"[${mat.getClass.getSimpleName}/diagReduce] run$run%4d check!")
      }
    }


    'mapReduce {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val nRows = rng.nextInt(128) + 1
        val nCols = rng.nextInt(128) + 1
        val L = IndexedSeq.tabulate(nRows, nCols) {
          (i,j) => -100 + rng.nextInt(200)
        }
        val M = mat.tabulate(nRows,nCols){L(_)(_)}

        val sum = M.mapReduce{
          (mat_ij,i,j) =>
            assert( mat_ij == L{i}{j} )
            mat_ij.toInt
        }(_ + _)
        assert( sum == L.map(_.sum).sum )

        println( f"[${mat.getClass.getSimpleName}/mapReduce] run$run%4d check!")
      }
    }


    'matMult {
      val rng = new Random(1337)
      def test( matR: MatFactory )
        = for( run <- 0 until 128 )
          {
            val A = mat .tabulate(rng.nextInt(128) + 1, rng.nextInt(128) + 1)( (_,_) => -10 + 20*rng.nextDouble )
            val B = matR.tabulate(A.nCols,              rng.nextInt(128) + 1)( (_,_) => -10 + 20*rng.nextDouble )
            val C = A *# B

            C.foreach{
              (Cij,i,j) =>
                val check = (0 until A.nCols).map{ k => A(i,k)*B(k,j) }.sum
                assert( isClose(Cij,check) )
            }

            println( f"[${mat.getClass.getSimpleName} *# ${matR.getClass.getSimpleName}] run$run%4d check!")
          }
      test(MatRM)
      test(MatCM)
    }


    'vecMult {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val A = mat.tabulate(rng.nextInt(128) + 1, rng.nextInt(128) + 1)( (_,_) => -10 + 20*rng.nextDouble )
        val b = Vec.tabulate(A.nCols)( _ => -10 + 20*rng.nextDouble )
        val c = A *# b

        assert( c.length == A.nRows )
        for( i <- 0 until A.nRows )
        {
          val c_i = (0 until A.nCols).map( j => A(i,j) * b(j) ).sum
          assert{ isClose(c(i), c_i) }
        }

        println( f"[${mat.getClass.getSimpleName} *# Vec] run$run%4d check!")
      }
    }


    'addInplace {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val A = mat.tabulate(rng.nextInt(128) + 1, rng.nextInt(128) + 1)( (_,_) => -10 + 20*rng.nextDouble )
        val B = mat.tabulate(A.nRows,              A.nCols             )( (_,_) => -10 + 20*rng.nextDouble )
        val C = A.copy

        C += B

        for( i <- 0 until A.nRows )
          for( j <- 0 until A.nCols )
            assert{ isClose( C(i,j), A(i,j) + B(i,j) ) }

        println( f"[${mat.getClass.getSimpleName} += ${mat.getClass.getSimpleName}] run$run%4d check!")
      }
    }
  }
}
object Mat_test extends Mat_test(Mat)
object Mat_extended_test extends TestSuite
{
  override def tests = this{
    'inner {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val nRows = rng.nextInt(128) + 1
        val nCols = rng.nextInt(128) + 1
        val u = Vec.tabulate(nRows)( _ => -10 + 20*rng.nextDouble )
        val v = Vec.tabulate(nCols)( _ => -10 + 20*rng.nextDouble )
        val A = Mat.tabulate(nRows,nCols){ (_,_) => -10 + 20*rng.nextDouble }

        var inner = Mat.inner(u,A,v)
        for(
          r <- 0 until nRows;
          c <- 0 until nCols
        ) inner -= u(r) * A(r,c) * v(c)

        assert( isClose(inner,0) )

        println( f"Mat/inner run$run%4d check!")
      }
    }
  }
}