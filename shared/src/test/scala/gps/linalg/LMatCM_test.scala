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
  * Created by dtub on 19.06.17.
  */
object LMatCM_test extends LMat_test(LMatCM)
object LMatCM_extended_test extends TestSuite
{
  override def tests = this{


    'symMatProdTrace {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128)+1

        val A = Array.tabulate(size){ i => Vec.tabulate(size){ j => if( i >= j ) rng.nextDouble()*200-100 else 0 }}
        val B = Array.tabulate(size){ i => Vec.tabulate(size){ j => if( i >= j ) rng.nextDouble()*200-100 else 0 }}
        for( row <- 0 until size;
             col <- 0 until row )
        {
          A(col)(row) = A(row)(col)
          B(col)(row) = B(row)(col)
        }
        for( row <- 0 until size;
             col <- 0 until size )
        {
          assert( A(row)(col) == A(col)(row) )
          assert( B(row)(col) == B(col)(row) )
        }

        val a = LMatCM.tabulate(size){A(_)(_)}
        val b = LMatCM.tabulate(size){B(_)(_)}

        val trRef = (0 until size).map( i => A(i) ⋅ B(i) ).sum
        val tr = LMatCM.symMatProdTrace(a,b)

        assert{ 1e-6 > (tr - trRef).abs }

        println( f"[${LMatCM.getClass.getSimpleName}/symMatProdTrace(_,_)] run$run%4d check!")
      }
    }


    'tabulateCM {
      val rng = new Random(1337)
      for (run <- 0 until 128)
      {
        val size = rng.nextInt(128) + 1
        val L = IndexedSeq.tabulate(size,size){ (_,_) => rng.nextDouble }
        LMatCM.tabulateCM(size)( j => i => L(i)(j) ) foreach { (Lij,i,j) => assert( Lij == L(i)(j) ) }
      }
    }


  }
}