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

package gps.util

import gps.util.Seqs.binSearch
import utest._

import scala.util.Random

object Seqs_test extends TestSuite
{
  override def tests = this{

    'binSearch {
      val rng = new Random(1337)
      import rng._

      for( _ <- 0 until 1024 )
      {
        val len = nextInt(1024)
        val seq = IndexedSeq.iterate( start = nextInt{128}-64, len=len ){ nextInt(8) + _ }
        val min   = if(0==len) -1 else seq.min-1
        val range = if(0==len)  2 else seq.max+1 - min

        for( _ <- 0 until 1024 )
        {
          val key = nextInt(range) + min
          var i = binSearch(seq,key)
          if( i >= 0 )
            assert(seq(i) == key)
          else {
            i = ~i
            assert( i == 0          || seq(i-1) < key )
            assert( i == seq.length || seq(i)   > key )
          }
        }
      }
    }

  }
}