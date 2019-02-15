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

import scala.annotation.tailrec

object Seqs
{
  /** An implementation of binary search for Scala indexed sequences.
    *
    * @param seq
    * @param key
    * @tparam E
    * @return
    */
  def binSearch[@specialized E]( seq: IndexedSeq[E], key: E )( implicit order: Ordering[E] ): Int =
  {
    @tailrec def binSearch( from: Int, to: Int ): Int
      = if( from > to )
          -(from+1)
        else {
          val mid = (from + to) / 2
          val c = order.compare(seq(mid),key)
               if( c < 0 ) binSearch(     mid+1,to)
          else if( c > 0 ) binSearch(from,mid-1   )
          else mid
        }
    binSearch(0,seq.size-1)
  }
}
