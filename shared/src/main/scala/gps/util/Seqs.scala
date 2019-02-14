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
