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

package gps

import gps.linalg.Vec
import scala.util.Random

/** A collection of Design of Experiments utility methods.
  *
  * Created by Dirk Toewe on 10.08.17.
  */
object doe
{
  /** Returns randomly sampled points using the latin hypercube sampling method.
    * The method guarantees that, when the points projected to each feature axis,
    * the points are always evenly spaced/distributed.
    *
    * @param nSamples
    * @param nFeatures
    * @param rng
    * @return
    */
  def lhs( nSamples: Int, nFeatures: Int, rng: Random = new Random() ): Array[Vec] =
  {
    assert( 1 < nSamples )
    val R = Array.tabulate(nSamples){ i => Vec.full(nFeatures)( i / (nSamples-1.0) ) }
    assert( R(0)         .forall(_ == 0.0) )
    assert( R(nSamples-1).forall(_ == 1.0) )

    // https://en.wikipedia.org/wiki/Fisher-Yates_shuffle
    var i = nSamples
    while( i > 0 )
    {
      i -= 1
      var j = nFeatures
      while( j > 0 )
      {
        j -= 1
        val k = rng.nextInt(i+1)
        val tmp = R(i)(j); R(i)(j) = R(k)(j); R(k)(j) = tmp
      }
    }
    R
  }
}
