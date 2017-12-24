package gps.linalg

import utest._

import scala.util.Random

/**
  * Created by dtub on 19.06.17.
  */
object LMatRM_test extends LMat_test(LMatRM)
object LMatRM_extended_test extends TestSuite
{
  override def tests = this {


    'tabulateRM {
      val rng = new Random(1337)
      for (run <- 0 until 128) {
        val size = rng.nextInt(128) + 1
        val L = IndexedSeq.tabulate(size,size){ (_,_) => rng.nextDouble }
        LMatRM.tabulateRM(size)( i => j => L(i)(j) ) foreach { (Lij,i,j) => assert( Lij == L(i)(j) ) }
      }
    }


    'choleskyAddRowTest {
      val rng = new Random(1337)
      for (run <- 0 until 128) {
        val size = rng.nextInt(128) + 1

        val L = Array.tabulate(size){
          i =>    Vec.tabulate(size){
          j =>
            if (i > j)
              rng.nextDouble() * 200 - 100
            else if (i == j)
              rng.nextDouble() * 950 + 50
            else
              0
        }}

        val LLT = Array.tabulate(size){ i => Vec.tabulate(size){L(i) ⋅ L(_)}}

        val ref = LMatRM.tabulate(size){ LLT(_)(_) }
        ref.choleskyDecomp()

        val mat = LMatRM.tabulate(0)(null)
        for (s <- 0 until size)
          mat.choleskyRowAppend( Vec.tabulate(s + 1){ LLT(s)(_) } )

        for (
          r <- 0 until size;
          c <- 0 until size
        )
          assert(ref(r, c) == mat(r, c))
        println( f"[${LMatRM.getClass.getSimpleName}/choleskyDecomp] run$run%4d check!")
      }
    }
  }
}