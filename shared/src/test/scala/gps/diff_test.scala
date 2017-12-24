package gps

import java.lang.Math.{pow, sqrt}

import gps.linalg.{Mat, MatRM, Vec}
import utest._

import scala.collection.immutable.{IndexedSeq => ISeq}
import scala.util.Random

/**
  * Created by Dirk Toewe on 16/07/17.
  */
object diff_test extends TestSuite
{
  def tests = this{
    'gradientTest {
      val rng = new Random(1337)

      for( run <- 0 until 128 )
      {
        val maxExp = 3
        val nVars = 1+rng.nextInt(4)
        val nCoefs = Iterator.fill(nVars)(maxExp).product // 3**nVars (dumb implementation i know)
        val coefs = Array.fill(nCoefs)( rng.nextFloat()*2 - 1.0 )

        // creates a random full factorial polynomial
        val exponents = {
          def exponents( nVars: Int ): Iterable[ISeq[Int]]
            = if( 0 == nVars )
                Iterable(Vector.empty)
              else for(
                vec <- exponents(nVars-1);
                i <- 0 until maxExp
              ) yield vec :+ i
          exponents(nVars)
        }

        assert( nCoefs == exponents.size )

        def f( x: Vec ): Double = {
          for( (e,c) <- exponents zip coefs )
            yield c * (x zip e).map{ case (b,e) => pow(b,e) }.product
        }.sum ensuring x.length == nVars

        def fGrad( x: Vec )
          = Vec.tabulate(nVars){ i => {
              for( (e,c) <- exponents zip coefs )
                yield c * e(i) * e.updated(i, e(i)-1 max 0).zip(x).map{ case (e,b) => pow(b,e) }.product
            }.sum } ensuring x.length == nVars

        val fDiff = diff.gradient(f)

        for( test <- 0 until 128 )
        {
          val x = Vec.tabulate(nVars)( _ => rng.nextFloat()*4 - 2.0 )
          val dist = linalg.distance(fGrad{x}, fDiff{x}) / sqrt(nVars)
          assert{ 1e-3 > dist }
        }

        println(f"[diff/gradient] run$run%4d check!")
      }
    }
    'jacobianTest {
      val rng = new Random(1337)

      for( run <- 0 until 128 )
      {
        val maxExp = 3
        val nVars = 1+rng.nextInt(4)
        val nCoefs = Iterator.fill(nVars)(maxExp).product // 3**nVars (dumb implementation i know)
        val nOutputs = 1+rng.nextInt(16)
        val coefs = Array.fill(nOutputs,nCoefs)( rng.nextFloat()*2 - 1.0 )

        // creates a random full factorial polynomial
        val exponents = {
          def exponents( nVars: Int ): Iterable[ISeq[Int]]
            = if( 0 == nVars )
                Iterable(Vector.empty)
              else for(
                vec <- exponents(nVars-1);
                i <- 0 until maxExp
              ) yield vec :+ i
          exponents(nVars)
        }

        assert( nCoefs == exponents.size )

        def f( x: Vec ) = Vec.tabulate(nOutputs){
          i => {
            for( (e,c) <- exponents zip coefs(i) )
              yield c * (x zip e).map{ case (b,e) => pow(b,e) }.product
          }.sum ensuring x.length == nVars
        }

        def fJac( x: Vec ): Mat
          = MatRM.tabulate(nOutputs,nVars){
              (i,j) => {
                for( (e,c) <- exponents zip coefs(i) )
                  yield c * e(j) * e.updated(j, e(j)-1 max 0).zip(x).map{ case (e,b) => pow(b,e) }.product
              }.sum ensuring x.length == nVars
            }

        val fDiff = diff.jacobian(f)

        def matDist( a: Mat, b: Mat ): Double =
        {
          var result = 0.0
          val nEntries = a.nRows * b.nCols
          a.foreach{ (a_ij,i,j) => val dx = b(i,j) - a_ij; result += dx*dx }
          sqrt(result / nEntries)
        }

        for( test <- 0 until 128 )
        {
          val x = Vec.tabulate(nVars)( _ => rng.nextFloat()*4 - 2.0 )
          val dist = matDist(fJac{x}, fDiff{x})
          assert{ 1e-3 > dist }
        }

        println(f"[diff/jacobian] run$run%4d check!")
      }
    }
  }
}
