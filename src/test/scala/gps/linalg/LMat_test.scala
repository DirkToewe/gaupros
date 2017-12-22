package gps.linalg

import utest._
import LMat.MAX_SIZE

import scala.util.Random

/**
  * Created by dtult on 15/07/17.
  */
object LMat_test extends LMat_test(LMat) {}
class LMat_test( lMat: LMatFactory ) extends TestSuite
{
  override def tests = this{

    'tabulate {
      for( size <- (1 until 128)/* ++ (MAX_SIZE-16 to MAX_SIZE) */ )
      {
        val mat = lMat.tabulate(size)(size+1 + size*_ + _)
        for(
          r <- 0 until size;
          c <- 0 until size
        ) if( c > r )
          assert( 0 == mat(r,c) )
        else
          assert( size+1 + size*r + c  ==  mat(r,c) )
        println( f"[${lMat.getClass.getSimpleName}/tabulate] run${size}%4d check!")
      }
    }


    'modifyAll {
      for( size <- 1 until 128 )
      {
        val mat = lMat.tabulate(size)(11 + 10*_ + _)
        mat.modify( (x,i,j) => 2*x + 11 + 10*i + j )
        for(
          r <- 0 until size;
          c <- 0 until size
        ) if( c > r )
            assert( 0 == mat(r,c) )
          else
            assert( 33 + 30*r + 3*c  ==  mat(r,c) )
        println( f"[${lMat.getClass.getSimpleName}/modifyAll] run${size}%4d check!")
      }
    }


    'foreach {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128) + 1
        val L = IndexedSeq.tabulate(size, size) {
          (i,j) => if( i >= j ) rng.nextInt(100) else 0L
        }
        val mat = lMat.tabulate(size){L(_)(_)}

        var sum = 0L
        mat foreach {
          (mat_ij,i,j) =>
            sum += mat_ij.toInt
            assert( mat_ij == L{i}{j} )
        }
        assert( sum == L.map(_.sum).sum )

        println( f"[${lMat.getClass.getSimpleName}/mapReduce] run$run%4d check!")
      }
    }


    'mapReduce {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128) + 1
        val L = IndexedSeq.tabulate(size, size) {
          (i,j) => if( i >= j ) rng.nextInt(100) else 0L
        }
        val mat = lMat.tabulate(size){L(_)(_)}

        def mapper( x: Double, i: Int, j: Int )
          = x.toLong ensuring x == L(i)(j)
        def reducer( a: Long, b: Long ) = a+b

        assert( L.map(_.sum).sum == mat.mapReduce(mapper)(reducer) )

        println( f"[${lMat.getClass.getSimpleName}/mapReduce] run$run%4d check!")
      }
    }


    'diag {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128) + 1
        val L = IndexedSeq.tabulate(size, size) {
          (i,j) => if( i >= j ) rng.nextInt(200) - 100 else 0
        }
        val diag = lMat.tabulate(size){L(_)(_)}.diag()
        for( i <- 0 until size )
          assert( diag(i) == L(i)(i) )
        println( f"[${lMat.getClass.getSimpleName}/diag] run$run%4d check!")
      }
    }


    'diagForeach {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128) + 1
        val L = IndexedSeq.tabulate(size, size) {
          (i,j) => if( i >= j ) rng.nextInt(200) - 100 else 0
        }
        val mat = lMat.tabulate(size){L(_)(_)}
        var sum = 0
        var i = 0
        mat.diagForeach{
          x =>
            sum += x.toInt
            i += 1
        }
        assert( sum == (0 until size).map{i => L(i)(i)}.sum )
        assert( i == size )
        println( f"[${lMat.getClass.getSimpleName}/diagForeach] run$run%4d check!")
      }
    }


    'diagForeachIndexed {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128) + 1
        val L = IndexedSeq.tabulate(size, size) {
          (i,j) => if( i >= j ) rng.nextDouble()*200 - 100 else 0
        }
        val mat = lMat.tabulate(size){L(_)(_)}
        val visited = Array.fill(size)(false)
        mat.diagForeach{
          (x,i) =>
            assert( ! visited(i) )
            assert( x == L(i)(i) )
            visited(i) = true
        }
        assert( visited.forall{_ == true} )
        println( f"[${lMat.getClass.getSimpleName}/diagForeach] run$run%4d check!")
      }
    }


    'triSolve {
      val rng = new Random(1337)
      def test( mat: MatFactory )
        = for( _ <- 0 until 128 )
          {
            val size = rng.nextInt(64)+1
            val nCols= rng.nextInt(128)+1
            val L = lMat.tabulate(size){
              (i,j) =>
                if(i==j) (1+9*rng.nextDouble)*{1-2*rng.nextInt(2)}
                else     2-4*rng.nextDouble
            }
            val rhs = mat.tabulate(size,nCols)( (_,_) => 10-20*rng.nextDouble )
            val solution = rhs.copy
            L.triSolve(solution)
            val k = 0 until size
            val check = Mat.tabulate(rhs.nRows,rhs.nCols){ (i,j) => k.map( k => L(i,k)*solution(k,j) ).sum }
            check.foreach{ (x,i,j) => assert{ isClose( x, rhs(i,j) ) } }
          }
      test(MatRM)
      test(MatCM)
    }


    'triInvert {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128)+1
        val L = Array.tabulate(size){
          i =>    Vec.tabulate(size){
          j => if( i > j )
            rng.nextDouble()*200-100
          else if( i == j )
            rng.nextDouble()*950+ 50
          else
            0
        }}
        val mat = lMat.tabulate(size)( L(_)(_) )
        mat.triInvert()
        val L_inv = Array.tabulate(size)( i => Vec.tabulate(size){mat(_,i)} )

        val I = Array.tabulate(size,size){
          (i,j) => L{i} ⋅ L_inv{j}
        }

        for( row <- 0 until size;
             col <- 0 until size )
          if( row == col ) assert( 1e-6 > {I(row)(col) - 1}.abs )
          else             assert( 1e-6 > {I(row)(col)    }.abs )
        println( f"[${lMat.getClass.getSimpleName}/invert] run$run%4d check!")
      }
    }


    'choleskyDecomp {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128) + 1
        val L = Array.tabulate(size){
          i =>    Vec.tabulate(size){
          j => if( i > j )
            rng.nextDouble() * 200 - 100
          else if( i == j )
            rng.nextDouble() * 950 + 50
          else
            0
        }}
        val LLT = Array.tabulate(size){ i => Vec.tabulate(size){ L(i) ⋅ L(_) }}
        val mat = lMat.tabulate(size){ LLT(_)(_) }
        mat.choleskyDecomp()
        for (
          i <- 0 until size;
          j <- 0 to i
        )
          assert( 1e-6 >= { mat(i, j) - L(i)(j) }.abs )
        println( f"[${lMat.getClass.getSimpleName}/choleskyDecomp] run$run%4d check!")
      }
    }


    'choleskySolve {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128)+1
        val L = Array.tabulate(size){
          i =>    Vec.tabulate(size){
          j => if( i > j )
            rng.nextDouble()*200-100
          else if( i == j )
            rng.nextDouble()*950+ 50
          else
            0
        }}
        val LLT = Array.tabulate(size){ i => Vec.tabulate(size){ L(i) ⋅ L(_) }}
        val mat = lMat.tabulate(size)( LLT(_)(_) )
        mat.choleskyDecomp()
        for( _ <- 0 until 128 )
        {
          val u  = Vec.tabulate(size)( _ => rng.nextDouble()*200 - 1 )
          val v = u.clone
          mat.choleskySolve(v)
          for( r <- 0 until size )
          {
            val u_r = v ⋅ LLT{r}
            assert( {u_r - u(r)}.abs <= 1e-6 )
          }
        }
        println( f"[${lMat.getClass.getSimpleName}/choleskySolve] run$run%4d check!")
      }
    }


    'choleskySolveMat {
      val rng = new Random(1337)
      def test( mat: MatFactory )
        = for( run <- 0 until 128 )
          {
            val size = rng.nextInt(64)+1
            val nCols= rng.nextInt(128)+1
            var A = lMat.tabulate(size){
              (i,j) =>
                if(i==j) 0.5 + 9.5*rng.nextDouble
                else     1.0 - 2.0*rng.nextDouble
            }
            val k = 0 until size
            A = lMat.tabulate(size){ (i,j) => k.map( k => A(i,k)*A(j,k) ).sum }
            val L = A.copy
            L.choleskyDecomp()
            val rhs = mat.tabulate(size,nCols)( (_,_) => 10-20*rng.nextDouble )
            val solution = rhs.copy
            L.choleskySolve(solution)
            val check = Mat.tabulate(rhs.nRows,rhs.nCols){ (i,j) => k.map( k => A.symGet(i,k)*solution(k,j) ).sum }
            check.foreach{ (x,i,j) => assert{ isClose( x, rhs(i,j) ) } }
          }
      test(MatRM)
      test(MatCM)
    }


    'choleskyInvert {
      val rng = new Random(1337)

      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128)+1
        val L = Array.tabulate(size){
          i =>    Vec.tabulate(size){
          j => if( i > j )
            rng.nextDouble()*200-100
          else if( i == j )
            rng.nextDouble()*950+ 50
          else
            0
        }}
        val LLT = Array.tabulate(size){ i => Vec.tabulate(size){ L(i) ⋅ L(_) }}
        val mat = lMat.tabulate(size)( LLT(_)(_) )
        mat.choleskyDecomp()
        mat.choleskyInvert()
        val LLT_inv = Array.tabulate(size){
          i =>          Vec.tabulate(size){
          j => if( j > i ) mat(j,i) else mat(i,j)
        }}
        val I = Array.tabulate(size,size){ LLT(_) ⋅ LLT_inv(_) }
        for( row <- 0 until size;
             col <- 0 until size )
          if( row == col ) assert( 1e-6 > {I(row)(col) - 1}.abs )
          else             assert( 1e-6 > {I(row)(col)    }.abs )

        println( f"[${lMat.getClass.getSimpleName}/choleskyInvert] run$run%4d check!")
      }
    }


    'symRowForeach {
      val rng = new Random(1337)

      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128)+1

        val A = collection.mutable.IndexedSeq.tabulate(size,size){
          (i,j) => if( i >= j ) rng.nextDouble()*200-100 else 0
        }

        for( row <- 0 until size;
             col <- 0 until row )
          A(col)(row) = A(row)(col)

        for( row <- 0 until size;
             col <- 0 until size )
          assert( A(row)(col) == A(col)(row) )

        val a = lMat.tabulate(size){A(_)(_)}

        for( row <- 0 until size )
        {
          val visited = Array.fill(size)(false)
          a.symRowForeach(row){
            (x,col) =>
              assert( ! visited(col) )
              assert( x == A(row)(col) )
              visited(col) = true
          }
          assert( visited.forall(_ == true) )
        }

        println( f"[${lMat.getClass.getSimpleName}/symRowForeach] run$run%4d check!")
      }
    }


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

        val a = lMat.tabulate(size){A(_)(_)}

        val trRef = (0 until size).map{ i => A(i) ⋅ B(i) }.sum
        val tr = a symMatProdTrace {B(_)(_)}

        assert{ 1e-6 > (tr - trRef).abs }

        println( f"[${lMat.getClass.getSimpleName}/symMatProdTrace(_,_)] run$run%4d check!")
      }
    }


    'symSquare {
      val rng = new Random(1337)
      for( run <- 0 until 256 )
      {
        val size = rng.nextInt(128)+1

        val A = Array.tabulate(size){ i => Vec.tabulate(size){ j => if( i >= j ) rng.nextDouble()*200-100 else 0 }}
        for( row <- 0 until size;
             col <- 0 until row )
        {
          A(col)(row) = A(row)(col)
        }
        for( row <- 0 until size;
             col <- 0 until size )
        {
          assert( A(row)(col) == A(col)(row) )
        }

        val a = lMat.tabulate(size){A(_)(_)}
        val aa = a.symSquare

        for( row <- 0 until size;
             col <- 0 to row )
        {
          val result = aa(row,col)
          val ref = A{row} ⋅ A{col}
          assert( 1e-7 > { result - ref }.abs )
        }

        println( f"[${lMat.getClass.getSimpleName}/symSquare] run$run%4d check!")
      }
    }


    'symSquareLike {
      val rng = new Random(1337)
      for( run <- 0 until 32 )
      {
        val size = rng.nextInt(128)+1

        val A = Array.tabulate(size){ i => Vec.tabulate(size){ j => if( i >= j ) rng.nextDouble()*200-100 else 0 }}
        for( row <- 0 until size;
             col <- 0 until row )
        {
          A(col)(row) = A(row)(col)
        }
        for( row <- 0 until size;
             col <- 0 until size )
        {
          assert( A(row)(col) == A(col)(row) )
        }

        val a = lMat.tabulate(size){A(_)(_)}
        val aa = a symSquareLike {
          (a_ik, a_kj, i,k,j) =>
            assert( a_ik == A(i)(k) )
            assert( a_kj == A(k)(j) )
            a_ik * a_kj
        }

        for( row <- 0 until size;
             col <- 0 to row )
        {
          val result = aa(row,col)
          val ref = A{row} ⋅ A{col}
          assert( 1e-7 > { result - ref }.abs )
        }

        println( f"[${lMat.getClass.getSimpleName}/symSquareLike] run$run%4d check!")
      }
    }


    'symInner {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128)+1
        val A = lMat.tabulate(size)( (_,_) => rng.nextDouble*64 - 32 )

        for( _ <- 0 until 128 )
        {
          val u = Vec.tabulate(size)( _ => rng.nextDouble*64 - 32 )
          val v = Vec.tabulate(size)( _ => rng.nextDouble*64 - 32 )

          val uAv = LMat.symInner(u,A,v)

          val check = (0 until size).map(
            i => (0 until size).map(
              j => u(i) * A.symGet(i,j) * v(j)
            ).sum
          ).sum

          assert{ isClose(uAv,check) }
        }

        println( f"[${lMat.getClass.getSimpleName}/symInner] run${run}%4d check!")
      }
    }


    'symMatMul {
      val rng = new Random(1337)
      for( run <- 0 until 128 )
      {
        val size = rng.nextInt(128)+1
        val A = lMat.tabulate(size)( (_,_) => rng.nextDouble*64 - 32 )

        for( _ <- 0 until 128 )
        {
          val u = Vec.tabulate(size)( _ => rng.nextDouble*64 - 32 )
          val Au = A.symMatMul(u)

          val check = Vec.tabulate(size) {
            i => (0 until size).map{ j => A.symGet(i,j) * u(j) }.sum
          }

          assert{ allClose(Au,check) }
        }

        println( f"[${lMat.getClass.getSimpleName}/symInner] run${run}%4d check!")
      }
    }


    'unaryMinus {
      val rng = new Random(1337)
      for( run <- 0 until 1024 )
      {
        val size = rng.nextInt(128)+1
        val A = lMat.tabulate(size)( (_,_) => rng.nextDouble*64 - 32 )
        val negA = -A

        A foreach {
          (Aij,i,j) => assert{ -Aij == negA(i,j) }
        }
      }
    }
  }
}