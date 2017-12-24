package gps.kernel

object Kernel_experiments
{
  def main( args: Array[String] ) =
  {
    val deltaPow = AbsDelta.pow('norm_p)
    val argSum = Noise('var_noise) + 'var_func * Exp(
      - VecIdx(0, deltaPow / 't0)
      - VecIdx(1, deltaPow / 't1)
      - VecIdx(2, deltaPow / 't2)
    )
    println(argSum)
    println(argSum.params)
    println(argSum pDiff 'var_noise)
    println()
    println(argSum pDiff 'var_func)
    println()
    println(argSum pDiff 't0)
    println()
    println(argSum pDiff 't1)
    println()
    println(argSum pDiff 't2)
    println()
    println(argSum pDiff 'norm_p)
    println()
    println(
      argSum pDiff 'norm_p subs (
        't0 -> 1.0, 't1 -> 1.0, 't2 -> 1.0, 'var_noise -> 1.0, 'var_func -> 1.0, 'norm_p -> 2.05
      ) apply (Vec(0,0,0),Vec(0,0,0))
    )

    println("\n\n========\n\n")

    val h = 1 / 'x / 'x
    println(h)
    println(h pDiff 'x)
    println(h pDiff 'x pDiff 'x)
    println(h pDiff 'x subs 'x -> 2)
    println( 1 / (2 * 'x) )

    println("\n\n========\n\n")

    val i = 'x*'y*3*'z*'x
    println(i)
    println(i.params)
    println(i pDiff 'x)
    println(i pDiff 'x pDiff 'y)
    println()
    println('x + 'y + 'x + 'z + 3*'z + 'z)

    println("\n\n========\n\n")

    println(math.log(0.0) + math.log(0.0))
    println()
    println( Sin('x) pDiff 'x subs 'x -> 3 )
    println( math.cos(3) )
  }
}
