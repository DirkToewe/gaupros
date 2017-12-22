package gps.kernel

/**
  * Created by dtitx on 31.07.17.
  */
private[kernel] object Difference
{
  def apply[X]( minuend: Kernel[X], subtrahend: Kernel[X] ): Kernel[X]
    = minuend + -subtrahend
}