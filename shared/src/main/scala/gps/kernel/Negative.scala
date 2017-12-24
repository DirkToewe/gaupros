package gps.kernel

/**
  *
  * Created by Dirk Toewe on 31.07.17.
  */
private[kernel] object Negative
{
  def apply[X]( kernel: Kernel[X] ): Kernel[X]
    = -1 * kernel
}