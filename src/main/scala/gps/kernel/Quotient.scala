package gps.kernel

/**
  * Created by dtitx on 31.07.17.
  */
private[kernel] object Quotient
{
  def apply[X]( dividend: Kernel[X], divisor: Kernel[X] ): Kernel[X]
    = dividend * (divisor pow -1)
}