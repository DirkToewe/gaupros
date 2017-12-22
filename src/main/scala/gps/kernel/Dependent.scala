package gps.kernel

/** Extractor object for Kernels that depend on the input x. For
  *
  * Created by Dirk Toewe on 23.08.17.
  */
private[kernel] object Dependent
{
  def unapply[X]( kernel: Kernel[X] )
    = if( kernel.argDependent ) Some(kernel) else None
}
