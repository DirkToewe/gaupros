package gps.kernel

/**
  * Created by Dirk Toewe on 23.08.17.
  */
private[kernel] object Independent
{
  def unapply[@specialized X]( kernel: Kernel[X] )
    = if( kernel.argDependent )
        None
      else
        Some(kernel.asInstanceOf[Kernel[Any]])
}