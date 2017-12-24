package gps.opt

class LineSearchFailure( val α: Double, val msg: String ) extends IllegalArgumentException(msg) {}
object LineSearchFailure
{
  def unapply( lineSearchFailure: LineSearchFailure )
    = Some(lineSearchFailure.α)
}