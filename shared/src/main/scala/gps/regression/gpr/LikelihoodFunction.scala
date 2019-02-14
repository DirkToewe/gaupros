package gps.regression.gpr

import gps.opt.ObjectiveWithGradient

/**
  * Created by Dirk Toewe on 20.08.17.
  */
abstract class LikelihoodFunction extends ObjectiveWithGradient
{
  def params: IndexedSeq[Symbol]
}
