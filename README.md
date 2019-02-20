# Overview
GauProS is a lightweight, memory-friendly implementation of Gaussian Processes for Scala and ScalaJS.
Aside from the Scala Standard Library there are no deployment dependencies, particularly no native
dependencies. An interactive example of `GauProS` running in the Browser can be found [here](
  https://dirktoewe.github.io/gaupros_sandbox/gaupros_sandbox-0.1.0.html#QxfmC0GxF/xDsMFRP/MzM0JrYVZCRryQQ/GDQkK0bpBDc8HQQs6U2UNlX/VDHjRsQy0nP0OXA35D3T8wQ6gCx0P2vh5DualsQ/Db50P1zkZDnnzoRAhN9UL+yH5EG+qsQtosf0Qp+NpDrDdoRCyWR0PAe3pEMdEjQ8w/50RJBlBEBS1qREJ8vkQOqBhEbkmrQ4b0DURqCdlDlf1DRHphxkO2BcNEhYaxQ50uMUSKl7VDiOof
).

# GP without Hyperparameters
The Kernel function is a functional trait of type `(xi: A, xj: A) => Double`. If the Kernel only has fixed/known
hyperparameters, it can simply be instantiated using a lambda expression.

```scala
import gps.regression._

// SOURCE: https://data.bls.gov/timeseries/APU0000711211
val banana_prices = Vec(
  0.629, 0.641, 0.634, 0.629, 0.622, 0.617, 0.616, 0.611, 0.605, 0.598, 0.561, 0.571, 0.586, 0.587, 0.575, 0.580, 0.571, 0.577, 0.583, 0.576, 0.573, 0.580, 0.581, 0.587,
  0.596, 0.625, 0.621, 0.621, 0.617, 0.614, 0.611, 0.606, 0.607, 0.607, 0.598, 0.599, 0.604, 0.603, 0.607, 0.603, 0.599, 0.605, 0.604, 0.595, 0.597, 0.602, 0.600, 0.606,
  0.609, 0.611, 0.609, 0.597, 0.603, 0.603, 0.602, 0.595, 0.601, 0.586, 0.587, 0.592, 0.595, 0.599, 0.597, 0.597, 0.603, 0.607, 0.606, 0.608, 0.606, 0.582, 0.592, 0.585,
  0.583, 0.591, 0.593, 0.597, 0.582, 0.574, 0.581, 0.580, 0.580, 0.577, 0.575, 0.580, 0.581, 0.573, 0.586, 0.574, 0.570, 0.569, 0.567, 0.562, 0.575, 0.570, 0.573, 0.576,
  0.573, 0.573, 0.571, 0.571, 0.564, 0.572, 0.565, 0.561, 0.544, 0.548, 0.551, 0.558, 0.568, 0.574, 0.582, 0.578, 0.575, 0.575, 0.575, 0.569, 0.573, 0.574, 0.566, 0.577,
  0.576
)
val months = Array.tabulate(banana_prices.length)(-banana_prices.length + 1.0 + _)

val σ_f = 1.0 / 64
val l   = 1.0
val mean_price = banana_prices.sum / banana_prices.length
val kernel: Kernel[Double] = (xi,xj) => σ_f*σ_f * Math.exp( -(xi-xj)*(xi-xj) / (2*l*l) )

val gpr = GPR( x=months, y=banana_prices, y_shift=mean_price, kernel=kernel )

gps.regression.gpr.plot1d(gpr, Vec(months: _*), banana_prices, marginAbs=0.5, marginRel=0.0)
```

# GP Regression with Hyperparameter Optimization
`GauProS` allows the easy composition and training of kernels with hyperparameters.
Trainable parameters are identified by a `Symbol`. `GauProS` uses an (un)boxed variant
of a BFGS optimizer. Said optimizer does not yet support starting from the box borders.
The initial hyperparameter values have to be well inside and not on the boundary of
the box.

```scala
// http://www.gaussianprocess.org/gpml/chapters/RW2.pdf#page=13
val kernel = Noise('σ_n * 'σ_n)  +  'σ_f * 'σ_f * Exp( - AbsDelta.pow(2.0) / (2 * 'l * 'l) )

val gpr = GPR.fit_shifted(
  x=months       ,
  y=banana_prices,
  y_shift = 'y_shift,
  kernel=kernel,
  likelihood_fn = GPR.logp_loo[Double], // <- The loss/likelihood function used for optimization.
                                        //    Built-in: GPR.logp_marginal and GPR.logp_loo
  param_init = Map(
    'σ_n -> 1,
    'σ_f -> 1,
    'l -> 4
  ),
  param_min = Map.empty[Symbol,Double] withDefaultValue Double.NegativeInfinity,
  param_max = Map.empty[Symbol,Double] withDefaultValue Double.PositiveInfinity
)
println(gpr.params) // <- the trained hyperparameter values

gps.regression.gpr.plot1d(gpr, Vec(months: _*), banana_prices, marginAbs=0.5, marginRel=0.0)
```
