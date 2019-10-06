# Actually running the model

## Forward model

## Inverse model

This takes a bit of care.

- only run it when the forward model runs flawlessly
- think about your priors
- first do MCMC inference of the logprior only and check its
  convergence.  To speed up convergence
  - decrease the number of parameters
- only then run the MCMC over the logposterior.  Tipically use about
  3x (?) as many samples.
