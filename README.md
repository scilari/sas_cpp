# Systematic Alias Sampling C++ implementation

The Alias Method by Walker offers an elegant way to sample from a discrete distribution in constant time.
A good explanation of the method can be found [here](http://www.keithschwarz.com/darts-dice-coins). 
The code here implements a systematic version of the Alias Method described in paper 

Vallivaara et al. 
*"Systematic Alias Sampling: an efficient and low-variance way to sample from a discrete distribution"*

ACM Transactions on Mathematical Software (TOMS)
Volume 43 Issue 3, August 2016
Article No. 18

http://dl.acm.org/citation.cfm?id=2935745

The Scala implementation of the method achieves 5-20X speed up compared to Apache Commons Math NormalDistribution.sample() when sampling in batches. 
The empirical distribution of the batches also have significantly better goodness-of-fit according to Cramer-Von-Mises statistic.

## Notes:
* This is a very rudimentary implementation, not including e.g. the Golden Ratio version of the method.
* Only basic properties are tested.
* The more comprehensive Scala implementation can be found [here](https://github.com/scilari/sas).
* The C++ performance is not tested here. However, initial performance tests for a C++ implementation by Bolong Zhang
 suggest that the performance gains are similar to the Scala version. 
 The tests can be found [here](https://github.com/bolongz/Systematic-Alias-Sampling). 
