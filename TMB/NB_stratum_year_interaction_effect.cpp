// Negative binomial example:
// Parameterized by mean 'mu' and precision parameter 'r':
// Stratum effect is included.
// Mean : 'mu = stratum_effect + year_effect + stratum_year_effect
// Variance : mu + (mu^2)/r

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Variable declarations:
  DATA_VECTOR(y);                        // Vector of count observations.
  DATA_FACTOR(stratum);                  // Stratum factor level.
  DATA_FACTOR(year);                     // Year factor level.
  PARAMETER(alpha);                      // Log-scale global mean parameter.
  PARAMETER_VECTOR(stratum_effect);      // stratum effect.
  PARAMETER(log_sigma_stratum);          // Log-scale stratum effect standard error.
  PARAMETER_VECTOR(year_effect);         // Year effect.
  PARAMETER(log_sigma_year);             // Log-scale year effect standard error.
  PARAMETER_VECTOR(stratum_year_effect); // Stratum-year interaction effect.
  PARAMETER(log_sigma_stratum_year);     // Log-scale stratum-year effect standard error.
  PARAMETER(log_r);                      // Log-scale precision parameter.

  // Transformed parameters:
  Type r = exp(log_r);
  Type sigma_stratum = exp(log_sigma_stratum);
  Type sigma_year = exp(log_sigma_year);
  Type sigma_stratum_year = exp(log_sigma_stratum_year);

  // Number of observations and factor levels:
  int nobs = y.size();
  int nstratum = stratum_effect.size();
  int nyear = year_effect.size();

  // Negative log-likelihood accumulator:
  Type res = 0;
  
  // Prior over stratum effects:
  for (int j=0; j < nstratum; j++){
     res -= dnorm(stratum_effect[j], Type(0.0), sigma_stratum, true);
  }

  // Prior over year effects:
  for (int j=0; j < nyear; j++){
     res -= dnorm(year_effect[j], Type(0.0), sigma_year, true);
  }

  // Prior over stratum-year interaction effects:
  for(int i=0; i < nstratum; i++){
     for(int j=0; j < nyear; j++){
        res -= dnorm(stratum_year_effect[i*nyear+j], Type(0.0), sigma_stratum_year, true);
     }
  }
  
  // Likelihood evaluation:
  Type mu;
  for(int i = 0; i < nobs; i++){
     mu = exp(alpha + stratum_effect[stratum[i]-1] + year_effect[year[i]-1] + stratum_year_effect[(stratum[i]-1)*nyear+(year[i]-1)]);
     res -= lgamma(y[i]+r) - lgamma(r) - lgamma(y[i]+1) + r*log(r) + y[i]*log(mu) - (r+y[i])*log(r+mu);
  }
  return res;
}
