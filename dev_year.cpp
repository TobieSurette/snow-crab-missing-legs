#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator()(){
   using namespace density;

   // Data declarations:
   DATA_MATRIX(z);                // Missing leg indicator observations.
   DATA_IVECTOR(maturity);        // Crab maturity.
   DATA_IVECTOR(size);            // Crab size.
   DATA_IVECTOR(year);            // Survey year.

   // Parameter declarations:
   PARAMETER(alpha);              // Intercept parameters.
   PARAMETER_VECTOR(leg_effect);  // Leg random effect.
   PARAMETER(log_sigma_leg);      // Leg random effect error.
   PARAMETER(beta_maturity);      // Maturity parameter.
   PARAMETER_VECTOR(size_effect); // Size random effect.
   PARAMETER_VECTOR(beta_size);   // Size scale parameters.

   PARAMETER(xp_size);            // Size location parameter.
   PARAMETER_VECTOR(xp_size_year);            // Size location parameter.
   PARAMETER(log_sigma_xp_size_year);         // Size location error parameter.

   PARAMETER(logit_rho_size);     // Size effect correlation parameter.
   PARAMETER_VECTOR(crab_effect); // Crab random effect.
   PARAMETER(log_sigma_crab);     // Crab random effect error.
   PARAMETER_VECTOR(year_effect); // Year effect.
   PARAMETER(log_sigma_year);     // Year effect error.

   PARAMETER_MATRIX(year_size_effect); // Year x size random effect.
   PARAMETER(log_sigma_year_size);     // Year x size effect error.

   // Variable dimensions:
   int n_obs  = z.rows();                        // Number of observations.
   int n_legs = z.cols();                        // Number of legs per crab.
   int n_year = year_effect.size();              // Number of survey years.

   // Transformed
   Type res = 0;                                       // Log-likelihood accumulator.
   Type rho_size = 1.0 / (1.0 + exp(-logit_rho_size)); // Size effect correlation parameter.

   // Random effects:
   res -= sum(dnorm(leg_effect, Type(0), exp(log_sigma_leg), true));   // Leg random effect.
   res += AR1(rho_size)(size_effect);                                  // Size random effect.
   res -= sum(dnorm(crab_effect, Type(0), exp(log_sigma_crab), true)); // Crab random effect.
   res -= sum(dnorm(year_effect, Type(0), exp(log_sigma_year), true)); // Year random effect.

   res -= sum(dnorm(xp_size_year, Type(0), exp(log_sigma_xp_size_year), true));
   res -= sum(dnorm(year_size_effect, Type(0), exp(log_sigma_year_size), true)); // Year random effect.

   // Likelihood model:
   matrix<Type> p(n_obs,n_legs);
   for (int i = 0; i < n_obs; i++){
      for (int j = 0; j < n_legs; j++){
         // Logit-linear binomial mean:
         Type mu = alpha +
                   leg_effect[j] +
                   beta_maturity * maturity[i] +
                   beta_size[0] * (size[i] - xp_size) * (size[i] - xp_size) +
                   beta_size[1] * size_effect[size[i]] +
                   year_effect[year[i]] +
                   year_size_effect(year[i],size[i]) +
                   crab_effect[i];

         // Calculate binomial probability:
         p(i,j) = 1.0 / (1.0 + exp(-mu));

         // Binomial likelihood:
         res -= z(i,j) * log(p(i,j)) + (1-z(i,j)) * log(1-p(i,j));
      }
   }

   // Export variables:
   REPORT(p);

   return res;
}

