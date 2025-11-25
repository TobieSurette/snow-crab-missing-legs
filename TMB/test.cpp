// Simple linear regression.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() (){
   // Data section:
   DATA_MATRIX(z);                         // Missing leg matrix pattern.
   //DATA_FACTOR(width);                   // Crab carapace width.
   //DATA_FACTOR(maturity);                // Crab maturity.
   //DATA_FACTOR(year);                    // Survey year.

   // Parameter section:
   PARAMETER_VECTOR(leg_effect);           // Leg effect parameters.
   PARAMETER(log_sigma_leg);               // Leg effect error parameter.

   PARAMETER_VECTOR(L_eps);                // L-matrix for correlated error structure.
   PARAMETER_VECTOR(log_scale_eps);        // Scale parameters for error covariance structure.
   PARAMETER(log_sigma_L_eps);
   PARAMETER(log_sigma_scale_eps);

   // Parameter variable declarations:
   //PARAMETER_VECTOR(width_effect);          // Crab size effect effect.
   //PARAMETER_VECTOR(year_effect);           // Survey year effect.
   //PARAMETER(maturity_effect);              // Maturity effect.
   //PARAMETER(log_sigma_width);              // Error parameter for length effects.
   //PARAMETER(log_sigma_year);               // Error parameter for year effects.

   // Initialize or transform variables:
   Type nll = 0;                           // Initialize negative log-likelihood variable.
   int n_obs  = z.rows();                  // Number of observations.
   int n_legs = z.cols();                  // Number of dimensions.
   // int n_width = width_effect.size();   // Number carapace width categories.
   // int n_year = year_effect.size();     // Number carapace width categories.

   // Random effect likelihood contributions:
   nll = -sum(dnorm(leg_effect, Type(0), exp(log_sigma_leg), true));
   nll = -sum(dnorm(L_eps, Type(0), exp(log_sigma_L_eps), true));
   nll = -sum(dnorm(log_scale_eps, Type(0), exp(log_sigma_scale_eps), true));

   // nll = -sum(dnorm(width_effect, Type(0), exp(log_sigma_width), true));
   // nll = -sum(dnorm(year_effect, Type(0), exp(log_sigma_year), true));

   // Define unstructured covariance matrix for between-leg correlations:
   int k = 0;
   matrix<Type> L(n_legs,n_legs);
   matrix<Type> D(n_legs,n_legs);
   L.fill(0); D.fill(0);
   for (int i = 0; i < n_legs; i++) L(i,i) = 1.0;
   for (int i = 1; i < n_legs; i++){
      for (int j = 0; j < i; j++){
         L(i,j) = L_eps[k]; k++;
      }
   }
   for (int i = 0; i < n_legs; i++){
      D(i,i) = sqrt(exp(log_scale_eps[i]));
   }
   matrix<Type> Sigma_eps = (D * (L * L.transpose())) * D;
   using namespace density;
   MVNORM_t<Type> density_eps(Sigma_eps);

   // Calculate log-likelihood:
   vector<Type> eta(n_legs);     // Logit-linear mean.
   vector<Type> p(n_legs);       // Binomial probabilities.
   matrix<Type> u(n_obs,n_legs); // Binomial log-likelihood.
   for (int i = 0; i < n_obs; i++){
      nll += density_eps(u.row(i));                                // Correlated errors.
      for (int j = 0; j < n_legs; j++){
         eta[j] = leg_effect[j] + u(i,j);                          // Logit-linear mean.
         p[j]= 1.0 / (1.0 + exp(-eta[j]));                         // Binomial probabilities.
         nll -= (1.0-z(i,j)) * log(1.0-p[j]) + z(i,j) * log(p[j]); // Binomial log-likelihood.
      }
   }

   // ADREPORT(exp(2*logSigma));

   return nll;
}

