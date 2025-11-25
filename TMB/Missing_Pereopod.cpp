// Single-year binomial regression analysis of missing crab pereopods: 
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() (){
   // Data variable declarations:
   DATA_VECTOR(z);                          // Indicator vector of missing pereopods. 
   DATA_FACTOR(width);                      // Vector of observed carapace widths IDs.
   DATA_VECTOR(parity);                     // Vector indicating the side of the crab. 
   DATA_FACTOR(position);                   // Vector indicating the position of the pereopod.
   DATA_VECTOR(maturity);                   // Vector of crab maturity IDs.
   DATA_FACTOR(carapace);                   // Vector of carapace condition IDs.
   DATA_FACTOR(grid);                       // Vector of sampling grid IDs.
 
   // Parameter variable declarations:
   PARAMETER(log_alpha);                    // Log-scale average missing pereopod rate.
   PARAMETER_VECTOR(width_effect);          // Crab size effect effect.
   PARAMETER_VECTOR(position_effect);       // Pereopod position effect.
   PARAMETER(parity_effect);                // Crab side effect.
   PARAMETER(maturity_effect);              // Maturity effect.
   PARAMETER_VECTOR(carapace_effect);       // Carapace condition effect.
   PARAMETER_VECTOR(grid_effect);           // Sampling grid identifier.
   PARAMETER(log_sigma_width);              // Log-scale error parameter for length effects.
   PARAMETER(log_sigma_position);           // Log-scale error parameter for pereopod position effects.
   PARAMETER(log_sigma_carapace);           // Log-scale error parameter for carapace condition effects.
   PARAMETER(log_sigma_grid);               // Log-scale error parameter for grid effects.
  
   // Interaction effects:
   PARAMETER_VECTOR(position_maturity_effect); // Position x maturity interaction term.
   PARAMETER(log_sigma_position_maturity);     // Log-scale error parameter for position x maturity effects.
   PARAMETER_VECTOR(width_position_effect); // Carapace width x leg-position interaction term.
   PARAMETER(log_sigma_width_position);     // Log-scale error parameter for carapace width x position effects.
   PARAMETER_VECTOR(width_maturity_effect); // Carapace width x maturity interaction term.
   PARAMETER(log_sigma_width_maturity);     // Log-scale error parameter for carapace width x maturity effects. 
   PARAMETER_VECTOR(grid_maturity_effect);  // Carapace grid x maturity interaction term.
   PARAMETER(log_sigma_grid_maturity);      // Log-scale error parameter for grid x maturity effects.
   PARAMETER_VECTOR(grid_position_effect);  // Carapace grid x position interaction term.
   PARAMETER(log_sigma_grid_position);      // Log-scale error parameter for grid x position effects.
      
   // Calculated variables:
   Type res = 0;                            // Negative log-likelihood accumulator.
   int n_obs = z.size();                    // Number of length-frequency categories. 
   int n_width = width_effect.size();       // Number carapace width categories. 
   int n_position = position_effect.size(); // Number of crab legs per side. 
   int n_carapace = carapace_effect.size(); // Number of carapace condition categories. 
   int n_grid = grid_effect.size();         // Number of sampling grids.
      
   // Prior over carapace width effects:
   Type sigma_width = exp(log_sigma_width); 
   for (int j = 0; j < n_width; j++){
      res -= dnorm(width_effect[j], Type(0.0), sigma_width, true);
   }
  
   // Prior over position effects:
   Type sigma_position = exp(log_sigma_position); 
   for (int j = 0; j < n_position; j++){
      res -= dnorm(position_effect[j], Type(0.0), sigma_position, true);
   }
   
   // Prior over carapace condition effects:
   Type sigma_carapace = exp(log_sigma_carapace); 
   for (int j = 0; j < n_carapace; j++){
      res -= dnorm(carapace_effect[j], Type(0.0), sigma_carapace, true);
   }
   
   // Prior over grid effects:
   Type sigma_grid = exp(log_sigma_grid);  
   for (int j = 0; j < n_grid; j++){
      res -= dnorm(grid_effect[j], Type(0.0), sigma_grid, true);
   }
   
   // Position x Maturity interaction effects:
   Type sigma_position_maturity = exp(log_sigma_position_maturity); 
   for (int j = 0; j < n_position; j++){
      res -= dnorm(position_maturity_effect[j], Type(0.0), sigma_position_maturity, true);
   }
   
   // Width x Position interaction effects:
   Type sigma_width_position = exp(log_sigma_width_position); 
   for (int i = 0; i < n_width; i++){
      for (int j = 0; j < n_position; j++){
         res -= dnorm(width_position_effect[i * n_position + j], Type(0.0), sigma_width_position, true);
      }
   }
    
   // Prior over carapace width x maturity effects:
   Type sigma_width_maturity = exp(log_sigma_width_maturity); 
   for (int j = 0; j < n_width; j++){
      res -= dnorm(width_maturity_effect[j], Type(0.0), sigma_width_maturity, true);
   }
   
   // Prior over carapace grid x maturity effects:
   Type sigma_grid_maturity = exp(log_sigma_grid_maturity); 
   for (int j = 0; j < n_grid; j++){
      res -= dnorm(grid_maturity_effect[j], Type(0.0), sigma_grid_maturity, true);
   }
   
   // Grid x Position interaction effects:
   Type sigma_grid_position = exp(log_sigma_grid_position); 
   for (int i = 0; i < n_grid; i++){
      for (int j = 0; j < n_position; j++){
         res -= dnorm(grid_position_effect[i * n_position + j], Type(0.0), sigma_grid_position, true);
      }
   }
   
   // Likelihood functions:
   Type mu = 0;
   Type p = 0;
   for (int i = 0; i < n_obs; i++){  
      // Logit-linear prediction:
      mu = log_alpha + 
           width_effect[width[i]-1] +  
           position_effect[position[i]-1] + 
           carapace_effect[carapace[i]-1] + 
           grid_effect[grid[i]-1] +
           parity_effect * parity[i] + 
           maturity_effect * maturity[i] + 
           position_maturity_effect[position[i]-1] * maturity[i] +
           width_maturity_effect[width[i]-1] * maturity[i] + 
           grid_maturity_effect[grid[i]-1] * maturity[i] + 
           width_position_effect[(width[i]-1) * n_position + (position[i]-1)] + 
           grid_position_effect[(grid[i]-1) * n_position + (position[i]-1)];  
          
      // Convert to probability:
      p = 1 / (1 + exp(-mu));
      
      // Bernouilli log-likelihood:
      res -= (1-z[i]) * log(1-p) + z[i] * log(p);
   }   
  
   return res;
}  
  