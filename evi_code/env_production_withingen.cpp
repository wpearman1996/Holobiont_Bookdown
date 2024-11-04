#include <Rcpp.h>
#include <numeric> // for std::accumulate
using namespace Rcpp;


double calculate_mean_fitness(NumericVector fit_value, NumericVector host_counts, int envpoolsize) {
  // Calculate the product of corresponding elements of microbe_probs_death and host_column
  NumericVector product = fit_value * host_counts;
  
  // Calculate the sum of the products
  double sum_of_products = sum(product);
  
  // Calculate the mean fitness
  double mean_fitness = sum_of_products / envpoolsize;
  
  return mean_fitness;
}

// Fitness function

// [[Rcpp::export]]
double fitness_func_cpp(double selection_parameter, double optima, double trait) {
  double q = exp(pow((trait - optima), 2) / -selection_parameter);
  return q;
}

// [[Rcpp::export]]
NumericVector rmultinom_cpp(int n, NumericVector prob) {
  SEXP res = Rf_allocVector(INTSXP, prob.size());
  prob = prob / sum(prob); // Adjust probabilities to sum to 1
  rmultinom(n, prob.begin(), prob.size(), INTEGER(res));
  return NumericVector(res);
}

// [[Rcpp::export]]
NumericMatrix process_microbe_probs(double env_cond_val,NumericVector fixed_envpool, NumericMatrix HostPopulation, double N_Microbes, double envpoolsize, double selection_parameter_env, NumericVector XY, NumericVector traitpool_microbes, NumericVector env_used) {
  int N_Species = traitpool_microbes.length();
  
  NumericMatrix var_microbeprobs(N_Species, 5); // Matrix for microbe probabilities
  
  // Compute microbe probabilities based on environmental conditions
    for (int i = 0; i < N_Species; ++i) {
      var_microbeprobs(i, 0) = env_used[i] / N_Microbes;
    }
    for (int i = 0; i < N_Species; ++i) {
      var_microbeprobs(i, 2) = fitness_func_cpp(selection_parameter_env, env_cond_val, traitpool_microbes[i]);
    }
    for (int i = 0; i < var_microbeprobs.nrow(); ++i) {
      if (NumericVector::is_na(var_microbeprobs(i, 2))) {
        var_microbeprobs(i, 2) = 0; // Replace NA with 0
      }
    }
    
    for (int i = 0; i < N_Species; ++i) {
      if (NumericVector::is_na(var_microbeprobs(i, 0))) {
        Rcout << "NA detected in var_microbeprobs(" << i << ", 0)" << std::endl;
      }
      if (NumericVector::is_na(var_microbeprobs(i, 2))) {
        Rcout << "NA detected in var_microbeprobs(" << i << ", 2)" << std::endl;
      }
      if (NumericVector::is_na(var_microbeprobs(i, 1))) {
        Rcout << "NA detected in var_microbeprobs(" << i << ", 1)" << std::endl;
      }
    }
    var_microbeprobs(_, 1) = (var_microbeprobs(_, 2) * var_microbeprobs(_, 0)) / sum(var_microbeprobs(_, 2) * var_microbeprobs(_, 0));
    for (int i = 0; i < var_microbeprobs.nrow(); ++i) {
      if (NumericVector::is_na(var_microbeprobs(i, 1))) {
        var_microbeprobs(i, 1) = 0; // Replace NA with 0
      }
    }
      // Matrix for final microbe probabilities
  NumericMatrix microbe_probs(N_Species, 8);
  
  // Compute final microbe probabilities
  microbe_probs(_, 0) = rowSums(HostPopulation, true) / sum(HostPopulation); // Relative abundance of microbes in host population
  microbe_probs(_, 1) = fixed_envpool / N_Microbes; // Relative abundance of microbes in fixed environment
  microbe_probs(_, 2) = var_microbeprobs(_, 1); // Fitness based sampling probability of environmental microbes based on fitness and relative abundance
  microbe_probs(_, 3) = XY[2] * microbe_probs(_, 0); // Sampling probability of host microbes
  microbe_probs(_, 4) = (1 - XY[1] - XY[2]) * microbe_probs(_, 1); // fixed environmental sampling probability
  microbe_probs(_, 5) = XY[1] * microbe_probs(_, 2); // environmental sampling probability
  microbe_probs(_, 6) = rowSums(microbe_probs(_, Range(3, 5)));
  for (int i = 0; i < N_Species; ++i) {
    if (NumericVector::is_na(microbe_probs(i, 2))) {
      Rcout << "NA detected in microbe_probs(" << i << ", 2)" << std::endl;
    }
    if (NumericVector::is_na(microbe_probs(i, 3))) {
      Rcout << "NA detected in microbe_probs(" << i << ", 3)" << std::endl;
    }
    if (NumericVector::is_na(microbe_probs(i, 4))) {
      Rcout << "NA detected in microbe_probs(" << i << ", 4)" << std::endl;
    }
    if (NumericVector::is_na(microbe_probs(i, 5))) {
      Rcout << "NA detected in microbe_probs(" << i << ", 5)" << std::endl;
    }
  }
  
  // Replace NA values with 0
  for (int i = 0; i < microbe_probs.nrow(); ++i) {
    for (int j = 0; j < microbe_probs.ncol(); ++j) {
      if (NumericVector::is_na(microbe_probs(i, j))) {
        microbe_probs(i, j) = 0;
      }
    }
  }
  
  // Sampled counts
  microbe_probs(_, 7) = rmultinom_cpp(envpoolsize, microbe_probs(_, 6));
  
  return microbe_probs;
}
