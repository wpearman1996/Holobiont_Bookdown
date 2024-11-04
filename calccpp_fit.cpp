#include <Rcpp.h>
using namespace Rcpp;

// Fitness function
double fitness_func_cpp(double selection_parameter, double optima, double trait) {
  double q = exp(pow((trait - optima), 2) / -selection_parameter);
  return q;
}

// [[Rcpp::export]]
NumericVector calculate_host_fitness_cpp(NumericMatrix HostPopulation, NumericVector traitpool_microbes, NumericVector microbiome_importances, NumericVector host_preference_optima, double selection_parameter_hosts, double env_cond_val_used) {
  
  int n = HostPopulation.ncol();
  NumericVector hostfitness(n);

  for(int host = 0; host < n; host++) {
    NumericVector pop = HostPopulation(_, host);
    double mean_microbial_trait_val = sum(traitpool_microbes * pop) / sum(pop);
    double composite_host_trait = ((mean_microbial_trait_val * microbiome_importances[host]) + (host_preference_optima[host] * (1 - microbiome_importances[host])));
    
    hostfitness[host] = fitness_func_cpp(selection_parameter_hosts, env_cond_val_used, composite_host_trait);
  }
  
  return hostfitness;
}

