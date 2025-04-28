#include <Rcpp.h>
#include <numeric> // for std::accumulate
using namespace Rcpp;

// [[Rcpp::export]]
double fitness_func_bacgen_cpp(double selection_parameter, double optima1, double optima2, double trait,double weighting) {
  // Calculate the mean of optima1 and optima2
  double meanval = (weighting * optima1) + ((1-weighting)*optima2);
  
  // Calculate fitness function
  double q = exp(pow((trait - meanval), 2) / -selection_parameter);
  
  return q;
}

// [[Rcpp::export]]
double calculate_mean_fitness(NumericVector microbe_probs_death_column, NumericVector host_column, int N_Microbes) {
  // Calculate the product of corresponding elements of microbe_probs_death and host_column
  NumericVector product = microbe_probs_death_column * host_column;
  
  // Calculate the sum of the products
  double sum_of_products = sum(product);
  
  // Calculate the mean fitness
  double mean_fitness = sum_of_products / N_Microbes;
  
  return mean_fitness;
}

// [[Rcpp::export]]
NumericVector rmultinom_cpp(int n, NumericVector prob) {
  SEXP res = Rf_allocVector(INTSXP, prob.size());
  prob = prob / sum(prob); // Adjust probabilities to sum to 1
  rmultinom(n, prob.begin(), prob.size(), INTEGER(res));
  return NumericVector(res);
}

// [[Rcpp::export]]
NumericMatrix process_host_cpp(NumericMatrix HostPopulation, double selection_parameter_microbes, NumericVector host_microbe_optima, double env_condition, NumericVector traitpool_microbes, double N_Microbes, double self_seed_prop, NumericVector env_used, double weighting) {
  int N_Species = traitpool_microbes.length();
  int N_Hosts = HostPopulation.ncol();
  
  NumericMatrix microbe_probs_death(N_Species, 7); // Create microbe_probs_death matrix
  
  for (int host = 0; host < N_Hosts; ++host) {
    NumericVector host_column = HostPopulation(_, host); 
    microbe_probs_death(_, 5) = host_column / N_Microbes; // rel abunds in parent
    microbe_probs_death(_, 6) = env_used / N_Microbes; // rel abunds in environment
    
    double pref_used = env_condition;
    double microbe_used = host_microbe_optima[host];
    
    microbe_probs_death(_, 0) = host_column;
    for (int i = 0; i < N_Species; ++i) {
        microbe_probs_death(i, 2) = fitness_func_bacgen_cpp(selection_parameter_microbes, microbe_used, pref_used, traitpool_microbes[i],weighting);
      
    //  produce fitnesses
    }
  //  double mean_fitness = calculate_mean_fitness(microbe_probs_death(_, 2), host_column, N_Microbes);
  //  microbe_probs_death(_, 2) = (microbe_probs_death(_, 2));// / mean_fitness);
    for (int i = 0; i < microbe_probs_death.nrow(); ++i) {
      if (NumericVector::is_na(microbe_probs_death(i, 2))) {
        microbe_probs_death(i, 2) = 0; // Replace NA with 0
      }
    }
    NumericVector prob = (microbe_probs_death(_, 2) * microbe_probs_death(_, 5)); // multiple fitness by rel abunds in parent
    NumericVector env_prob = (microbe_probs_death(_, 2) * microbe_probs_death(_, 6)); // multiply fitnesses by rel abunds in environment
    NumericVector replacements_from_prevgen = rmultinom_cpp(N_Microbes * self_seed_prop, prob);
    NumericVector replacements_from_env = rmultinom_cpp(N_Microbes * (1 - self_seed_prop), env_prob);
    NumericVector replacements = replacements_from_env + replacements_from_prevgen;
    for (int i = 0; i < replacements.length(); ++i) {
      HostPopulation(i, host) = replacements[i];
    }
  }
 return HostPopulation;
 }
  