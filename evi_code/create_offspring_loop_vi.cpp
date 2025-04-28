#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double fitness_func_bacgen_cpp(double selection_parameter, double optima1, double optima2, double trait,double weighting) {
  // Calculate the mean of optima1 and optima2
  double meanval = (weighting * optima1) + ((1-weighting)*optima2);
  
  // Calculate fitness function
  double q = exp(pow((trait - meanval), 2) / -selection_parameter);
  
  return q;
}

// Function to simulate multinomial distribution in C++
// This is a simple version, replace with your actual implementation if needed

// [[Rcpp::export]]
NumericVector rmultinom_cpp(int n, NumericVector prob) {
  SEXP res = Rf_allocVector(INTSXP, prob.size());
  prob = prob / sum(prob); // Adjust probabilities to sum to 1
  rmultinom(n, prob.begin(), prob.size(), INTEGER(res));
  return NumericVector(res);
}



// [[Rcpp::export]]
List offspring_loopfunc(NumericMatrix host_pop, NumericVector ENV_sampling_probability, NumericVector host_microbe_optima,
                                 double X, double selection_parameter_microbes, double env_condition, NumericVector microbe_trait_list,
                                 int n_micro, StringVector microbe_names, double weighting){ //, NumericVector fitness_microbes, NumericVector weighted_samplingprob) {
//  Rcpp::Rcout << "weighting: " << X << std::endl;
  int N_Species = host_pop.nrow();
  int n_hosts = host_pop.ncol();
  NumericVector fitness_microbes(n_hosts);
  NumericVector weighted_samplingprob(n_hosts);
  //Rcpp::Rcout << "Starting loop\n"; // Debug print
  for (int i = 0; i < n_hosts; i++) {

    
    NumericVector parents = host_pop(_, i); // Selects a parent to reproduce
    // Rcpp::Rcout << "parents: " << parents << std::endl;
    NumericMatrix microbe_probs(N_Species, 10);

//  colnames(microbe_probs) = CharacterVector::create("parprob", "envprob", "ParProvided", "EnvProvided", "sampling_probability",
 //            "MicrobialFitness", "MicrobialFitnessAbs", "EnvDistance", "Env_Cond", "sampled_counts");
    rownames(microbe_probs) = microbe_names;
    
    microbe_probs(_, 0) = parents / n_micro; // parprob

    
    microbe_probs(_, 1) = ENV_sampling_probability; // envprob
    // Rcpp::Rcout << "Host Pop Element: ";
    // 
    // for (int i = 0; i < microbe_probs.nrow(); i++) {
    //   Rcpp::Rcout << microbe_probs(i, 1) << " ";
    // }
    // Rcpp::Rcout << std::endl;
    microbe_probs(_, 2) = X * microbe_probs(_, 0); // ParProvided
    microbe_probs(_, 3) = (1 - X) * microbe_probs(_, 1); // EnvProvided
    microbe_probs(_, 4) = microbe_probs(_, 3)  + microbe_probs(_, 2); //rowSums(microbe_probs(_, Range(2, 3)), true); // sampling_probability
//    microbe_probs(_, Range(0, 9)).attr("dimnames") = List::create(microbe_names, colnames(microbe_probs));
    for (int z = 0; z < N_Species; ++z) {
      microbe_probs(z, 5) = fitness_func_bacgen_cpp(selection_parameter_microbes, host_microbe_optima[i], env_condition, microbe_trait_list[z], weighting);
    }
    
    // microbe_probs(_, 5) = microbe_probs(_, 6); // MicrobialFitness (assuming no normalization)
    
    // Normalize the sampling probability by microbial fitness
    NumericVector sampling_probability = microbe_probs(_, 5) * microbe_probs(_, 4);

    
    sampling_probability = sampling_probability / sum(sampling_probability);
    
    microbe_probs(_, 6) = sampling_probability;
    NumericVector sampled_counts = rmultinom_cpp(n_micro, sampling_probability);
    microbe_probs(_, 9) = sampled_counts; // sampled_counts
    host_pop(_, i) = sampled_counts; // Update host_pop with sampled_counts
    
    // Calculate fitness_microbes and weighted_samplingprob
    fitness_microbes[i] = sum(microbe_probs(_, 5) * microbe_probs(_, 9)) / sum(microbe_probs(_, 9));
    weighted_samplingprob[i] = sum(microbe_probs(_, 6) * microbe_probs(_, 9)) / sum(microbe_probs(_, 9));
  }
  // Create and return the list
  return List::create(Named("host_pop") = host_pop,
                      Named("fitness_microbes") = fitness_microbes,
                      Named("weighted_samplingprob") = weighted_samplingprob);
}



// [[Rcpp::export]]
List offspring_loopfunc_vi(NumericMatrix host_pop, NumericVector ENV_sampling_probability, NumericVector host_microbe_optima,
                                 double X, double selection_parameter_microbes, double env_condition, NumericVector microbe_trait_list,
                                 int n_micro, StringVector microbe_names, double weighting){ //, NumericVector fitness_microbes, NumericVector weighted_samplingprob) {
//  Rcpp::Rcout << "weighting: " << X << std::endl;
  int N_Species = host_pop.nrow();
  int n_hosts = host_pop.ncol();
  NumericVector fitness_microbes(n_hosts);
  NumericVector weighted_samplingprob(n_hosts);
  //Rcpp::Rcout << "Starting loop\n"; // Debug print
  for (int i = 0; i < n_hosts; i++) {

    
    NumericVector parents = host_pop(_, i); // Selects a parent to reproduce
    // Rcpp::Rcout << "parents: " << parents << std::endl;
    NumericMatrix microbe_probs(N_Species, 10);

//  colnames(microbe_probs) = CharacterVector::create("parprob", "envprob", "ParProvided", "EnvProvided", "sampling_probability",
 //            "MicrobialFitness", "MicrobialFitnessAbs", "EnvDistance", "Env_Cond", "sampled_counts");
    rownames(microbe_probs) = microbe_names;
    
    microbe_probs(_, 0) = parents / n_micro; // parprob
    
     if (i == 0) {
      // Shift values in microbe_probs(_, 0) down by 200 positions
      NumericVector shifted_microbe_probs(N_Species);
      for (int z = 0; z < 200; ++z) {
        shifted_microbe_probs[z] = 0;
      }
      for (int z = 200; z < N_Species; ++z) {
        shifted_microbe_probs[z] = microbe_probs(z - 200, 0);
      }
      microbe_probs(_, 0) = shifted_microbe_probs;
    }
    
    microbe_probs(_, 1) = ENV_sampling_probability; // envprob
    // Rcpp::Rcout << "Host Pop Element: ";
    // 
    // for (int i = 0; i < microbe_probs.nrow(); i++) {
    //   Rcpp::Rcout << microbe_probs(i, 1) << " ";
    // }
    // Rcpp::Rcout << std::endl;
    microbe_probs(_, 2) = X * microbe_probs(_, 0); // ParProvided
    microbe_probs(_, 3) = (1 - X) * microbe_probs(_, 1); // EnvProvided
    microbe_probs(_, 4) = microbe_probs(_, 3)  + microbe_probs(_, 2); //rowSums(microbe_probs(_, Range(2, 3)), true); // sampling_probability
//    microbe_probs(_, Range(0, 9)).attr("dimnames") = List::create(microbe_names, colnames(microbe_probs));
    for (int z = 0; z < N_Species; ++z) {
      microbe_probs(z, 5) = fitness_func_bacgen_cpp(selection_parameter_microbes, host_microbe_optima[i], env_condition, microbe_trait_list[z], weighting);
    }
    
    // microbe_probs(_, 5) = microbe_probs(_, 6); // MicrobialFitness (assuming no normalization)
    
    // Normalize the sampling probability by microbial fitness
    NumericVector sampling_probability = microbe_probs(_, 5) * microbe_probs(_, 4);

    
    sampling_probability = sampling_probability / sum(sampling_probability);
    
    microbe_probs(_, 6) = sampling_probability;
    NumericVector sampled_counts = rmultinom_cpp(n_micro, sampling_probability);
    microbe_probs(_, 9) = sampled_counts; // sampled_counts
    host_pop(_, i) = sampled_counts; // Update host_pop with sampled_counts
    
    // Calculate fitness_microbes and weighted_samplingprob
    fitness_microbes[i] = sum(microbe_probs(_, 5) * microbe_probs(_, 9)) / sum(microbe_probs(_, 9));
    weighted_samplingprob[i] = sum(microbe_probs(_, 6) * microbe_probs(_, 9)) / sum(microbe_probs(_, 9));
  }
  // Create and return the list
  return List::create(Named("host_pop") = host_pop,
                      Named("fitness_microbes") = fitness_microbes,
                      Named("weighted_samplingprob") = weighted_samplingprob);
}

