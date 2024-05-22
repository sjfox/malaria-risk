#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
List subs_parms(List sub_parms, List ref_parms){
  StringVector ref_names = ref_parms.names();
  StringVector sub_names = sub_parms.names();
  for(int j = 0; j < ref_parms.size(); j++){
    for(int i =0; i < sub_parms.size(); i++){
      if(ref_names[j] == sub_names[i]){
        ref_parms[j] = sub_parms[i];
      }
    }
  }
  return(ref_parms);
}

// [[Rcpp::export]]
NumericVector find_rnot_ods(NumericVector rnot){
  // ## Takes in a vector of rnot values and returns vector of same length of dispersions
  // Dispersion parameter is drawn to replicate 80/20 rule when R0 = 0.5, as malaria reproduction numbers
  // all seem to be well below 1 
  NumericVector ods(rnot.size(), 0.169);
  return(ods);
}

// [[Rcpp::export]]
double offspring_size_llike(double sec_trans, double ods, double rnots){
  // Takes in likelihood parms, and calculates the likelihood of an outbreak of size sec_trans
  //exp(lgamma(k*j+j-1)-lgamma(k*j)-lgamma(j+1)+(j-1)*log(R0/k)-(k*j+j-1)*log(1+R0/k))
  double log_like = R::lgammafn( ods * sec_trans + sec_trans - 1) -
                      R::lgammafn( ods * sec_trans) -
                      R::lgammafn(sec_trans + 1) +
                      (sec_trans - 1) * log( (rnots / ods)) -
                      (ods * sec_trans+sec_trans - 1) * log(1 + rnots/ods);
  // Rcout << R::lgammafn(sec_trans + 1) << ", " << sec_trans << std::endl;
  return(log_like);

}

// [[Rcpp::export]]
double offspring_size_obs_llike(double sec_trans, double ods, double rnots, double reporting_rate){
  // Takes in likelihood parms, and calculates the likelihood of obsreving an outbreak of size sec_trans
  // issue is that sec_trans is defined as being 0 if outbreak is size 1, so now adding 1 to clear up
  // see `sim_test_likelihoood_fxns.R` file for likelihood comparison with simulation results
  sec_trans = sec_trans + 1;
  IntegerVector true_outbreak_sizes = seq(sec_trans, 50);

  double total_prob = 0.0;
  double prob_outbreak_one_size = 0.0;
  for(int i = 0; i < true_outbreak_sizes.size(); i++){
    prob_outbreak_one_size = exp(offspring_size_llike(true_outbreak_sizes[i], ods, rnots) + R::dbinom((sec_trans - 1), (true_outbreak_sizes[i] - 1), reporting_rate, true));
    total_prob += prob_outbreak_one_size;
    // Rcout << "Probability is" << R::dbinom((sec_trans - 1), (true_outbreak_sizes[i] - 1), reporting_rate, true) << std::endl <<
      // "Variables are: sec_trans=" << true_outbreak_sizes[i] << " ods=" <<ods<<" rnot=" << rnots << " rr=" << reporting_rate << std::endl;
  }
  // Rcout << "Probability is" << total_prob << std::endl <<
    // "Variables are: sec_trans=" << sec_trans << " ods=" <<ods<<" rnot=" << rnots << " rr=" << reporting_rate << std::endl;
  return(log(total_prob));
}


// Old dispersion parameter finding algorithm. Not needed because no longer have the table for searching
// // [[Rcpp::export]]
// NumericVector find_rnot_ods(NumericVector rnot, DataFrame disp_df){
//   // ## Takes in a vector of rnot values and returns
//   // ## The dispersions for nbinom distribution
//   // ## disp_dt must be sorted data table, with one column as rnot and other as ods
//   // ## Can be obtained by running the calc_dispersion_table.R script
//   IntegerVector indices(rnot.size());
//   NumericVector rnots = disp_df["rnots"];
//   for(int i =0; i < rnot.size(); i++){
//     indices[i] = std::distance(rnots.begin(), std::lower_bound(rnots.begin(), rnots.end(), rnot[i]));
//     if(indices[i] != 0){
//       if( std::abs(rnots[indices[i] - 1] - rnots[i]) < std::abs(rnots[indices[i]] - rnots[i])){
//         indices[i] --;
//       }
//     }
//   }
//   NumericVector ods = disp_df["ods"];
//   return(ods[indices]);
// }

// [[Rcpp::export]]
NumericVector intro_loglike(List parms) {
// ## Calculates the likelihood based on an R0 and number of introductions
  String dist = Rcpp::as<String>(parms["distribution"]);
  NumericVector rnots = Rcpp::as<NumericVector>(parms["rnot"]);
  NumericVector ods = Rcpp::as<NumericVector>(parms["overdispersion"]);
  NumericVector num_intros = Rcpp::as<NumericVector>(parms["num_intros"]);
  NumericVector sec_trans = Rcpp::as<NumericVector>(parms["secondary_trans"]);
  double reporting_rate = as<double>(parms["reporting_rate"]);
  NumericVector log_likes(rnots.size());
  if(dist =="pois"){
    for(int i =0; i < log_likes.size();i++){
      log_likes[i] = R::dpois(sec_trans[i], rnots[i], true) * num_intros[i];
    }
  } else if(dist =="nbinom"){
    for(int i =0; i < log_likes.size(); i++){

      // log_likes[i] = R::dnbinom_mu(sec_trans[i], ods[i], rnots[i], true) * num_intros[i];
      log_likes[i] = offspring_size_obs_llike(sec_trans[i], ods[i], rnots[i], reporting_rate) * num_intros[i];

    }
  } else{
    Rcpp::stop("Inadmissible distribution value");
  }
  return(log_likes);
}

// [[Rcpp::export]]
double scaling_loglike_cpp(double alpha, List params){
// ## Returns the log likelihood for set of parameters
  List parms = Rcpp::clone(params);
  NumericVector rnots = as<NumericVector>(parms["rnot"]);
  double reporting_rate = as<double>(parms["reporting_rate"]);

  if(any(!Rcpp::is_na(rnots))){
    // This used in mcmc sampling

    rnots = rnots * alpha; //* reporting_rate;
    // Rcout << "The value of r0 : " << rnots << "\n";
    NumericVector ods = find_rnot_ods(rnots);
    parms = subs_parms(Rcpp::List::create(Rcpp::Named("rnot") = rnots,
                                          Rcpp::Named("overdispersion") = ods), parms);

    return(Rcpp::sum(intro_loglike(parms)));
  } else{
    // This used when sampling from distributions - no longer necessary
    NumericMatrix rnot_dist = internal::convert_using_rfunction(as<DataFrame>(parms["rnot_dist"]), "as.matrix");
    rnot_dist = rnot_dist * alpha * reporting_rate;
    NumericVector rnots = internal::convert_using_rfunction(rnot_dist, "as.numeric");
    NumericVector log_likes(rnot_dist.ncol());
    NumericVector ods = find_rnot_ods(rnots);
    // ods.attr("dim") = Dimension(rnot_dist.nrow(), rnot_dist.ncol());

    for(int i =0; i < rnot_dist.ncol(); i++){
      IntegerVector indices = Rcpp::seq(0, (rnot_dist.nrow()-1));
      for(int j =0; j < indices.size(); j++){
        indices[j] = indices[j] + i * (rnot_dist.nrow());
      }
      parms = subs_parms(Rcpp::List::create(Rcpp::Named("rnot") = rnots[indices],
                                            Rcpp::Named("overdispersion")=ods[indices]), parms);
      log_likes[i] = -Rcpp::sum(intro_loglike(parms));
    }
    return(Rcpp::mean(log_likes));
  }
}

NumericVector range(double min, double max, int N) {
  NumericVector range(N);
  double delta = (max-min) / double(N-1);
  for(int i=0; i<N; i++) {
    range[i] = min + i*delta;
  }
  return range;
}


// [[Rcpp::export]]
DataFrame get_alpha_likes_cpp(List parms){
  // # Returns likelihood values for a variety of alphas, so that distributions can be calculated post-hoc

  CharacterVector date = as<CharacterVector>(parms["date"]);
  std::string date_str = as<std::string>(date[0]);
  int n = 1000;
  NumericVector alphas = range(0.0, 1.0, n);
  NumericVector nllikes(n);
  for(int i =0; i <n; i++){
    // std::cout << alphas[i] << std::endl;
    nllikes[i] = scaling_loglike_cpp(alphas[i], parms);
  }

  nllikes = Rcpp::exp(-nllikes);

  return Rcpp::DataFrame::create(Rcpp::Named("alpha") = alphas, Rcpp::Named(date_str) = nllikes);
}

// double get_secondary_above_20(double rnot){
//   // Takes in an rnot value and returns the probability of seeing
//   // > 20 secondary cases from that rnot
//   NumericVector p1 = NumericVector::create(0.425806451612903, 0.8458765530605259);
//   NumericVector p2 = NumericVector::create(4.341935483870967, 3.297197366921235);
//   double yint = - ((p2[1] - p1[1]) / (p2[0] - p1[0])) * p1[0] + p1[1];
//   if(rnot < yint){
//     // warning("R0 is low and returning zero") // Happens very often, so not worth warning
//     return(0.0);
//   }
//   return ((rnot - yint) / ((p2[1] - p1[1]) / (p2[0] - p1[0]))  / 100);
// }
//
// // [[Rcpp::export]]
// NumericVector compare_ps(NumericVector x, double prob_above, double rnot){
//   NumericVector ps(x.size());
//   for(int i =0; i < x.size(); i ++){
//     ps[i] = R::pnbinom_mu(20, x[i], rnot, 0, 0) - prob_above;
//   }
//   return(ps);
// }



//
// NumericVector get_alpha_likes_cpp(List parms){
//   // Returns likelihood values for a variety of alphas, so that distributions can be calculated post-hoc
//   NumericVector alphas = range(0.0, 1.0, 5000);
//   NumericVector nllikes(alphas.size());
//   Function loglike("scaling_loglike");
//
//   for(int i =0; i < alphas.size(); i++){
//     nllikes[i] = as<NumericVector>(loglike(Rcpp::wrap(parms)));
//   }
//
//   return(Rcpp::exp(-nllikes));
// }

// // [[Rcpp::export]]
// List find_overdispersion(double rnot){
// // Find the overdispersion parameter for a given R0
//   double prob_above = get_secondary_above_20(rnot);
//
//   Function uniroot = Environment("package:stats")["uniroot"];
//
//   // print(rnot)
//   // if(prob_above == 0){
//   //   if(rnot==0) {
//   //     return(1e-16);
//   //   }
//   //
//   //   NumericVector seq_lower = NumericVector::create(range(1e-16, 0.5, 1000));
//   //   double low = -1;
//   //   double prob_above = 1e-5;
//   //   NumericVector ps = compare_ps(seq_lower, prob_above, rnot);
//   //
//   //   while(Rcpp::max(ps) < 0){
//   //     prob_above = prob_above/10;
//   //     ps = compare_ps(seq_lower, prob_above, rnot);
//   //   }
//   //   low = seq_lower[Rcpp::which_max(ps)];
//   //
//   //
//   //   List overdisp = f(Rcpp::Named("f", compare_ps),
//   //                     Rcpp::Named("interval", NumericVector::create(low, 1)),
//   //                     Rcpp::Named("rnot", rnot),
//   //                     Rcpp::Named("prob_above", prob_above));
//   //   return(overdisp["root"]);
//   // } else {
//     NumericVector seq_lower = range(0, 0.5, 1000);
//     NumericVector ps = compare_ps(seq_lower, prob_above, rnot);
//     int ind = 0;
//     for(int i =0; i<ps.size();i++){
//       if(ps[i] >= 0){
//         break;
//       }
//       ind++;
//     }
//     double low = seq_lower[ind];
//     List overdisp = uniroot(Rcpp::Named("f") =compare_ps,
//                       Rcpp::Named("interval")= NumericVector::create(low, 10),
//                       Rcpp::Named("rnot")= rnot,
//                       Rcpp::Named("prob_above") = prob_above);
//     return(overdisp);
//   // }
// }



// // [[Rcpp::export]]
// List uniroot_c(){
//   Function uniroot = Environment("package:stats")["uniroot"];
//   List overdisp = uniroot(Rcpp::Named("f") = get_ys, Rcpp::Named("interval") = NumericVector::create(0,50));
//   return(overdisp);
// }
