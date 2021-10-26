#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector sortedIndex(NumericVector x){
  IntegerVector idx = seq_along(x) - 1;
  
  std::stable_sort(idx.begin(), idx.end(), [&](int i, int j){return x[i] < x[j];});
  
  return idx;
}

// [[Rcpp::export]]
IntegerVector compare_self(NumericVector x){
  int n_entry = x.size();
  IntegerVector match_self (n_entry);
  match_self[0] = 1;
  
  int idx = 1;
  
  for (int i = 1; i < (n_entry); i++) {
    if (x[i] != x[(i - 1)]) {
      match_self[idx] = 1;
    } else {
      match_self[idx] = 0;
    }
    idx++;
  }
  return match_self;
}

// [[Rcpp::export]]
IntegerVector compare_both(IntegerVector x, IntegerVector y){
  int n_entry = x.size();
  IntegerVector match_self (n_entry);
  match_self[0] = 1;
  
  int idx = 1;
  
  for (int i = 1; i < (n_entry); i++) {
    if ((x[i] != x[(i - 1)]) || (y[i] != y[(i - 1)])) {
      match_self[idx] = 1;
    } else {
      match_self[idx] = 0;
    }
    idx++;
  }
  match_self.push_back(1);
  return match_self;
}

// [[Rcpp::export]]
IntegerVector which_notzero(IntegerVector x){
  IntegerVector notzero (x.size());
  int idx = 0;
  
  for (int i = 0; i < x.size(); i++) {
    if (x[i] != 0) {
      notzero[idx] = i;
      idx++;
    }
  }
  IntegerVector keep_loc = seq(0, (idx - 1));
  notzero = notzero[keep_loc];
  return notzero;
}

// [[Rcpp::export]]
int kendall_discordant(IntegerVector x, IntegerVector y){
  double sup = 1 + max(y);
  
  IntegerVector arr(sup, 0);
  double i = 0;
  double k = 0;
  int n = x.size();
  int idx = 0;
  int dis = 0;
  
  while (i < n){
    while ((k < n) && (x[i] == x[k])) {
      dis = dis + i;
      idx = y[k];
      while (idx != 0) {
        dis = dis - arr[idx];
        idx = idx & (idx - 1);
      }
      k++;
    }
    while (i < k) {
      idx = y[i];
      while (idx < sup) {
        arr[idx] = arr[idx] + 1;
        idx = idx + (idx & (-1*idx));
      }
      i++;
    }
  }
  return dis;
}

// [[Rcpp::export]]
NumericVector count_rank_tie(IntegerVector ranks){
  
  LogicalVector dup_ranks(ranks.size());
  dup_ranks = duplicated(ranks);
  IntegerVector ranks2 = ranks[dup_ranks];
  IntegerVector number_tied;
  number_tied = table(ranks2) + 1;
  
  NumericVector counts(3);
  counts(0) = sum(number_tied * (number_tied - 1)) / 2;
  counts(1) = sum(number_tied * (number_tied - 1) * (number_tied - 2)) / 2;
  counts(2) = sum(number_tied * (number_tied - 1) * (2 * number_tied + 5));
  counts.names() = CharacterVector({"ntie", "t0", "t1"});
  
  return counts;
}


inline double signC(double x) {
  if (x > 0) {
    return 1.0;
  } else if (x == 0) {
    return 0.0;
  } else {
    return -1.0;
  }
}

//' Calculates ici-kendall-tau
//' 
//' @param x numeric vector
//' @param y numeric vector
//' @param perspective should we consider the "local" or "global" perspective?
//' @param alternative what is the alternative for the p-value test?
//' @param continuity logical: if true, a continuity correction is used
//' @param output used to control reporting of values for debugging
//' 
//' @details Calculates the information-content-informed Kendall-tau correlation measure.
//'   This correlation is based on concordant and discordant ranked pairs, like Kendall-tau,
//'   but also includes missing values (as NA). Missing values are assumed to be *primarily* due
//'   to lack of detection due to instrumental sensitivity, and therefore encode *some* information.
//'   
//'   For more details see the ICI-Kendall-tau vignette:
//'   
//'   \code{browseVignettes("ICIKendallTau")}
//' 
//' @examples 
//' x = sort(rnorm(100))
//' y = x + 1
//' y2 = y
//' y2[1:10] = NA
//' ici_kt(x, y)
//' ici_kt(x, y2, "global")
//' ici_kt(x, y2)
//' 
//' @importFrom Rcpp sourceCpp
//' @export
//' @useDynLib ICIKendallTau
//' @return kendall tau correlation, p-value, max-correlation
// [[Rcpp::export]]
NumericVector ici_kt(NumericVector x, NumericVector y, String perspective = "local", String alternative = "two.sided", bool continuity = false, String output = "simple") {
  
  if (x.length() != y.length()) {
    throw std::range_error("X and Y are not the same length!");
    exit(-1);
  }
  
  
  NumericVector z_b (1);
  NumericVector p_value (1);
  
  LogicalVector matching_na;
  //double n_matching_na;
  
  if (perspective == "local") {
    matching_na = is_na(x) & is_na(y);
    //n_matching_na = sum(matching_na);
    x = x[!matching_na];
    y = y[!matching_na];
  }
  
  NumericVector x2 = clone(x);
  NumericVector y2 = clone(y);
  
  int n_na_x = sum(is_na(x));
  int n_na_y = sum(is_na(y));
  
  if ((n_na_x == x.size()) || (n_na_y == y.size())) {
    NumericVector na_res (2);
    na_res(0) = NA_REAL;
    na_res(1) = NA_REAL;
    na_res.names() = CharacterVector({"tau", "pvalue"});
    return na_res;
  }
  
  x2 = x[!is_na(x)];
  y2 = y[!is_na(y)];
  
  double min_x = min(x2) - 0.1;
  double min_y = min(y2) - 0.1;
  
  x2 = clone(x);
  y2 = clone(y);
  x2[is_na(x)] = min_x;
  y2[is_na(y)] = min_y;
  
  
  int64_t n_entry = x2.size();
  //Rprintf("n_entry: %i\n", n_entry);
  
  if (n_entry < 2) {
    NumericVector na_res (2);
    na_res(0) = NA_REAL;
    na_res(1) = NA_REAL;
    na_res.names() = CharacterVector({"tau", "pvalue"});
    return na_res;
  }
  
  IntegerVector perm_y = sortedIndex(y2);
  x2 = x2[perm_y];
  y2 = y2[perm_y];
  IntegerVector y3 = compare_self(y2);
  IntegerVector y4 = cumsum(y3);
  //return y4;
  
  IntegerVector perm_x = sortedIndex(x2);
  x2 = x2[perm_x];
  y4 = y4[perm_x];
  IntegerVector x3 = compare_self(x2);
  IntegerVector x4 = cumsum(x3);
  
  //return x4;
  IntegerVector obs = compare_both(x4, y4);
  //return obs;
  int64_t sum_obs = sum(obs);
  IntegerVector cnt = diff(which_notzero(obs));
  int64_t dis = kendall_discordant(x4, y4);
  
  long double ntie = sum((cnt * (cnt - 1)) / 2);
  // three values should be read as:
  // xtie, x0, and x1, and then same for y
  NumericVector x_counts = count_rank_tie(x4);
  double xtie = x_counts[0];
  double x0 = x_counts[1];
  double x1 = x_counts[2];
  
  NumericVector y_counts = count_rank_tie(y4);
  double ytie = y_counts[0];
  double y0 = y_counts[1];
  double y1 = y_counts[2];
  
  int64_t tot = (n_entry * (n_entry - 1)) / 2;
  
  //Note that tot = con + dis + (xtie - ntie) + (ytie - ntie) + ntie
  //              = con + dis + xtie + ytie - ntie
  //Therefore con - dis = tot - xtie - ytie + ntie - 2 * dis
  //And concordant - discordant is the numerator for kendall-tau.
  //con + dis = tot - xtie - ytie + ntie
  //This is the maximum theoretical number of concordant values given
  //the data.
  
  NumericVector k_res(3);
  k_res.names() = CharacterVector({"tau", "pvalue", "tau_max"});
  if ((xtie == tot) || (ytie == tot)) {
    k_res(0) = NA_REAL;
    k_res(1) = NA_REAL;
    k_res(3) = NA_REAL;
    return k_res;
  }
  
  long double con_minus_dis = tot - xtie - ytie + ntie - 2 * dis;
  long double tau = con_minus_dis / sqrt((tot - xtie) * (tot - ytie));
  long double con_plus_dis = tot - xtie - ytie + ntie;
  long double tau_max = con_plus_dis / sqrt((tot - xtie) * (tot - ytie));
  if (tau > 1) {
    tau = 1;
  } else if (tau < -1) {
    tau = -1;
  }
  
  int64_t m = n_entry * (n_entry - 1);
  //Rprintf("m: %f\n", m);
  long double var = ((m * (2 * n_entry + 5) - x1 - y1) / 18 +
                (2 * xtie * ytie) / m + x0 * y0 / (9 * m * (n_entry - 2)));
  //Rprintf("var: %f\n", var);
  long double s_adjusted = tau * sqrt(((m / 2) - xtie) * ((m / 2) - ytie));
  if (continuity) {
    long double adj_s2 = signC(s_adjusted) * (std::abs(s_adjusted) - 1);
    s_adjusted = adj_s2;
  }
  long double z_b_0 = s_adjusted / sqrt(var);
  z_b[0] = s_adjusted / sqrt(var);
  
  if (alternative == "less") {
    k_res[1] = pnorm(z_b, 0.0, 1.0)[0];
  } else if (alternative == "greater") {
    k_res[1] = pnorm(z_b, 0.0, 1.0, false, false)[0];
  } else if (alternative == "two.sided") {
    NumericVector p_res (2);
    p_res[0] = pnorm(z_b, 0.0, 1.0)[0];
    p_res[1] = pnorm(z_b, 0.0, 1.0, false)[0];
    k_res[1] = 2 * min(p_res);
  }
  k_res[0] = tau;
  k_res[2] = tau_max;
  
  //Rprintf("n_entry: %f\n", n_entry);
  
  if (output != "simple") {
    std::string report_st = "min_x: " + std::to_string(min_x) + "\n" +
      "min_y: " + std::to_string(min_y) + "\n" +
      "n_entry: " + std::to_string(n_entry) + "\n" +
      "tot: " + std::to_string(tot) + "\n" +
      "sum_obs: " + std::to_string(sum_obs) + "\n" +
      "dis: " + std::to_string(dis) + "\n" +
      "con_minus_dis (k_numerator): " + std::to_string(con_minus_dis) + "\n" +
      "n_tie: " + std::to_string(ntie) + "\n" +
      "m: " + std::to_string(m) + "\n" +
      "x_tie: " + std::to_string(xtie) + "\n" +
      "y_tie: " + std::to_string(ytie) + "\n" +
      "s_adjusted: " + std::to_string(s_adjusted) + "\n" +
      "var: " + std::to_string(var) + "\n" +
      "z_b: " + std::to_string(z_b_0) + "\n" +
      "tau: " + std::to_string(tau) + "\n" +
      "tau_max:" + std::to_string(tau_max) + "\n";
      "pvalue: " + std::to_string(k_res[1]) + "\n";
    Rcout << report_st;
  }
  
  return k_res;
}


// [[Rcpp::export]]
NumericVector ici_kt_pairs(NumericVector x, NumericVector y, String perspective = "local", String alternative = "two.sided", String output = "simple") {
  
  if (x.length() != y.length()) {
    throw std::range_error("X and Y are not the same length!");
    exit(-1);
  }
  
  double sum_concordant = 0;
  double sum_discordant = 0;
  double sum_x_ties = 0;
  double sum_y_ties = 0;
  //double sum_tied_x = 0;
  //double sum_tied_y = 0;
  double sum_tied_x_na = 0;
  //double sum_tied_y_na = 0;
  double sum_all_na = 0;
  double k_numerator;
  double k_denominator;
  double k_tau;
  //bool reject_concordant;
  //bool reject_discordant;
  
  // for generating the p-value
  //double xties = 0;
  //double yties = 0;
  double t_0 = 0;
  double s_adjusted = 0;
  double x_tied_sum_t1;
  double y_tied_sum_t2;
  double v_0_sum = 0;
  double v_t_sum = 0;
  double v_u_sum = 0;
  double v_t1_sum = 0;
  double v_t2_sum = 0;
  double s_adjusted_variance = 0;
  NumericVector z_b (1);
  NumericVector p_value (1);
  
  LogicalVector matching_na;
  //double n_matching_na;
  
  if (perspective == "local") {
    matching_na = is_na(x) & is_na(y);
    //n_matching_na = sum(matching_na);
    x = x[!matching_na];
    y = y[!matching_na];
  }
  
  NumericVector x2 = clone(x);
  NumericVector y2 = clone(y);
  
  int n_na_x = sum(is_na(x));
  int n_na_y = sum(is_na(y));
  
  if ((n_na_x == x.size()) || (n_na_y == y.size())) {
    return 0.0;
  }
  
  x2 = x[!is_na(x)];
  y2 = y[!is_na(y)];
  
  double min_value = min(NumericVector::create(min(x2), min(y2)));
  double na_value = min_value - 0.1;
  
  x2 = clone(x);
  y2 = clone(y);
  x2[is_na(x)] = na_value;
  y2[is_na(y)] = na_value;
  
  
  double n_entry = x2.size();
  
  if (n_entry < 2) {
    return 0.0;
  }
  
  
  for (int i = 0; i < (n_entry - 1); i++) {
    for (int j = (i+1); j < n_entry; j++) {
      sum_concordant+= (signC(x2[i] - x2[j]) * signC(y2[i] - y2[j])) > 0;
      sum_discordant+= (signC(x2[i] - x2[j]) * signC(y2[i] - y2[j])) < 0;
    }
  }

  // if (perspective == "global") {
  //   sum_x_ties = sum_tied_x + sum_tied_x_na + half_sum_na_ties;
  //   sum_y_ties = sum_tied_y + sum_tied_y_na + half_sum_na_ties;
  // } else {
  //   sum_x_ties = sum_tied_x;
  //   sum_y_ties = sum_tied_y;
  // }
  // 
  k_numerator = sum_concordant - sum_discordant;
  
  LogicalVector dup_x(x2.size());
  dup_x = duplicated(x2);
  NumericVector x3 = x2[dup_x];
  NumericVector x_tied_values_t1;
  x_tied_values_t1 = table(x3) + 1;
  
  LogicalVector dup_y(y2.size());
  dup_y = duplicated(y2);
  NumericVector y3 = y2[dup_y];
  NumericVector y_tied_values_t2;
  y_tied_values_t2 = table(y3) + 1;
  
  t_0 = n_entry * (n_entry - 1) / 2;
  x_tied_sum_t1 = sum(x_tied_values_t1 * (x_tied_values_t1 - 1)) / 2;
  //Rprintf("x_tied_sum_t1: %f\n", x_tied_sum_t1);
  y_tied_sum_t2 = sum(y_tied_values_t2 * (y_tied_values_t2 - 1)) / 2;
  
  k_denominator = sqrt((t_0 - x_tied_sum_t1) * (t_0 - y_tied_sum_t2));
  
  if (k_denominator == 0) {
    k_tau = 0;
  } else {
    k_tau = k_numerator / k_denominator;
  }
  
  
  // p-value calculation
  
  
  s_adjusted = k_tau * sqrt((t_0 - x_tied_sum_t1) * (t_0 - y_tied_sum_t2));
  v_0_sum = n_entry * (n_entry - 1) * (2 * n_entry + 5);
  v_t_sum = sum(x_tied_values_t1 * (x_tied_values_t1 - 1) * (2 * x_tied_values_t1 + 5));
  v_u_sum = sum(y_tied_values_t2 * (y_tied_values_t2 - 1) * (2 * y_tied_values_t2 + 5));
  v_t1_sum = sum(x_tied_values_t1 * (x_tied_values_t1 - 1)) * sum(y_tied_values_t2 * (y_tied_values_t2 - 1));
  v_t2_sum = sum(x_tied_values_t1 * (x_tied_values_t1 - 1) * (x_tied_values_t1 - 2)) * sum(y_tied_values_t2 * (y_tied_values_t2 - 1) * (y_tied_values_t2 - 2));

  s_adjusted_variance = (v_0_sum - v_t_sum - v_u_sum) / 18 +
    v_t1_sum / (2 * n_entry * (n_entry - 1)) +
    v_t2_sum / (9 * n_entry * (n_entry - 1) * (n_entry - 2));

  double s_adjusted2 = signC(s_adjusted) * (std::abs(s_adjusted) - 1);
  z_b[0] = s_adjusted2 / sqrt(s_adjusted_variance);
  if (alternative == "less") {
    p_value[0] = pnorm(z_b, 0.0, 1.0)[0];
  } else if (alternative == "greater") {
    p_value[0] = pnorm(z_b, 0.0, 1.0, false, false)[0];
  } else if (alternative == "two.sided") {
    NumericVector p_res (2);
    p_res[0] = pnorm(z_b, 0.0, 1.0)[0];
    p_res[1] = pnorm(z_b, 0.0, 1.0, false)[0];
    p_value[0] = 2 * min(p_res);
  }

  
  // debugging
  double con_plus_dis = sum_concordant + sum_discordant;
  if (output != "simple") {
    Rprintf("min_value: %f \n", min_value);
    Rprintf("na_value: %f \n", na_value);
    
    Rprintf("n_entry: %f \n", n_entry);
    //Rprintf("x_ties: %f \n", sum_tied_x);
    //Rprintf("x_na_ties: %f \n", sum_tied_x_na);
    Rprintf("sum_x_ties: %f \n", sum_x_ties);
    Rprintf("sum_y_ties: %f \n", sum_y_ties);
    Rprintf("sum_na_ties: %f \n", sum_all_na);
    //Rprintf("half_sum_na_ties: %f \n", half_sum_na_ties);
    Rprintf("sum_concordant: %f \n", sum_concordant);
    Rprintf("sum_discordant: %f \n", sum_discordant);
    Rprintf("k_numerator: %f \n", k_numerator);
    Rprintf("con_plus_dis: %f\n", con_plus_dis);
    Rprintf("k_denominator: %f \n", k_denominator);
    Rprintf("t_0: %f \n", t_0);
    Rprintf("x_tied_sum_t1: %f \n", x_tied_sum_t1);
    Rprintf("y_tied_sum_t2: %f \n", y_tied_sum_t2);
    Rprintf("s_adjusted: %f \n", s_adjusted2);
    Rprintf("s_adjusted_variance: %f \n", s_adjusted_variance);
    Rprintf("k_tau: %f \n", k_tau);
    Rprintf("pvalue: %f \n", p_value[0]);
  }
  
  NumericVector out_values = {k_tau, p_value[0]};
  out_values.names() = CharacterVector({"tau", "pvalue"});
  
  return out_values;
}
