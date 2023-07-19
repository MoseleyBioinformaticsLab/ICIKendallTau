#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame pairwiseComparisons(CharacterVector entries, int n_cores, bool include_self = false){
  int64_t n_entries = entries.length();
  int64_t n_todo = (n_entries * (n_entries - 1)) / 2;
  Rprintf("n_todo: %i\n", n_todo);
  
  if (include_self) {
    n_todo = n_todo + n_entries;
  }
  //Rprintf("n_todo: %i\n", n_todo);
  
  long double n_each = (long double)n_todo / n_cores;
  Rprintf("n_each: %f\n", n_each);
  long double n_each2 = ceil(n_each);
  Rprintf("n_each 2: %f\n", n_each2);
  
  CharacterMatrix compare_matrix(n_todo, 3);
  
  int entry_loc = 0;
  int list_loc = 0;
  int check_each = 0;
  
  for (int i_entry = 0; i_entry < n_entries; i_entry++) {
    for (int j_entry = (include_self) ? i_entry : i_entry + 1; j_entry < n_entries; j_entry++) {
      //Rprintf("i: %i\nj: %i\nentry: %i\nlist: %i\ncheck: %i\n", i_entry, j_entry, entry_loc, list_loc, check_each);
      
      compare_matrix(entry_loc, 0) = entries[i_entry];
      compare_matrix(entry_loc, 1) = entries[j_entry];
      compare_matrix(entry_loc, 2) = (char)list_loc;
      
      entry_loc++;
      check_each++;
      
      if (check_each == n_each2) {
        list_loc++;
        check_each = 0;
      }
    }
  }
  
  return DataFrame::create(Named("C1") = compare_matrix(_, 0),
                           Named("C2") = compare_matrix(_, 1),
                           Named("Group") = compare_matrix(_, 2));
}
