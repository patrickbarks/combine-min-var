#include <Rcpp.h>
using namespace Rcpp;


int get_group_n(int x1, int x2, IntegerVector vec)
{
  int sum = 0;
  for (int i = x1; i <= x2; ++i) sum += vec[i];
  return sum;
}


CharacterVector get_group_members(int x1, int x2, CharacterVector vec_labels)
{
  int len = x2 - x1 + 1;
  CharacterVector labels(len);
  for (int i = 0, j = x1; j <= x2; ++i, ++j) labels[i] = vec_labels[j];
  return labels;
}


// [[Rcpp::export]]
List combine_min_var(IntegerVector vec, CharacterVector vec_labels, int n_groups, float min_var_init = 9999999.9)
{

  // check inputs
  if (n_groups > vec.size() || n_groups < 2) stop("n_groups must be > 1 and <= length(vec)");
  if (vec.size() != vec_labels.size()) stop("Arguments vec and vec_labels must be of same length");
  if (any(is_na(vec))) stop("Missing values are not permitted within vec");


  // set relevant variables
  int n_total = vec.size();                 // total number of initial groups
  int n_cuts = n_groups - 1;                // number of desired cut-points
  int count = 0;                            // permutation counter

  float variance;                           // variance at given combination of cut-points
  float min_var = min_var_init;             // minimum variance

  IntegerVector cuts(n_cuts);               // vector of cut-points
  IntegerVector maxima(n_cuts);             // vector of cut-point maxima
  IntegerVector group_ns(n_groups);         // vector of group sample sizes
  IntegerVector cuts_minVar(n_cuts);        // cut-points yielding minimum variance
  IntegerVector group_ns_minVar(n_groups);  // group sample sizes at cut-points yielding minimum variance

  bool exit = false;                        // exit status (break out of main loop if true)
  bool ties = false;                        // returns true if there are ties for min variance


  // set initial cut-points
  for (int i = 0; i < n_cuts; ++i) cuts[i] = i;


  // set cut point maxima
  for (int i = 0, j = n_total - n_cuts - 1; i < n_cuts; ++i, ++j) maxima[i] = j;


  // main loop
  while (!exit) { 


    // test for exit condition
    if (cuts[0] == maxima[0]) exit = true;


    // calculate sample size within each group
    group_ns[0] = get_group_n(0, cuts[0], vec);                                // first group
    group_ns[n_groups-1] = get_group_n(cuts[n_groups-2]+1, n_total-1, vec);    // last group
    if(n_groups > 2)                                                           // all other groups, if applicable
      for(int i = 1; i < n_groups-1; ++i)
        group_ns[i] = get_group_n(cuts[i-1]+1, cuts[i], vec);


    // calculate among-group variance in sample size
    variance = var(group_ns);


    // check for tied or minimum variance
    if (abs(variance - min_var) < 0.0000001)     // if tied for minimum variance
    {
      ties = true;
    } else if (variance < min_var)               // else, if minimum variance, set relevant vars
    {
      ties = false;                                                          // reset variable tracking ties
      min_var = variance;                                                    // set var_min to given variance
      for (int i = 0; i < n_groups; ++i) group_ns_minVar[i] = group_ns[i];   // set group sample sizes yielding min var
      for (int i = 0; i < n_cuts; ++i) cuts_minVar[i] = cuts[i];             // set cutoffs yielding min var
    }


    // increment cut-points, from highest (rightmost) to lowest (leftmost), until next combination reached
    int r = n_cuts-1;          // set cut-point index to max (highest/rightmost cut-point)

    while (r >= 0) {

      // increment given cut-point
      ++cuts[r];

      // reset higher cut-points (if any)
      if (r < n_cuts-1)
        for (int i = r+1; i < n_cuts; ++i)
          cuts[i] = cuts[i-1] + 1;

      // check whether given cut-point maxed
      if (cuts[r] > maxima[r])      // if so, move on to lower cut-point (decrement r) and continue
      {
        --r;
      } else {                      // else, break out of while loop (lower cut-points need not be changed)
        break;
      }
 
    }


    ++count;                        // increment cut-point combination counter
  }


  // ensure min_var < min_var_init
  if (min_var >= min_var_init) stop("Minimum variance !< min_var_init");


  // create group membership list given cut-points yielding minimum variance
  List group_members(n_groups);                                                                   // create list object
  group_members[0] = get_group_members(0, cuts_minVar[0], vec_labels);                            // first group
  group_members[n_groups-1] = get_group_members(cuts_minVar[n_cuts-1]+1, n_total-1, vec_labels);  // last group
  if(n_groups > 2)                                                                                // all other groups, if applicable
    for(int i = 1; i < n_groups-1; ++i)
      group_members[i] = get_group_members(cuts_minVar[i-1]+1, cuts_minVar[i], vec_labels);


  // list of output variables
  List out = List::create(count, ties, min_var, group_ns_minVar, group_members);
  out.names() = CharacterVector::create("combinations_searched", "ties", "minimum_variance", "group_sample_sizes", "group_membership");

  return out;

}


