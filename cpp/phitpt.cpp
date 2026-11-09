#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Saddlepoint approximation for the latent multinomial model for extended batch-marking data
  // Data objects
  DATA_VECTOR(obs_counts);
  DATA_SPARSE_MATRIX(linkmat);
  DATA_IMATRIX(true_hists);
  DATA_IMATRIX(true_prim_hists);
  DATA_IMATRIX(pp_ind);
  DATA_IMATRIX(eventone_ind);
  DATA_SCALAR(Nmin);
  DATA_INTEGER(nocca);
  DATA_INTEGER(nper);
  DATA_INTEGER(nhist);
  DATA_VECTOR(phi_lambda);
  DATA_VECTOR(p_lambda);
  DATA_VECTOR(beta0_lambda);
  DATA_INTEGER(penalty);
  
  // Parameters (after transformation for the purpose of unconstrained optimization)
  PARAMETER(Nt); 
  PARAMETER_VECTOR(pt);
  PARAMETER_VECTOR(phit); 
  PARAMETER_VECTOR(beta0t);
  
  // Random effects
  PARAMETER_VECTOR(svec);
  
  // Integers used in some for loops
  int i, j, k;
  
  // Original parameters
  Type N = Nmin + exp(Nt);
  vector<Type> p = exp(pt) / (1 + exp(pt));
  vector<Type> phi = exp(phit) / (1 + exp(phit));
  
  // Compute entry probabilities
  Type beta0 = 0.0;
  Type cumsum = 0.0;
  vector<Type> beta(nper);
  vector<Type> betat(nper);
  
  for(j=0;j < nper; j++){
    // Compute probability of entry conditional on not having entered yet
    if(j == nper - 1){
      beta0 = 1;
    }
    else{
      beta0 = exp(beta0t[j]) / (1 + exp(beta0t[j]));
    }
    // Compute marginal probability of entry
    beta[j] = beta0 * (1 - cumsum);
    betat[j] = log(beta[j]) - log(1 - beta[j]);
    // Update cumulative sum of entry probabilities
    cumsum += beta[j];
  }
  
  // ---------------------------------------------------------------------------------------------------
  // First of all, construct a probability matrix, pmat, which consists of nocca rows and 3 columns.
  // The element at the i-th row and j-th column of pmat is the probability of capture code (j-1) at time i.
  // For example, the first row of the matrix should be: (1-p[0], p[0], p[0]).
  // Probabilities of events 0 - 2 on each occasion in terms of model parameters
  // Pr(0 on occasion t) = 1 - p[t]
  // Pr(1 on occasion t) = p[t]
  // Pr(2 on occasion t) = p[t]
  matrix<Type> pmat(nocca, 4);
  pmat.col(0) = 1 - p;
  pmat.col(1) = p;
  pmat.col(2) = p;
  
  // REPORT(pmat);
  
  // How to use pmat for some convenience:
  // For 4 capture occasions, the cell probability (without survival) for history 1022 can be obtained:
  // Pr("1022") = pmat(0, 1) * pmat(1, 0) * pmat(2, 2) * pmat(3, 2).
  // A vector to store the probabilities:
  vector<Type> true_hists_probs(nhist);
  true_hists_probs.fill(1.0);    // Initialization
  
  // Define some variables for later use
  vector<int> true_hist(nocca);
  vector<Type> pvec(nper);
  int ppstart, ppend;
  Type tmp, tmpprob;
  
  // First calculate the probability of the null history, which is different from others
  // The null history lies in the first row
  vector<int> nullhist = true_hists.row(0);
  // Calculate the probabilities of event 0 in each primary period
  pvec.fill(1.0); // Initialize the vector
  for(j = 0; j < nper; j++){
    ppstart = pp_ind(j, 0);
    ppend   = pp_ind(j, 1);
    for(k = ppstart - 1; k < ppend; k++){
      pvec[j] *= pmat(k, nullhist[k]);
    }
  }
  // Probability that an individual entered the population in the last period and was not captured
  Type nullprob = beta[nper-1] * pvec[nper-1];
  // Calculate the null probability iteratively
  for(i = 0; i < nper - 1; i++){
    // An individual entered the population in period i+1 and was not captured
    tmpprob = beta[i] * pvec[i];
    // Deal with 0 entries at the end of the history
    tmp = 1 - phi[nper-2] + phi[nper-2] * pvec[nper-1];
    // If there are two or more 0 entries:
    if(i < nper - 2){
      for(k = nper - 2; k > i; k--) tmp = tmp * phi[k-1] * pvec[k] + 1 - phi[k-1];
    }
    tmpprob *= tmp;
    nullprob += tmpprob;
  }
  true_hists_probs[0] = nullprob;
  
  // Now calculate the probability of each true capture history
  for(i = 1; i < nhist; i++){
    true_hist = true_hists.row(i);
    // Probability of each event (0 or 1) in the corresponding true primary history
    pvec.fill(1.0);
    for(j = 0; j < nper; j++){
      ppstart = pp_ind(j, 0);
      ppend   = pp_ind(j, 1);
      for(k = ppstart - 1; k < ppend; k++){
        pvec[j] *= pmat(k, true_hist[k]);
      }
    }
    // REPORT(pvec);
    // Probability of all the 0 entries before the first 1 entry
    if(eventone_ind(i, 0) > 1) {
      // Only one 0 entry
      tmp = beta[eventone_ind(i, 0) - 1];
      tmp += beta[0] * pvec.head(eventone_ind(i, 0) - 1).prod() * phi.head(eventone_ind(i, 0) - 1).prod();
      // Two or more 0 entries
      if(eventone_ind(i, 0) > 2) {
        for(j = 1; j < eventone_ind(i, 0) - 1; j++) {
          tmp += beta[j] * (pvec.segment(j, eventone_ind(i, 0) - j -1)).prod() * (phi.segment(j, eventone_ind(i, 0) - j - 1)).prod();
        }
      }
      true_hists_probs[i] *= tmp;
    } else {
      true_hists_probs[i] *= beta[0];
    }
    // Then consider the strand between the first and last entry 1 of the primary history
    for(j = eventone_ind(i, 0) - 1; j < eventone_ind(i, 1); j++) {
      // Capture / non-capture probability
      true_hists_probs[i] *= pvec[j];
    }
    if(eventone_ind(i, 1) > eventone_ind(i, 0)) {
      // Survival between two 1 entries
      for(k = eventone_ind(i, 0) - 1; k < eventone_ind(i, 1) - 1; k++){
        true_hists_probs[i] *= phi[k];
      }
    }
    // Lastly consider 0 entries (if any) at the end of the primary history
    // e.g. the last two 0's of 10100
    if(eventone_ind(i, 1) < nper) {
      // Only one 0 entry at the end:
      tmp = 1 - phi[nper-2] + phi[nper-2] * pvec[nper-1];
      // If there are two or more 0 entries, calculate the probability iteratively:
      if(eventone_ind(i, 1) < nper - 1) {
        for(k = nper - 2; k > eventone_ind(i, 1) - 1; k--){
          tmp = tmp * phi[k-1] * pvec[k] + 1 - phi[k-1];
        }
      }
      true_hists_probs[i] *= tmp;
    }
  }
  REPORT(true_hists_probs);
  // Type sum_probs = true_hists_probs.sum();
  
  // Construct the cumulant generating function
  Type sum = 0.0;
  for(int i = 0; i < nhist; i++){
    sum += true_hists_probs[i] * exp((vector<Type> (linkmat.col(i)) * svec).sum());
  }
  Type CGF = N * log(sum);
  
  // Inner product of vectors obs_counts and svec
  Type obs_prod_svec = (obs_counts * svec).sum();
  
  // Compute penalty term
  Type pen = 0.0;
  if(penalty == 1){
    //pt
    for(j=0; j < nocca; j++){
      pen -= p_lambda[j] * pt[j] * pt[j];
    }
    //beta0t
    for(j=0; j < (nper - 1); j++){
      pen -= beta0_lambda[j] * beta0t[j] * beta0t[j];
    }
    //phit
    for(j=0; j< (nper-1); j++){
      pen -= phi_lambda[j] * phit[j] * phit[j];
    }
    //pen -= S_lambda * St * St;
  }

  // Our objective is CGF - marg_prod_svec:
  Type obj = CGF - obs_prod_svec + pen;
  
  

  vector<Type> surv = invlogit(phit);
  // ---------------------------------------------------------------------------
  // Number alive in each primary period: N_alive[k] = number alive in period k+1
  // (length nper; no confounding here, so the final period is identifiable too)
  //   N_alive[0] = N * beta_1
  //   N_alive[k] = phi_k * N_alive[k-1] + N * beta_{k+1}
  // ---------------------------------------------------------------------------
  vector<Type> N_alive(nper);
  N_alive[0] = N * beta[0];
  for (int k = 1; k < nper; k++) {
    N_alive[k] = surv[k - 1] * N_alive[k - 1] + N * beta[k];
  }
  
  
  // REPORT(obj);
  // REPORT(CGF);
  // REPORT(obs_prod_svec);
  // REPORT(pen);
  ADREPORT(Nt);
  ADREPORT(pt);
  ADREPORT(phit);
  ADREPORT(betat);
  
  ADREPORT(N_alive);
  
  
  return obj;
} // End of template
