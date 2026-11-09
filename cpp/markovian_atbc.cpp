
#include <TMB.hpp>

template<class Type>
matrix<Type> withinPrimaryTrans(Type captureProb) {
  
  // transitions within a primary occasion
  matrix<Type> x(5, 5);
  
  for (int i = 0; i < 5; ++i) {
    x(i,0) = (i == 0) ? 1 : 0;
    x(i,1) = 0;
    x(i,2) = 0;
    x(i,3) = (i == 3) ? 1 : 0;
    x(i,4) = (i == 4) ? 1 : 0; 
    
    if (i == 1) {
      x(i,1) = 1 - captureProb;
      x(i,2) = captureProb;
    }
    
    if (i == 2) {
      x(i,1) = 1 - captureProb;
      x(i,2) = captureProb;
    }
  }
  return x;
}

template<class Type>
matrix<Type> BetweenPrimaryTrans(Type migProb,
                                 Type returnProb,
                                 Type entryProb,
                                 Type survProb,
                                 Type captureProb) {
  
  // transitions between the primary occasions
  matrix<Type> x(5, 5);
  
  for (int i = 1; i < 4; ++i) {
    x(i,0) = 0;
    x(i,1) = survProb * (1 - migProb) * (1 - captureProb);
    x(i,2) = survProb * (1 - migProb) * captureProb;
    x(i,3) = survProb * migProb;
    x(i,4) = 1 - survProb; 
  }
  
  
  x(0,0) = 1;
  x(0,1) = entryProb * (1 - captureProb);
  x(0,2) = entryProb * captureProb;
  x(0,3) = 0;
  x(0,4) = 0;
  
  x(3,1) = survProb * returnProb * (1 - captureProb);
  x(3,2) = survProb * returnProb * captureProb;
  x(3,3) = survProb * (1 - returnProb);
  
  x(4,0) = 0;
  x(4,1) = 0;
  x(4,2) = 0;
  x(4,3) = 0;
  x(4,4) = 1;
  
  return x;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  
  // Data objects
  DATA_SPARSE_MATRIX(linkmat);
  DATA_VECTOR(obs_counts);
  DATA_IMATRIX(true_hists);
  DATA_SCALAR(Nmin);
  DATA_INTEGER(nocca); 
  DATA_INTEGER(nper);
  DATA_INTEGER(nhist);
  DATA_IVECTOR(primindend);    // indices of the last secondary occasion within each primary period
  DATA_IVECTOR(subprimindend); // indices of the last secondary occasion within each primary period (excluding the first one)
  DATA_IVECTOR(primindstart);  // indices of the first secondary occasion within each primary period
  
  PARAMETER(Nt); 
  PARAMETER_VECTOR(alphaprimt); 
  PARAMETER(betat); 
  PARAMETER_VECTOR(pt);
  PARAMETER_VECTOR(phiprimt);
  PARAMETER_VECTOR(gamma0t);
  PARAMETER_VECTOR(svec);
  
  DATA_VECTOR(phi_lambda);
  DATA_VECTOR(p_lambda);
  DATA_VECTOR(gamma0_lambda);
  DATA_VECTOR(a_lambda);
  DATA_SCALAR(b_lambda);
  DATA_INTEGER(penalty);
  
  Type N = Nmin + exp(Nt);
  vector<Type> alpha(nocca);
  alpha.fill(1);
  int b = 0;
  for (int index : primindend) {
    alpha[index - 1] = exp(alphaprimt[b])/(1+exp(alphaprimt[b]));
    b++;
  }
  
  vector<Type> beta(nocca);
  beta.fill(1);
  for (int index : subprimindend) {
    beta[index - 1] = exp(betat)/(1+exp(betat));
  }
  
  // different p across all capture occasions
  vector<Type> p = exp(pt) / (1 + exp(pt));
  
  // different survival across the primary periods
  vector<Type> phi(nocca);
  phi.fill(1);
  int d = 0;
  for (int index : primindend) {
    phi[index - 1] = exp(phiprimt[d])/(1+exp(phiprimt[d]));
    d++;
  }
  
  Type gamma0 = 0.0;
  Type cumsum = 0.0;
  vector<Type> gamma(nper);
  vector<Type> gammat(nper);
  
  for(int j=0;j < nper; j++){
    // Compute probability of entry conditional on not having entered yet
    if(j == nper - 1){
      gamma0 = 1;
    }
    else{
      gamma0 = exp(gamma0t[j]) / (1 + exp(gamma0t[j]));
    }
    // Compute marginal probability of entry
    gamma[j] = gamma0 * (1 - cumsum);
    gammat[j] = log(gamma[j]) - log(1 - gamma[j]);
    // Update cumulative sum of entry probabilities
    cumsum += gamma[j];
  }
  
  // gamma vector used in calculations
  vector<Type> gammafull(nocca);
  gammafull.fill(1);
  int ii = 0;
  for (int index : primindstart) {
    gammafull[index - 1] = gamma[ii];
    ii++;
  }
  
  vector<int> true_hist(nocca);
  vector<int> ArrivalInd(nhist);
  ArrivalInd.fill(-1); 
  
  for (int i = 0; i < nhist; i++) {
    true_hist = true_hists.row(i); 
    
    for (int j = 0; j < nocca; j++) {
      if (true_hist[j] == 1 || true_hist[j] == 2) {
        ArrivalInd[i] = j;
        break;
      } 
    }
  }
  
  vector<int> ArrivalState(nhist);
  vector<Type> probs(nhist);
  probs.fill(1);
  matrix<Type> temp(5, 5);
  
  // Iterate over each history
  for (int i = 0; i < nhist; i++) {
    true_hist = true_hists.row(i);
    ArrivalState[i] = true_hist[ArrivalInd[i]];
    
    if (ArrivalState[i] == -1) {
      probs[i] = 0;
    } else {
      for (int j = 0; j < nocca - 1; j++) {
        
        
        if (std::find(primindend.begin(), primindend.end(), j+1) != primindend.end()) {
          Type survProbSub = phi[j];
          Type captureProbSub = p[j + 1];
          Type alphaprobsub = alpha[j];
          Type betaprobsub = beta[j];
          Type entryprobsub = gammafull[j+1];
          
          temp = BetweenPrimaryTrans(alphaprobsub, betaprobsub, entryprobsub, 
                                     survProbSub, captureProbSub);
          
          int temp_from = true_hist[j];
          int temp_to = true_hist[j + 1];
          
          probs[i] *= temp(temp_from, temp_to);
          
        } else {
          Type captureProbSub = p[j + 1];
          temp = withinPrimaryTrans(captureProbSub);
          int temp_from = true_hist[j];
          int temp_to = true_hist[j + 1];
          
          probs[i] *= temp(temp_from, temp_to);
        }
      }}
    
    
    bool S1 = (true_hist[0] == 1);
    bool S2 = (true_hist[0] == 2);
    bool S3 = (true_hist[0] == 3);
    bool S4 = (true_hist[0] == 4);
    
    if (S1 || S2) {
      
      Type probarrival = (true_hist[0] == 1) ? gammafull[0] * (1 - p[0]) : gammafull[0] * p[0];
      probs[i] *= probarrival;
    } else if (S3 || S4) {
      probs[i] = 0;
    }
  }
  
  
  // ---------------------------------------------------------------------------------------------------  
  REPORT(probs); 
  Type sum_probs = probs.sum();
  //Construct the cumulant generating function
  Type sum = 0.0;
  for(int i = 0; i < nhist; i++){
    sum += probs[i] * exp((vector<Type> (linkmat.col(i)) * svec).sum());
  }
  REPORT(sum);
  REPORT(N);
  REPORT(Nt);
  Type CGF = N * log(sum);
  
  // Inner product of vectors obs_counts and svec
  Type obs_prod_svec = (obs_counts * svec).sum();
  
  //Compute penalty term
  Type pen = 0.0;
  if(penalty == 1){
    //pt
    for(int j=0; j < nocca; j++){
      pen -= p_lambda[j] * pt[j] * pt[j];
    }
    //beta0t
    for(int j=0; j < (nper - 1); j++){
      pen -= gamma0_lambda[j] * gamma0t[j] * gamma0t[j];
    }
    //phit
    for(int j=0; j< (nper-1); j++){
      pen -= phi_lambda[j] * phiprimt[j] * phiprimt[j];
    }
    
    for(int j=0; j< (nper-1); j++){
      pen -= a_lambda[j] * alphaprimt[j] * alphaprimt[j];
    }
    
    pen -= b_lambda * betat * betat;
    
  }
  
  Type obj = CGF - obs_prod_svec + pen;
  REPORT(obj);
  REPORT(CGF);
  REPORT(obs_prod_svec);
  
  vector<Type> a = invlogit(alphaprimt);
  Type ret_prob = invlogit(betat);
  vector<Type> surv = invlogit(phiprimt);
  
  // ---------------------------------------------------------------------------
  // Number available for detection: N_avail[k] = N_{k+1} (length nper)
  // Recursion (U = alive but outside the study area, U_1 = 0):
  //   N_{k+1} = phi_k * [(1 - alpha_k) * N_k + beta_k * U_k] + N * gamma_{k+1}
  //   U_{k+1} = phi_k * [alpha_k * N_k + (1 - beta_k) * U_k]
  // ---------------------------------------------------------------------------
  vector<Type> N_avail(nper);
  N_avail[0] = N * gamma[0];
  Type U_cur = Type(0.0);                 // U_1 = 0 (entrants are available on entry)
  for (int k = 1; k <= nper - 1; k++) {      // transition from primary period k to k+1
    Type alpha_k = a[k - 1];
    Type beta_k  = ret_prob;
    Type N_next = surv[k - 1] * ((Type(1.0) - alpha_k) * N_avail[k - 1] + beta_k * U_cur)
                + N * gamma[k];
    Type U_next = surv[k - 1] * (alpha_k * N_avail[k - 1] + (Type(1.0) - beta_k) * U_cur);
    N_avail[k] = N_next;
    U_cur = U_next;
  }

  // ---------------------------------------------------------------------------
  // Number alive: N_alive[j] = number alive in primary period j+2 (period 1 equals N_avail[0])
  // length nper-1 without confounding, nper-2 with confounding (final period not identifiable)
  // ---------------------------------------------------------------------------
  vector<Type> N_alive(nper - 1);
  Type alive_prev = N_avail[0];
  for (int j = 0; j < nper - 1; j++) {
    alive_prev = surv[j] * alive_prev + N * gamma[j + 1];
    N_alive[j] = alive_prev;
  }
  
  
  ADREPORT(Nt);
  ADREPORT(pt);
  ADREPORT(phiprimt);
  ADREPORT(alphaprimt);
  ADREPORT(betat);
  ADREPORT(gammat); 
  
  ADREPORT(N_avail);
  ADREPORT(N_alive);
  


  
  return obj;
  
} 
