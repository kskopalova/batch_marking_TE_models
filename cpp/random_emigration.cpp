

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
  
  x(4,0) = 0;
  x(4,1) = 0;
  x(4,2) = 0;
  x(4,3) = 0;
  x(4,4) = 1;
    
  return x;
}

template<class Type>
matrix<Type> LastTwoPer(Type migProb,
                                 Type entryProb,
                                 Type survProb,
                                 Type captureProb,
                                 Type migProd) {
  matrix<Type> x(5, 5);
  
  for (int i = 1; i < 4; ++i) {
    x(i,0) = 0;
    x(i,1) = migProd * (1 - captureProb);
    x(i,2) = migProd * captureProb;
    x(i,3) = (1-migProd)/2;
    x(i,4) = (1-migProd)/2;
  }
  
  
  x(0,0) = 1;
  x(0,1) = entryProb * (1 - captureProb);
  x(0,2) = entryProb * captureProb;
  x(0,3) = 0;
  x(0,4) = 0;
  
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
  DATA_IVECTOR(primindend_par);
  DATA_IVECTOR(primindstart);  // indices of the first secondary occasion within each primary period
  
  PARAMETER(Nt); 
  PARAMETER_VECTOR(alphaprimt);
  PARAMETER_VECTOR(pt);
  PARAMETER_VECTOR(phiprimt);
  PARAMETER_VECTOR(gamma0t);
  PARAMETER(alphaprodt);
  PARAMETER_VECTOR(svec);
  
  DATA_VECTOR(phi_lambda);
  DATA_VECTOR(p_lambda);
  DATA_VECTOR(gamma0_lambda);
  DATA_VECTOR(a_lambda);
  DATA_SCALAR(alphaprod_lambda);
  DATA_INTEGER(penalty);
  
  
  Type N = Nmin + exp(Nt);
  Type alphaprod = exp(alphaprodt)/(1+exp(alphaprodt));

  //alpha (needs to be of length 'nocca' for calculations)
  vector<Type> alpha(nocca-1);
  alpha.fill(1);
  int b = 0;
  for (int index : primindend_par) {
    alpha[index - 1] = exp(alphaprimt[b])/(1+exp(alphaprimt[b]));
    b++;
  }

  // different p across all capture occasions
  vector<Type> p = exp(pt) / (1 + exp(pt));
  
  // different survival across the primary periods
  vector<Type> phi(nocca-1);
  phi.fill(1);
  int d = 0;
  for (int index : primindend_par) {
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
    int secondLastElement = primindend[primindend.size() - 1]-1;

    if (ArrivalState[i] == -1) {
      probs[i] = 0;
    } else {
      for (int j = 0; j < nocca - 1; j++) {

        if (j == secondLastElement) {
          Type survProbSub = phi[j];
          Type captureProbSub = p[j + 1];
          Type alphaprobsub = alpha[j];
          Type entryprobsub = gammafull[j+1];
          
          temp = LastTwoPer(alphaprobsub, entryprobsub, 
                                     survProbSub, captureProbSub, alphaprod);
          
          int temp_from = true_hist[j];
          int temp_to = true_hist[j + 1];
          
          probs[i] *= temp(temp_from, temp_to);
          
        } else if (std::find(primindend.begin(), primindend.end(), j+1) != primindend.end()) {
          Type survProbSub = phi[j];
          Type captureProbSub = p[j + 1];
          Type alphaprobsub = alpha[j];
          Type entryprobsub = gammafull[j+1];
          
          temp = BetweenPrimaryTrans(alphaprobsub, entryprobsub, 
                                     survProbSub, captureProbSub);
          
          int temp_from = true_hist[j];
          int temp_to = true_hist[j + 1];
          
          probs[i] *= temp(temp_from, temp_to);
        } 
        else {
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
    for(int j=0; j< (nper-2); j++){
      pen -= phi_lambda[j] * phiprimt[j] * phiprimt[j];
    }
    for(int j=0; j< (nper-2); j++){
      pen -= a_lambda[j] * alphaprimt[j] * alphaprimt[j];
    }
    pen -= alphaprod_lambda * alphaprodt * alphaprodt;
  }
  
  // Our objective is CGF - marg_prod_svec:
  Type obj = CGF - obs_prod_svec + pen;
    
  // product term is phi*(1-alpha)
  vector<Type> a = invlogit(alphaprimt);
  Type prod_term = invlogit(alphaprodt);
  vector<Type> surv = invlogit(phiprimt);
  
  // ---------------------------------------------------------------------------
  // Number alive: alive_all[k] = number alive in primary period k+1, k = 0..nper-2
  // (period nper cannot be obtained because phi_{nper-1} is confounded with alpha_{nper-1})
  // ---------------------------------------------------------------------------
  vector<Type> alive_all(nper - 1);
  alive_all[0] = N * gamma[0];
  for (int k = 1; k < nper - 1; k++) {
    alive_all[k] = surv[k - 1] * alive_all[k - 1] + N * gamma[k];
  }
  // N_alive[j] = number alive in primary period j+2 (length nper-2)
  vector<Type> N_alive(nper - 2);
  for (int j = 0; j < nper - 2; j++) {
    N_alive[j] = alive_all[j + 1];
  }

  // ---------------------------------------------------------------------------
  // Number available for detection: N_avail[k] = N_{k+1} (length nper)
  // Completely random emigration: N_{k+1} = phi_k * (1 - alpha_k) * (alive at k) + N * gamma_{k+1}
  // ---------------------------------------------------------------------------
  vector<Type> N_avail(nper);
  N_avail[0] = N * gamma[0];
  for (int k = 1; k < nper - 1; k++) {
    N_avail[k] = surv[k - 1] * (Type(1.0) - a[k - 1]) * alive_all[k - 1]
               + N * gamma[k];
  }
  // final transition uses the confounded product phi_{nper-1} * (1 - alpha_{nper-1})
  N_avail[nper - 1] = prod_term * alive_all[nper - 2]
                    + N * gamma[nper - 1];
  
  
  REPORT(obj);
  REPORT(CGF);
  REPORT(obs_prod_svec);
  ADREPORT(Nt);
  ADREPORT(pt);
  ADREPORT(phiprimt);
  ADREPORT(alphaprimt);
  ADREPORT(alphaprodt);
  ADREPORT(gammat);
  
  ADREPORT(N_avail);
  ADREPORT(N_alive);
  
  return obj;
  
} 
