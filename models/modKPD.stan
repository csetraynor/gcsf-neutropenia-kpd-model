// --------------------------------------------------------------------------
// Kinetic-Pharmacodynamic (KPD) Model for Docetaxel-Induced Neutropenia
// --------------------------------------------------------------------------
// Friberg-style 3-transit-compartment myelosuppression model with virtual
// kinetic input representing docetaxel-driven proliferation inhibition.
// Fitted using Torsten Stan (pmx_solve_group_rk45).
//
// Reference: Serra Traynor C, Kekic M, Boulton D, Zhou D.
//   PAGE 34 (2026) Abstr 12119.
// --------------------------------------------------------------------------

functions {
  
  /**
  * Raise each element of x to the power of y
  *
  * @param x Vector
  * @param y Real, the power to raise to
  * @return vector
  */
  vector pow_vec(vector x, real y) {
    int N = rows(x);
    vector[N] res;
    for (n in 1:N)
      res[n] = pow(x[n], y);
    return res;
  }
  
  /**
  * Apply Box-Cox transformation to each element of x with lambda y
  *
  * @param x Vector
  * @param y Real, the power to raise to
  * @return vector
  */
  vector boxcox_vec(vector x, real y) {
    int N = rows(x);
    vector[N] res;
    for (n in 1:N)
      res[n] = (pow(x[n], y) - 1.0) / y;
    return res;
  }
  
  /**
  * Hyperbolic equation
  *
  * @param EMAX real
  * @param EC50 real
  * @param CP real
  * @return real
  */
  real kmEQ(real EMAX, real EC50, real CP) {
    real a = fmax(CP,0.0);
    return EMAX*(a/(EC50 + a));
  }
  
  /**
  * Ensures each element of x is > 0
  *
  * @param x Vector
  * @return vector
  */
  vector strictly_positive(vector x) {
    int N = rows(x);
    vector[N] res;
    for (n in 1 : N) {
      res[n] = fmax(x[n], machine_precision());
    }
    return res;
  }
  
  // define ODE system KPD
  vector ModelKPD(real t, vector x, array[] real parms,
                      array[] real rdummy, array[] int idummy) {
    // Return object (derivative)
    vector[6] y; // 
    
    real circ0         = parms[1];
    real slope         = parms[2];
    real MTT           = parms[3];
    real gamma         = parms[4];
    // real MTT           = 4.2;
    // real gamma         = 0.12;
    
    
    real t12blood = 7.0/24.0; 
    real ktr   = (4.0/MTT);
    real ken   = 2.377; // log(2)/t12blood;
    
    // KPD parameters
    // real KDE   = 0.5;
    real KDE   = parms[5];
    real VDP   = 1.0;
    
    // G-CSF parameters
    real kaG = 7.2;
    real keG = 0.72;
    
    // states
    real A1DP = x[1];
    
    // -----------PD ----------
    real prol0 = fmax(machine_precision(), (ken*circ0) / ktr);
    real prol     = x[2] + prol0;
    real transit1 = x[3] + prol0;
    real transit2 = x[4] + prol0;
    real transit3 = x[5] + prol0;
    real circ = fmax(machine_precision(), x[6] + circ0); // Device for implementing a modeled initial condition
    
    // ----------Drug effect ---
    real VIRDP   =  A1DP / VDP; // Virtual concentration
    
    real Edrug = slope*VIRDP;
    real FB = pow(circ0/circ, gamma); // Feedback 
    
    // Virtual Kinetic component 
    y[1] = - KDE * A1DP;
    
    // PD component of the ODE system
    y[2] = ktr * prol * FB * (1.0 - Edrug) - ktr * prol;
    y[3] = ktr * (prol - transit1);
    y[4] = ktr * (transit1 - transit2);
    y[5] = ktr * (transit2 - transit3);
    y[6] = ktr * transit3 - ken*circ;
    
    return y;
  }
  
}

data {
  int<lower=1> nId;
  int<lower=1> nt;
  int<lower=1> nObsPD; // number of observations
  array[nObsPD] int<lower=1> iObsPD; // index of observation
  
  int<lower=1> nIdGCSF;
  
  // NONMEM data
  array[nt] int<lower=1> cmt;
  array[nId] int<lower=1> start;
  array[nId] int<lower=1> end;
  
  array[nt] int evid;
  array[nt] real amt;
  array[nt] real time;
  array[nt] real<lower=0> rate;
  array[nt] real<lower=0> ii;
  array[nt] int<lower=0> addl;
  
  vector<lower=0>[nObsPD] PDObs; // observed PTX concentration (dependent variable)
  
  // Model WI priors 
  real circ0Prior;
  real gammaPrior;
  real mttPrior;
  real KDEPrior;
  
  real circ0PriorCV;
  real gammaPriorCV;
  real mttPriorCV;
  real KDEPriorCV;
  
  real ALPHAPrior;
  real ALPHAPriorCV;
  
  real priorSigmaPD;
}

transformed data {
  
  array[nt] int<lower=0> ss = rep_array(0, nt);
  int<lower=1> nCmt = 6;
  array[nt, nCmt] real biovar;
  array[nt, nCmt] real tlag;
  real lambda = 0.2;
  
  vector[nObsPD] logPDObs = log(PDObs);
  // vector[nObsPD] BCPDObs = boxcox_vec(PDObs, lambda);
  
  int<lower=1> nParam  = 5;
  int<lower=1> nRandom = 2;
  array[nId] int len;

  int nPred = 1;

  for (j in 1 : nId) {
    len[j] = end[j] - start[j] + 1;
  }
  
  for (j in 1 : nt) {
    for (i in 1 : nCmt) {
      biovar[j, i] = 1;
      tlag[j, i] = 0;
    }
  }
  
}

parameters {
  
  // model parameters
  real<lower=0> circ0Hat;
  real<lower=0> ALPHAHat;
  
  real<lower=0> gammaHat;
  real<lower=0> mttHat;
  real<lower=0> KDEHat;
  
  // real lambda;
  real<lower=0> sigmaPD;
  
  matrix[nRandom, nId] etaStd;
  cholesky_factor_corr[nRandom] L;
  vector<lower=0, upper=1>[nRandom] omega;
  
}

transformed parameters {
  // IIV
  // vector<lower=0>[nRandom] thetahat = to_vector({circ0Hat, mttHat, gammaHat, ALPHAHat});
  vector<lower=0>[nRandom] thetahat = to_vector({circ0Hat, ALPHAHat});
  matrix<lower=0>[nId, nRandom] theta;
  
  theta = (rep_matrix(thetahat, nId) .* exp(diag_pre_multiply(omega, L * etaStd)))';
}

model {
  
  vector[nObsPD] PDhatObs; // predicted PD
  vector[nt] PD; // predicted PD
  matrix[nCmt, nt] x;
  
  array[nId, nParam] real parms;
  
  array[nId] real mtt;
  array[nId] real circ0;
  array[nId] real gamma;
  array[nId] real KDE;
  array[nId] real slope;
  
  for (j in 1 : nId ) {
    circ0[j] = theta[j, 1];
    slope[j] = theta[j, 2];
    mtt[j]   = mttHat;
    gamma[j] = gammaHat;
    KDE[j] = KDEHat;
    
    parms[j,  : ] = {circ0[j], slope[j],  mtt[j], gamma[j], KDE[j]};
  }
  
  // print("GCSF: ", parmsGCSF);
   // print("rest: ", parms[(nIdGCSF+1):nId, :]);
  x = pmx_solve_group_rk45(ModelKPD, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, parms, 1e-6, 1e-6, 1e6);
  
  for (j in 1 : nId) {
    PD[start[j] : end[j]] = x[6, start[j] : end[j]]' + circ0[j];
  }
  
  PDhatObs = strictly_positive(PD[iObsPD]);
  
  // wip
  circ0Hat ~ lognormal(log(circ0Prior), circ0PriorCV);
  ALPHAHat ~ lognormal(log(ALPHAPrior), ALPHAPriorCV);
  mttHat ~ lognormal(log(mttPrior), mttPriorCV);
  gammaHat ~ lognormal(log(gammaPrior), gammaPriorCV);
  KDEHat ~ lognormal(log(KDEPrior), KDEPriorCV);
  
  // SIGMA
  sigmaPD ~ exponential(priorSigmaPD);
  
  // OMEGA
  omega ~ normal(0, 0.5);
  L ~ lkj_corr_cholesky(1);
  
  // Inter-individual variability
  to_vector(etaStd) ~ normal(0, 1);
  
  // observed data likelihood
  // BCPDObs ~ normal(boxcox_vec(PDhatObs, lambda), sigmaPD);
  logPDObs ~ normal(log(PDhatObs), sigmaPD);
  
}

generated quantities {
  
  vector[nId] mtt_IPRED;
  vector[nId] circ0_IPRED;
  vector[nId] gamma_IPRED;
  vector[nId] ALPHA_IPRED;
  vector[nId] KDE_IPRED;
 
  vector[nPred] mtt_PRED;
  vector[nPred] circ0_PRED;
  vector[nPred] gamma_PRED;
  vector[nPred] ALPHA_PRED;
  vector[nPred] KDE_PRED;
  
  // log likelihood
  vector[nObsPD] log_lik;

 //  // Variables for IIV
  matrix[nRandom, nPred] etaStdPred;
  matrix<lower=0>[nPred, nRandom] thetaPredM;
  corr_matrix[nRandom] rho;
  
 // 
  rho = L * L';

  for (i in 1 : nPred) {
    for (j in 1 : nRandom) {
      etaStdPred[j, i] = normal_rng(0, 1);
    }
  }


  thetaPredM = (rep_matrix(thetahat, 1) .* exp(diag_pre_multiply(omega, L * etaStdPred)))';
                

  for (j in 1 : nId) { // Fix issue should be from j in 1:nId
    circ0_IPRED[j]      = theta[j, 1];
    mtt_IPRED[j]        = mttHat;
    gamma_IPRED[j]      = gammaHat;
    ALPHA_IPRED[j]      = theta[j, 2];
    KDE_IPRED[j]        = KDEHat;
  }
  
  for (j in 1 : nPred) {
    circ0_PRED[j]      = thetaPredM[j, 1];
    mtt_PRED[j]        = mttHat;
    gamma_PRED[j]      = gammaHat;
    ALPHA_PRED[j]      = thetaPredM[j, 2];
    KDE_PRED[j]        = KDEHat;
    
  }
  
  // log - lik calculations
  {
    array[nId, nParam] real parms_IPRED;
    vector[nObsPD] PDhatObs_PRED; // predicted PD
    vector[nt] PD_PRED; // predicted PD
    matrix[nCmt, nt] x_PRED;
    
    for (j in 1 : nId ) {
      parms_IPRED[j,  : ] = {circ0_IPRED[j], ALPHA_IPRED[j], mtt_IPRED[j], gamma_IPRED[j], KDE_IPRED[j]};
    }
  
    x_PRED = pmx_solve_group_rk45(ModelKPD, nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, parms_IPRED, 1e-6, 1e-6, 1e6);
  
  
    for (j in 1 : nId) {
      PD_PRED[start[j] : end[j]] = x_PRED[6, start[j] : end[j]]' + circ0_IPRED[j];
   }
  
    PDhatObs_PRED = strictly_positive(PD_PRED[iObsPD]);
    
    for (n in 1:nObsPD) {
      log_lik[n] = normal_lpdf(logPDObs[n] | log(PDhatObs_PRED[n]),  sigmaPD);
    }
    
  }

}


