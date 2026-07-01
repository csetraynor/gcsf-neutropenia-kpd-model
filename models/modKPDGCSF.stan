// --------------------------------------------------------------------------
// Kinetic-Pharmacodynamic (KPD) Model with G-CSF for Docetaxel-Induced
// Neutropenia
// --------------------------------------------------------------------------
// Friberg-style 3-transit-compartment myelosuppression model augmented with
// G-CSF stimulation of neutrophil proliferation and maturation.
// Virtual kinetic input represents docetaxel-driven proliferation inhibition.
// Patients receiving G-CSF are solved with ModelGCSF; others with ModelKPD.
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
    vector[8] y; // 
    
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
    real GSQ  = x[2];
    real G    = x[3];
    
    // -----------PD ----------
    real prol0 = fmax(machine_precision(), (ken*circ0) / ktr);
    real prol     = x[4] + prol0;
    real transit1 = x[5] + prol0;
    real transit2 = x[6] + prol0;
    real transit3 = x[7] + prol0;
    real circ = fmax(machine_precision(), x[8] + circ0); // Device for implementing a modeled initial condition
    
    // ----------Drug effect ---
    real VIRDP   =  A1DP / VDP; // Virtual concentration
    
    real Edrug = slope*VIRDP;
    real FB = pow(circ0/circ, gamma); // Feedback 
    
    // Virtual Kinetic component 
    y[1] = - KDE * A1DP;
    
      // Dummy G-CSF KPD
    y[2] = -kaG * GSQ;
    y[3] =  kaG * GSQ - keG * G;
    
    // PD component of the ODE system
    y[4] = ktr * prol * FB * (1.0 - Edrug) - ktr * prol;
    y[5] = ktr * (prol - transit1);
    y[6] = ktr * (transit1 - transit2);
    y[7] = ktr * (transit2 - transit3);
    y[8] = ktr * transit3 - ken*circ;
    
    return y;
  }
  
  // define ODE system with G-CSF effect
  vector ModelGCSF(real t, vector x, array[] real parms,
                      array[] real rdummy, array[] int idummy) {
    // Return object (derivative)
    vector[8] y; // 
    
    real circ0         = parms[1];
    real slope         = parms[2];
    real MTT           = parms[3];
    real gamma         = parms[4];
    // real MTT           = 4.2;
    // real gamma         = 0.12;
    
    // Fixed parameters 
    real t12blood = 7.0/24.0; 
    real ktr0   = (4.0/MTT);
    real ken   = 2.377; // log(2)/t12blood;
    
    // KPD parameters
    // real KDE   = 0.5;
    real KDE   = parms[5];
    real VDP   = 1.0;
    
    // G-CSF parameters
    real Emax_prol     = parms[6];
    // real Emax_prol     = 1.0;
    real Emax_mat      = 1.0;
    real kaG = 7.2;
    real keG = parms[7];
    // real keG = 0.72;
        
    real VPDG = 1.0;
    real EC50_prol = 0.5;
    real EC50_mat  = 0.5;
    
    // states
    real A1DP = x[1];
    real GSQ  = x[2];
    real G    = x[3];
    
    // Effective transit rate (maturation stimulation)
    real CG   = G / VPDG;
    real Eprol = kmEQ(Emax_prol, EC50_prol, CG); 
    real Emat  = kmEQ(Emax_mat, EC50_mat, CG); 
    
    real ktr = ktr0 * (1.0 + Emat);
    
    // -----------PD ----------
    real prol0 = fmax(machine_precision(), (ken*circ0) / ktr0);
    real prol     = x[4] + prol0;
    real transit1 = x[5] + prol0;
    real transit2 = x[6] + prol0;
    real transit3 = x[7] + prol0;
    real circ = fmax(machine_precision(), x[8] + circ0); // Device for implementing a modeled initial condition
    
    // ----------Drug effect ---
    real VIRDP   =  A1DP / VDP; // Virtual concentration
    
    real Edrug = slope*VIRDP;
    real FB = pow(circ0/circ, gamma);
    
    // Virtual Kinetic component 
    y[1] = - KDE * A1DP;
    
      // Dummy G-CSF KPD
    y[2] = -kaG * GSQ;
    y[3] =  kaG * GSQ - keG * G;
    
    // PD component of the ODE system
    y[4] = ktr * prol * FB * (1.0 - Edrug) * (1.0 + Eprol) - ktr * prol;
    y[5] = ktr * (prol - transit1);
    y[6] = ktr * (transit1 - transit2);
    y[7] = ktr * (transit2 - transit3);
    y[8] = ktr * transit3 - ken*circ;
    
    return y;
  }
}

data {
  int<lower=1> nId;
  int<lower=1> nt;
  int<lower=1> nObsPD; // number of observations
  array[nObsPD] int<lower=1> iObsPD; // index of observation
  
  int<lower=1> nIdGCSF;
  int<lower=1> ntGCSF;
  
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
  
  real Emax_prolPrior;
  real keG_Prior;

  real Emax_prolPriorCV;
  real keG_PriorCV;
    
  real priorSigmaPD;
}

transformed data {
  
  array[nt] int<lower=0> ss = rep_array(0, nt);
  int<lower=1> nCmt = 8;
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
  
  real<lower=0> Emax_prolHat;
  real<lower=0> keGHat;
  
  // real lambda;
  real<lower=0> sigmaPD;
  
  matrix[nRandom, nId] etaStd;
  cholesky_factor_corr[nRandom] L;
  vector<lower=0, upper=1>[nRandom] omega;
  
  matrix[2, nIdGCSF] etaStdGCSF;
  cholesky_factor_corr[2] LGCSF;
  vector<lower=0, upper=1>[2] omegaGCSF;
  
}

transformed parameters {
  // IIV
  // vector<lower=0>[nRandom] thetahat = to_vector({circ0Hat, mttHat, gammaHat, ALPHAHat});
  vector<lower=0>[nRandom] thetahat = to_vector({circ0Hat, ALPHAHat});
  matrix<lower=0>[nId, nRandom] theta;
  
  vector<lower=0>[2] thetahatGCSF = to_vector({Emax_prolHat, keGHat});
  matrix<lower=0>[nIdGCSF, 2] thetaGCSF;
  
  theta = (rep_matrix(thetahat, nId) .* exp(diag_pre_multiply(omega, L * etaStd)))';
  thetaGCSF = (rep_matrix(thetahatGCSF, nIdGCSF) .* exp(diag_pre_multiply(omegaGCSF, LGCSF * etaStdGCSF)))';
}

model {
  
  vector[nObsPD] PDhatObs; // predicted PD
  vector[nt] PD; // predicted PD
  matrix[nCmt, nt] x;
  
  array[nIdGCSF, nParam+2] real parmsGCSF;
  array[nId, nParam] real parms;
  
  array[nId] real mtt;
  array[nId] real circ0;
  array[nId] real gamma;
  array[nId] real KDE;
  array[nId] real slope;
  
  array[nId] real Emax_prol;
  array[nId] real keG;
  
  for (j in 1 : nIdGCSF) {
    circ0[j] = theta[j, 1];
    slope[j] = theta[j, 2];
    mtt[j]   = mttHat;
    gamma[j] = gammaHat;
    KDE[j] = KDEHat;
    
    Emax_prol[j] = thetaGCSF[j, 1];
    keG[j]       = thetaGCSF[j, 2];
    
    parmsGCSF[j,  : ] = {circ0[j], slope[j], mtt[j], gamma[j], KDE[j], Emax_prol[j], keG[j]};
  }
  
  for (j in 1 : nId ) {
    circ0[j] = theta[j, 1];
    slope[j] = theta[j, 2];
    mtt[j]   = mttHat;
    gamma[j] = gammaHat;
    KDE[j] = KDEHat;
    
    parms[j,  : ] = {circ0[j], slope[j],  mtt[j], gamma[j], KDE[j]};
  }
  
  // print("GCSF: ", parmsGCSF);
  x[:, 1:ntGCSF] = pmx_solve_group_rk45(ModelGCSF, nCmt, len[1:nIdGCSF], time[1:ntGCSF], amt[1:ntGCSF], rate[1:ntGCSF], ii[1:ntGCSF], evid[1:ntGCSF], cmt[1:ntGCSF], addl[1:ntGCSF], ss[1:ntGCSF], parmsGCSF, 1e-6, 1e-6, 1e6);
  
   // print("rest: ", parms[(nIdGCSF+1):nId, :]);
  x[:, (ntGCSF+1):nt ] = pmx_solve_group_rk45(ModelKPD, nCmt, len[(nIdGCSF+1):nId], time[(ntGCSF+1):nt], amt[(ntGCSF+1):nt], rate[(ntGCSF+1):nt], ii[(ntGCSF+1):nt], evid[(ntGCSF+1):nt], cmt[(ntGCSF+1):nt], addl[(ntGCSF+1):nt], ss[(ntGCSF+1):nt], parms[(nIdGCSF+1):nId, :], 1e-6, 1e-6, 1e6);
  
  for (j in 1 : nId) {
    PD[start[j] : end[j]] = x[8, start[j] : end[j]]' + circ0[j];
  }
  
  PDhatObs = strictly_positive(PD[iObsPD]);
  
  // wip
  circ0Hat ~ lognormal(log(circ0Prior), circ0PriorCV);
  ALPHAHat ~ lognormal(log(ALPHAPrior), ALPHAPriorCV);
  mttHat ~ lognormal(log(mttPrior), mttPriorCV);
  gammaHat ~ lognormal(log(gammaPrior), gammaPriorCV);
  KDEHat ~ lognormal(log(KDEPrior), KDEPriorCV);
  
  Emax_prolHat ~ lognormal(log(Emax_prolPrior), Emax_prolPriorCV);
  keGHat ~ lognormal(log(keG_Prior), keG_PriorCV);

  // SIGMA
  sigmaPD ~ exponential(priorSigmaPD);
  
  // OMEGA
  omega ~ normal(0, 0.5);
  L ~ lkj_corr_cholesky(1);
  
  // Inter-individual variability
  to_vector(etaStd) ~ normal(0, 1);
  
   // OMEGA
  omegaGCSF ~ normal(0, 0.5);
  LGCSF ~ lkj_corr_cholesky(1);
  // 
  // // Inter-individual variability
  to_vector(etaStdGCSF) ~ normal(0, 1);
  
  
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
  vector[nIdGCSF] Emax_prol_IPRED;
  vector[nIdGCSF] keG_IPRED;
 
  vector[nPred] mtt_PRED;
  vector[nPred] circ0_PRED;
  vector[nPred] gamma_PRED;
  vector[nPred] ALPHA_PRED;
  vector[nPred] KDE_PRED;
  vector[nPred] Emax_prol_PRED;
  vector[nPred] keG_PRED;

 //  // Variables for IIV
  matrix[nRandom, nPred] etaStdPred;
  matrix[2, nPred] etaStdPredGCSF;
  matrix<lower=0>[nPred, nRandom] thetaPredM;
  matrix<lower=0>[nPred, 2] thetaPredMGCSF;
  corr_matrix[nRandom] rho;
  corr_matrix[2] rhoGCSF;
  
  // log likelihood
  vector[nObsPD] log_lik;
  
 // 
  rho = L * L';
  rhoGCSF = LGCSF * LGCSF';

  for (i in 1 : nPred) {
    for (j in 1 : nRandom) {
      etaStdPred[j, i] = normal_rng(0, 1);
    }
  }

  for (i in 1 : nPred) {
    for (j in 1 : 2) {
      etaStdPredGCSF[j, i] = normal_rng(0, 1);
    }
  }

  thetaPredM = (rep_matrix(thetahat, 1) .* exp(diag_pre_multiply(omega, L * etaStdPred)))';
                
  thetaPredMGCSF = (rep_matrix(thetahatGCSF, 1) .* exp(diag_pre_multiply(omegaGCSF, LGCSF * etaStdPredGCSF)))';

  for (j in 1 : nId) { // Fix issue should be from j in 1:nId
    circ0_IPRED[j]      = theta[j, 1];
    mtt_IPRED[j]        = mttHat;
    gamma_IPRED[j]      = gammaHat;
    ALPHA_IPRED[j]      = theta[j, 2];
    KDE_IPRED[j]        = KDEHat;
  }
  
  for (j in 1 : nIdGCSF) { // Fix issue should be from j in 1:nId
    Emax_prol_IPRED[j]  = thetaGCSF[j, 1];
    keG_IPRED[j]   = thetaGCSF[j, 2];
  }
  
  for (j in 1 : nPred) {
    circ0_PRED[j]      = thetaPredM[j, 1];
    mtt_PRED[j]        = mttHat;
    gamma_PRED[j]      = gammaHat;
    ALPHA_PRED[j]      = thetaPredM[j, 2];
    KDE_PRED[j]        = KDEHat;
    
    Emax_prol_PRED[j]  = thetaPredMGCSF[j, 1];
    keG_PRED[j]   = thetaPredMGCSF[j, 2];
  }
  
  
  // log - lik calculations
  {
    array[nIdGCSF, nParam+2] real parmsGCSF_IPRED;
    array[nId, nParam] real parms_IPRED;
    vector[nObsPD] PDhatObs_PRED; // predicted PD
    vector[nt] PD_PRED; // predicted PD
    matrix[nCmt, nt] x_PRED;
    
    for (j in 1 : nIdGCSF) {
      parmsGCSF_IPRED[j,  : ] = {circ0_IPRED[j], ALPHA_IPRED[j], mtt_IPRED[j], gamma_IPRED[j], KDE_IPRED[j], Emax_prol_IPRED[j], keG_IPRED[j]};
    }
  
    for (j in 1 : nId ) {
      parms_IPRED[j,  : ] = {circ0_IPRED[j], ALPHA_IPRED[j], mtt_IPRED[j], gamma_IPRED[j], KDE_IPRED[j]};
    }
  
    x_PRED[:, 1:ntGCSF] = pmx_solve_group_rk45(ModelGCSF, nCmt, len[1:nIdGCSF], time[1:ntGCSF], amt[1:ntGCSF], rate[1:ntGCSF], ii[1:ntGCSF], evid[1:ntGCSF], cmt[1:ntGCSF], addl[1:ntGCSF], ss[1:ntGCSF], parmsGCSF_IPRED, 1e-6, 1e-6, 1e6);
  
    x_PRED[:, (ntGCSF+1):nt ] = pmx_solve_group_rk45(ModelKPD, nCmt, len[(nIdGCSF+1):nId], time[(ntGCSF+1):nt], amt[(ntGCSF+1):nt], rate[(ntGCSF+1):nt], ii[(ntGCSF+1):nt], evid[(ntGCSF+1):nt], cmt[(ntGCSF+1):nt], addl[(ntGCSF+1):nt], ss[(ntGCSF+1):nt], parms_IPRED[(nIdGCSF+1):nId, :], 1e-6, 1e-6, 1e6);
  
    for (j in 1 : nId) {
      PD_PRED[start[j] : end[j]] = x_PRED[8, start[j] : end[j]]' + circ0_IPRED[j];
   }
  
    PDhatObs_PRED = strictly_positive(PD_PRED[iObsPD]);
    
    for (n in 1:nObsPD) {
      log_lik[n] = normal_lpdf(logPDObs[n] | log(PDhatObs_PRED[n]),  sigmaPD);
    }
    
  }

}


