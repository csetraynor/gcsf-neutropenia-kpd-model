// Created: 2020-01-06 15:13:31 UTC
// Author: Carlos Serra Traynor
$GLOBAL

double HillEQ(double x, double k, double tau) {
  double a = pow(std::max(x,0.0),k);
  return a/(pow(tau,k) + a);
}

double kmEQ(double VMAX, double KM, double CP) {
  double a = std::max(CP,0.0);
  return VMAX/(KM + a);
}

$PARAM @annotated

  // Docetaxel KDE
  KDE   : 5 : KDE Paclitaxe
  VDP  : 1   : Virtual Volume

  TVCIR       :  5.26     : Initial CIRC
  TVMTT       :  4.208    : MTT (days)
  TVHALF      :  0.2916   : Half-life blood (days)
  TVGAM       :  0.12     : Feeddback proliferation
  TVBET       :  0.08     : Feedback maturation   
  
  // Model of chemotherapy-induced myelosuppression with parameter consistency across drugs
  ALPHA       : 0.01  : Slope CBDCA
  
  ntrans  : 3 : n transit compartment

$OMEGA
  0
  0
  0
  0

$SIGMA
  0

$CMT

  A1DP // Docetaxel virtual compartment
    
  prol  // Proliferative
  transit1  // Transit compartment
  transit2  // Transit compartment
  transit3  // Transit compartment
  circ  // Circulation

$MAIN
  
  double VIRDP =  A1DP / VDP; //*KDEPTX;
  
  // pd model
  double CIR0  = TVCIR * exp(ETA(1));
  double MTT  = TVMTT * exp(ETA(2));
  double GAM   = TVGAM * exp(ETA(3));
  double ktr   = ( (double)ntrans + 1.0 ) / MTT;
  
  double BET   = TVBET;
  double ken   = log(2)/TVHALF;
  // double BET = 0.0;
  // double ken = ktr;
  
  prol_0       = (ken *CIR0) / ktr ;
  transit1_0   = (ken *CIR0) / ktr ;
  transit2_0   = (ken *CIR0) / ktr ;
  transit3_0   = (ken *CIR0) / ktr ;
  
  circ_0 = CIR0; 
  
$ODE
  
  double CIRCp = std::max( circ, 1e-15 );
  double BASNp = std::max( circ, 1e-15 );

  // double FBMN = pow((CIR0 - CIRCp)/CIR0, BET);
  double FB   = pow(CIR0/CIRCp, GAM);
  
  // double FN   = GAM * ((CIR0-CIRCp)/CIR0);
  // double EDrug = 1 / (1 + ALPHACOMBO*VIRPTX + ALPHACOMBO*VIRCBDCA);
  
  double Edrug = ALPHA*VIRDP;
  
  // K model
  dxdt_A1DP = - KDE * A1DP;
  
  // PD model
  dxdt_prol     = ktr * prol * FB * (1.0 - Edrug) - ktr * prol;
  dxdt_transit1 = ktr * (prol - transit1);
  dxdt_transit2 = ktr * (transit1 - transit2);
  dxdt_transit3 = ktr * (transit2 - transit3);
  dxdt_circ     = ktr * transit3 - ken*CIRCp;
  // dxdt_circ     = (ktr / MTT) * prol - ken * CIRCp;

$TABLE

  double VIRdp   =  A1DP/VDP;
  double edrug = ALPHA*VIRDP;
  double circp = std::max( circ, 1e-15 );
  double fn = GAM * ((CIR0 - circp)/CIR0);
  double ANC = std::max( circ, 1e-16 );

$CAPTURE
  VIRdp edrug ANC fn
  