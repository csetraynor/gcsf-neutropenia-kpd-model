// --------------------------------------------------------------------------
// Kinetic-Pharmacodynamic (KPD) Model with G-CSF for Docetaxel-Induced
// Neutropenia
// --------------------------------------------------------------------------
// Author: Carlos Serra Traynor
$GLOBAL

double HillEQ(double x, double k, double tau) {
  double a = pow(std::max(x,0.0),k);
  return a/(pow(tau,k) + a);
}

double kmEQ(double VMAX, double KM, double CP) {
  double a = std::max(CP,0.0);
  return VMAX*(a/(KM + a));
}

$PARAM @annotated

  // Docetaxel KDE
  KDE   : 0.5 : KDE Paclitaxe
  VDP  : 1   : Virtual Volume

  TVCIR       :  5.26     : Initial CIRC
  TVMTT       :  4.208    : MTT (days)
  TVHALF      :  0.2916   : Half-life blood (days)
  TVGAM       :  0.12     : Feeddback proliferation
  TVBET       :  0.08     : Feedback maturation   
  
  // Model of chemotherapy-induced myelosuppression with parameter consistency across drugs
  ALPHA       : 0.01  : Slope CBDCA
  
  ntrans  : 3 : n transit compartment
  
  // G-CSF KPD
  kaG   : 7.2     : absorption rate (day^-1) for SC G-CSF
  keG   : 0.72    : elimination rate (day^-1) for G-CSF KPD
  VPDG  : 1.0     : KPD 'volume' for virtual concentration (arbitrary units)
    
  // G-CSF PD parameters
  Emax_prol : 2.0     : max fold stimulation of proliferation
  EC50_prol : 0.5     : virtual conc for half-max proliferation stimulation
  Emax_mat  : 1.5     : max fold stimulation of maturation (ktr increase)
  EC50_mat  : 0.5     : virtual conc for half-max maturation stimulation
  
  USE_PROL    : 1         : GCSF effect on prol
  USE_MAT     : 1         : GCSF effect on ktr

$OMEGA
  0
  0
  0
  0

$SIGMA
  0

$CMT

  A1DP // Docetaxel virtual compartment
  
  GSQ  // G-CSF SUBCUT
  G    //G-CSF
    
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
  double ktr0   = ( (double)ntrans + 1.0 ) / MTT;
  
  // PD effects (saturable Emax)
  double CG   = G / VPDG;
  double Eprol = 0.0;
  double Emat  = 0.0;
  if(USE_PROL == 1) Eprol = kmEQ(Emax_prol, EC50_prol, CG); 
  if(USE_MAT  == 1) Emat  = kmEQ(Emax_mat, EC50_mat, CG); 
  
  // Effective transit rate (maturation stimulation)
  double ktr = ktr0 * (1.0 + Emat);
  
  double BET   = TVBET;
  double ken   = log(2)/TVHALF;
  // double BET = 0.0;
  // double ken = ktr;
  
  prol_0       = (ken *CIR0) / ktr0 ;
  transit1_0   = (ken *CIR0) / ktr0 ;
  transit2_0   = (ken *CIR0) / ktr0 ;
  transit3_0   = (ken *CIR0) / ktr0 ;
  
  circ_0 = CIR0; 
  
$ODE
  
  double CIRCp = std::max( circ, 1e-15 );
  // double BASNp = std::max( circ, 1e-15 );

  // double FBMN = pow((CIR0 - CIRCp)/CIR0, BET);
  double FB   = pow(CIR0/CIRCp, GAM);
  
  // double FN   = GAM * ((CIR0-CIRCp)/CIR0);
  // double EDrug = 1 / (1 + ALPHACOMBO*VIRPTX + ALPHACOMBO*VIRCBDCA);
  
  double Edrug = ALPHA*VIRDP;
  
  // K model
  dxdt_A1DP = - KDE * A1DP;
  
  // G-CSF KPD
  dxdt_GSQ = -kaG * GSQ;
  dxdt_G    =  kaG * GSQ - keG * G;
  
  // PD model
  dxdt_prol     = ktr * prol * FB * (1.0 - Edrug) * (1.0 + Eprol) - ktr * prol;
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
  double Vcg   = G / VPDG;

$CAPTURE
  VIRdp edrug ANC fn Eprol Emat Vcg KDE ktr CIR0 ken MTT GAM Emax_prol Emax_mat
  
