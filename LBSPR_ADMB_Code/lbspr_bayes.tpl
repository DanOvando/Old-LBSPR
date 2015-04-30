// TPL file for LBSPR estimation procedure
// GTG Model - Length Structured 
DATA_SECTION
  init_number logFMMean;
  init_number logFMSigma;
  init_number logSL50Mean;
  init_number logSL50Sigma;
  init_number logDeltaMean;
  init_number logDeltaSigma;
  
  init_number logMkMean;
  init_number logMkSigma;  
  init_number logLinfMean;
  init_number logLinfSigma;
  init_number CVLinflb; // lower bound on uniform distribution
  init_number CVLinfub; // upper bound on uniform distribution

  init_number logL50Mean;
  init_number logL50Sigma;
  init_number logMatDeltaMean;
  init_number logMatDeltaSigma;
  
  init_number logFecBMean;
  init_number logFecBSigma;
  
  init_number logMpowMean;
  init_number logMpowSigma;
  
  init_number NGTG;
  init_number MaxSD;
  init_number NLenMids;
  init_vector LenMids(1,NLenMids); 
  init_vector ObsLength(1,NLenMids); 
  init_vector LenBins(1,NLenMids+1); 
  

// ----------------------------------------------------------------------------
 PARAMETER_SECTION  
  init_number logFM;
  init_number logSL50;
  init_number logDelta;
  init_number logLinf;
  init_number logMk;
  init_number logL50;
  init_number logMatDelta;
  init_number logFecB;
  init_number logMpow;
  init_bounded_number CVLinf(CVLinflb,CVLinfub,1); 
  

  // Storage
  number Linf;
  number Mk;
  number L50; 
  number L95;
  number Mpow;
  number FecB;
  number Delta;
  number SDLinf; 
  number GTGLinf;
  number LinfdL;
  number MatDelta;

  vector DiffLinfs(1,NGTG);
  vector RecProbs(1,NGTG);
  
  matrix UnfishedMatrix(1,NGTG,1,NLenMids);
  matrix FishedMatrix(1,NGTG,1,NLenMids);
  vector PredLenComp(1,NLenMids);
  vector Vul(1,NLenMids+1);
  vector MkL(1,NLenMids+1);
  vector ZkL(1,NLenMids+1);
  
  vector Mat(1,NLenMids);
  vector Wght(1,NLenMids);
  vector Fec(1,NLenMids);
  
  vector EP0_gtg(1,NGTG);
  vector EPf_gtg(1,NGTG);
  number EP0;
  number EPf;
   
  
  vector PUnFished(1,NLenMids+1);
  vector PFished(1,NLenMids+1);
  vector NUnFished(1,NLenMids);
  vector NFished(1,NLenMids);
  
  sdreport_number FMpar;
  sdreport_number SPR;
  sdreport_number SL50; 
  sdreport_number SL95;
  
  objective_function_value obj_fun;
  
  !! cout << "Done Parameters" << endl;
  
// ---------------------------------------------------------------------------- 
  
PROCEDURE_SECTION
  int X;
  int Gtype, L;
  double pi = 3.14159265358979323844;
  
  // Maturity
  L50 = mfexp(logL50);
  MatDelta = mfexp(logMatDelta);
  L95 = L50 + MatDelta;

  Mat = 1.0/(1+mfexp(-log(19)*(LenMids-L50)/MatDelta));

  FecB = mfexp(logFecB);
  Fec = elem_prod(Mat, pow(LenMids, FecB));
  Fec = Fec/max(Fec);

  Linf = mfexp(logLinf); 
  SDLinf = CVLinf * Linf;
  LinfdL = ((Linf + MaxSD * SDLinf) - (Linf - MaxSD * SDLinf))/(NGTG-1);
  
  for (X=0;X<=NGTG;X++) {
   DiffLinfs(X+1) = (Linf - MaxSD * SDLinf) + X * LinfdL;  	 
 }
  

  
  RecProbs = 1/(sqrt(2*pi*SDLinf*SDLinf)) * mfexp(-(elem_prod((DiffLinfs-Linf),(DiffLinfs-Linf)))/(2*SDLinf*SDLinf));
  RecProbs = RecProbs/sum(RecProbs);
  
  FMpar = mfexp(logFM); // estimated F/M
  SL50 = mfexp(logSL50); 
  Delta = mfexp(logDelta); 
  Vul = 1.0/(1+mfexp(-log(19)*(LenBins-SL50)/Delta));
  
  Mk = mfexp(logMk);
  Mpow = mfexp(logMpow);
  MkL = Mk * pow(Linf/LenBins, Mpow);
  ZkL = MkL + (FMpar * Mk) * Vul;
  
  for (Gtype=1;Gtype<=NGTG;Gtype++) {

   PUnFished = 0; 
   PFished = 0; 
   NUnFished = 0; 
   NFished = 0 ;
   PUnFished(1) = RecProbs(Gtype);
   PFished(1) = RecProbs(Gtype);
   GTGLinf = DiffLinfs(Gtype); 
  
   for (L=2;L<=NLenMids+1;L++) {
    if (LenBins(L) < GTGLinf) {
      PUnFished(L) = PUnFished(L-1) * pow(((GTGLinf-LenBins(L))/(GTGLinf-LenBins(L-1))),MkL(L-1));
      PFished(L) = PFished(L-1) * pow(((GTGLinf-LenBins(L))/(GTGLinf-LenBins(L-1))),ZkL(L-1));
    }
    if (LenBins(L) >= GTGLinf) {
      PUnFished(L) = 0;
      PFished(L) = 0; 
    }
  }
 
  for (L=1;L<=NLenMids;L++) {
	NUnFished(L)  = (PUnFished(L) - PUnFished(L+1))/MkL(L);
	NFished(L) =  (PFished(L) - PFished(L+1))/ZkL(L);	
  }

 UnfishedMatrix(Gtype) =  NUnFished;
 FishedMatrix(Gtype) = NFished;  
 
 EP0_gtg(Gtype) = sum(elem_prod(NUnFished, Fec));
 EPf_gtg(Gtype) = sum(elem_prod(NFished, Fec));
 }
 
  EP0 = sum(EP0_gtg);
  EPf = sum(EPf_gtg);
  SPR =  EPf/EP0;
  
  PredLenComp = colsum(FishedMatrix);
  PredLenComp =  elem_prod(PredLenComp,  1.0/(1+mfexp(-log(19)*(LenMids-SL50)/Delta)));
  PredLenComp = PredLenComp/sum(PredLenComp);
  SL95 = SL50 + Delta;
 
  // obj_fun = -sum(elem_prod(ObsLength, log(PredLenComp))); // NLL 
  obj_fun = -sum(elem_prod(ObsLength, log(PredLenComp+0.00001))); // NLL 

  
  obj_fun+= 1/(2*square(logFMSigma)) * square(logFMMean - logFM);	// add NLL(FM)
  obj_fun+= 1/(2*square(logSL50Sigma)) * square(logSL50Mean - logSL50);	// add NLL(SL50)
  obj_fun+= 1/(2*square(logDeltaSigma)) * square(logDeltaMean - logDelta);	// add NLL(Delta)
  obj_fun+= 1/(2*square(logMkSigma)) * square(logMkMean - logMk);	// add NLL(MK)
  obj_fun+= 1/(2*square(logLinfSigma)) * square(logLinfMean - logLinf);   // add NLL(Linf)
  obj_fun+= 1/(2*square(logL50Sigma)) * square(logL50Mean - logL50);	// add NLL(L50)
  obj_fun+= 1/(2*square(logMatDeltaSigma)) * square(logMatDeltaMean - logMatDelta);	// add NLL(MatDelta)
  obj_fun+= 1/(2*square(logFecBSigma)) * square(logFecBMean - logFecB); // add NLL(FecB)
  obj_fun+= 1/(2*square(logMpowSigma)) * square(logMpowMean - logMpow);	// add NLL(Mpow)
 

 
// ---------------------------------------------------------------------------- 

// ---------------------------------------------------------------------------- 
REPORT_SECTION

 
 report << NLenMids << endl; // number of length classes
 report << PredLenComp << endl; // model fit
 report << ObsLength << endl; // original length data
 report << LenMids << endl; // length bins
 // report << SPR << endl; // estimate of SPR
 // report << FMpar << endl; // estimate of F/M 
 // report << SL50 << endl; // estimate of SL50
 // report << SL95  << endl; // estimate of SL95
 report << obj_fun << endl; // likelihood obj_fun
 report << objective_function_value::pobjfun->gmax; // maximum gradient value
 


 
// ---------------------------------------------------------------------------- 
TOP_OF_MAIN_SECTION
 arrmblsize = 500000000;
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
 gradient_structure::set_MAX_NVAR_OFFSET(5000);
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
  
