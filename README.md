# UnfoldingJPsiInJets
The code for 2D unfolding. It is adapted for the analysis "JPsi in jets", meaning the variables of the unfolding are z (J/Psi pt / jet pt) and jet pt

MC unfolding is in "mcUnfCode" directory, and the data one in "dataUnfCode". 

The first step is to create 4D transfer matrix : done with "prepareInputsMC.cc" or "prepareInputsData.cc" ;

The second step is to prepare the 4D transfer matrix for the unfolding (and also to prepare proper "measured" distribtuion for the MC closure test) : done with "createRooUnfoldResponseMC.cxx" or "createRooUnfoldResponseData.cxx" ;

The unfolding itself is done in "unfoldMCStep.cxx" or "unfoldDataStep.cxx". 
