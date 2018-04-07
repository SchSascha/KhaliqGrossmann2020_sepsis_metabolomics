function [bResult,cField] = getfield(cClass, sString)
    cField = [];
    bResult = 1;
    switch sString
       case 'p'
               cField = cClass.p;
               return;               
       case 'tau'
               cField = cClass.tau;
               return;               
       case 'vLam'
               cField = cClass.vLam;
               return;               
       case 'bWithAdaptation'
               cField = cClass.bWithAdaptation;
               return;               
       otherwise
               oBaseObject = cClass.metric;
               [bResult, cField] = getfield(oBaseObject, sString);
               return;
    end
    bResult = 0;
return 