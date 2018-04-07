function bConverged = CheckConvergence(oGTMTTModel)
vLLhood    = get(oGTMTTModel,'LLhood');
dTolerance = get(oGTMTTModel,'Tolerance');
bVerb      = get(oGTMTTModel,'verb');
iter       = length(vLLhood);
LLini      = vLLhood(1);
LLhood     = vLLhood(iter);
bConverged = false;
if (iter < 2)
    return;    
elseif (LLhood-LLini)<(1 + dTolerance)*(vLLhood(iter-1)-LLini) && iter>50
    if bVerb 
        fprintf('MSG: The algorithm has reached the convergence.\n'); 
    end;
    set(oGTMTTModel,'status','trained');
    bConverged = true;
elseif ~isfinite(LLhood)
    if bVerb 
        fprintf('MSG: The training has been early stopped becouse the log-likelihood tends to infinite\n'); 
    end;
    set(oGTMTTModel,'status','failed');
    bConverged = true;
end        
