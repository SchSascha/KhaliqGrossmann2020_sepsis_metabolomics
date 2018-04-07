function [cClassReturn] = setfield(cClass, sString, oObject)
    bResult = 1;
    cClassReturn = cClass;
    [bResult] = isfield(cClass,sString);
    if(bResult)
       cClass.(sString) = oObject;
       cClassReturn = cClass;
    end
return 