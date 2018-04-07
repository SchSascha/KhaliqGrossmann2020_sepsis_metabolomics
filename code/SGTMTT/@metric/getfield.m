function [bResult,cField] = getfield(cClass, sString)
    cField = [];
    bResult = 1;
    switch sString
        case 'descriptor'
               cField = cClass.descriptor;
               return;               
        case 'type'
               cField = cClass.type;
               return;               
    end
    bResult = 0;
return 