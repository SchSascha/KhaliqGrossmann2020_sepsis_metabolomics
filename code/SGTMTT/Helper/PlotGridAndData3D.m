function PlotGridAndData3D(oModel,iDim1,iDim2,iDim3,iTimePoint,T,D,bFlag)
    sLineColor = 'b';
    if(~bFlag)
        hold;
        sLineColor = 'r';
    end
    L=sqrt(get(oModel,'K')); % length of the map (square)
    Y=get(oModel,'y');
    X=get(oModel,'x');
    MY=Y(:,[iDim1,iDim2,iDim3]);
    if(exist('som_connection','file'))
        [r,c] = find(full(som_connection({'rect',[L L],'sheet'})));
        CoordGrid = sortrows([r c]);
    else % only for 3 x 3
        CoordGrid = [1 2;1 4; 2 3; 2 5; 3 6; 4 5; 4 7; 5 6; 5 8; 6 9; 7 8; 8 9];
    end
    line([MY(CoordGrid(1:end,1),1),MY(CoordGrid(1:end,2),1)]',[ MY(CoordGrid(1:end,1),2) MY(CoordGrid(1:end,2),2)]',[ MY(CoordGrid(1:end,1),3) MY(CoordGrid(1:end,2),3)]','Color',sLineColor)
    text(MY(:,1),MY(:,2),MY(:,3),cellstr(num2str([1:size(Y,1)]')))
    for(k=1:length(X))
        vCurrent = reshape(permute(X{k},[2,1,3]),T,D);
        plot3(vCurrent(iTimePoint,iDim1),vCurrent(iTimePoint,iDim2),vCurrent(iTimePoint,iDim3),'x','Color',sLineColor);
    end