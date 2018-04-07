function plot_data_grid_2d(Data, Labels, Y, Y_Labels, nF,bFigure,iDim1,iDim2)

N=numel(Labels);
[K,D]=size(Y);

%plot data
if(bFigure)
    figure
else
    cla
end
hold on
for i=1:10:N
    switch Labels(i)
        case 1
            plot(Data(i,iDim1),Data(i,iDim2),'+b');
        case 2
            plot(Data(i,iDim1),Data(i,iDim2),'og');
        case 3
            plot(Data(i,iDim1),Data(i,iDim2),'*r');
    end
end

%plot gtm grid in data space
A=grid_adjacency(sqrt(K),sqrt(K));
gplot(A,Y(:,[1 2]));

for k=1:K
    switch Y_Labels(k)
        case 1
            plot(Y(k,iDim1), Y(k,iDim2), '+b','LineWidth',2,'MarkerSize',12);
        case 2
            plot(Y(k,iDim1), Y(k,iDim2), 'og','LineWidth',2,'MarkerSize',12);
        case 3
            plot(Y(k,iDim1), Y(k,iDim2), '*r','LineWidth',2,'MarkerSize',12);
    end
end

title(['nX=', num2str(sqrt(K)), ', nF=', num2str(nF), ', D=', num2str(D)]);

hold off
pause(0.1)