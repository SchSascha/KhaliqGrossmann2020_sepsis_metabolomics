%% call me like this 
% oModel = vModels{1}.oModelR; plot_data_grid_2d(X, L, get(oModel,'y'), ones(size(get(oModel,'y'),1),1), '1',false,6,7)
function plot_data_grid_3d(Data,Labels,Y, Y_Labels, nF,bFigure)

N=numel(Labels);
[K,D]=size(Y);

%plot data
if(bFigure)
    figure
else
    cla
end
hold on
for i=1:N
    switch Labels(i)
        case 1
            plot3(Data(i,1),Data(i,2),Data(i,3),'or');
        case 2
            plot3(Data(i,1),Data(i,2),Data(i,3),'xb');
        case 3
            plot3(Data(i,1),Data(i,2),Data(i,3),'sg');
    end
end

%plot gtm grid in data space
A=grid_adjacency(sqrt(K),sqrt(K));
gplot(A,Y(:,1:3));

for k=1:K
    switch Y_Labels(k)
        case 1
            plot3(Y(k,1), Y(k,2), Y(k,3), 'or','LineWidth',2,'MarkerSize',12);
        case 2
            plot3(Y(k,1), Y(k,2), Y(k,3), 'xb','LineWidth',2,'MarkerSize',12);
        case 3
            plot3(Y(k,1), Y(k,2), Y(k,3), 'sg','LineWidth',2,'MarkerSize',12);
    end
end

title(['nX=', num2str(sqrt(K)), ', nF=', num2str(nF), ', D=', num2str(D)]);
hold off;
