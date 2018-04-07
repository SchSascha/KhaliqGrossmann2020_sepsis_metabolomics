%
% Draws a grid of the latent space 
% and the mod projections of the data points
% to the latent space including lines for
% the temporal trajectory
% 
function membmap(noLatVar, mod, targets,Path,bFigureFlag)
%membmap(noLatVar, mod, targets,Path)

if(size(mod,2)~=2) % two dimensional mode plot only
    error('Not implemented for non 2-D structures')
end
    
if(~exist('bFigureFlag','var'))
    bFigureFlag = 0;
end
error(nargchk(2, 5, nargin, 'struct'));

if nargin>=3
    if sum(targets.*targets==targets)<=1
        Path = targets;
        targets = [];
    else
        Path =[];
    end
elseif nargin==2
    targets = [];
    Path=[];
end

grid = rctg(noLatVar);

if isempty(targets)
    greathit = gtm_hit2(grid,mod,'quiet');
else
    greathit = gtm_hit2(grid,mod,targets,'quiet');
end

if(~bFigureFlag)
    figure;
end
set(gcf,'Color',[0 1 1],'Name','Cluster Membership Map');
if ~isempty(targets)
    subplot(1,2,1);
end
r1=reshape(greathit.rati1,sqrt(noLatVar),sqrt(noLatVar));
gtm_plane('rect',r1,greathit.hitprop);%.*reshape(greathit.rati1,8,8));
set(gca,'FontSize',7);
%title('MEMBERSHIP MAP')
colormap('gray');

if ~isempty(targets)
    clas1nam = '1';
    clas0nam = '0';
    
    textit4=cell(2,1);textit4{1,1}='        MAP with "unit size" scale        ';textit4{2,1}=['WHITE units for class ' clas1nam ' BLACK units for class ' clas0nam];     
    	  title(textit4);      
    subplot(1,2,2);
    r1=reshape(greathit.rati1,sqrt(noLatVar),sqrt(noLatVar));
    gtm_plane('rect',r1,greathit.hitpfake);
    set(gca,'FontSize',7);
    textit4=cell(2,1);textit4{1,1}='       MAP without "unit size" scale      ';textit4{2,1}=['WHITE units for class ' clas1nam ' BLACK units for class ' clas0nam];
    title(textit4);
    colormap('gray');
    
end


if ~isempty(Path)
    mod1(:,2) = -0.5*(sqrt(noLatVar) - 1)*(mod(:,2)-1) + 1;
    mod1(:,1) = 0.5*(sqrt(noLatVar) - 1)*(mod(:,1)+1) + 1;
    hold on;
    arrow(mod1(Path(1:end-1),:), mod1(Path(2:end),:))
    plot(mod1(Path(1),1),mod1(Path(1),2),'s');
    plot(mod1(Path(end),1),mod1(Path(end),2),'o');
end



