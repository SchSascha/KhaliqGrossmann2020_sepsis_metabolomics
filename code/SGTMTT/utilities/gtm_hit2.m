function hit=gtm_hit2(X,mod,datarg,q)
% Maps of fuzzy class-membership, or hit maps
% for the GTM.
% X - latent nodes
% mod - modes
% [datarg: optional argument] - data structure
% [q : optional argument]- quiet mode (no displays if string 'quiet':quiet=1)
% ======
% Returns the structure hit, containing the fields:[count],countot,empties,noempts,clusmemb.
% ======
% Alfredo Vellido (September 1999)
% Revised and modified (June 2001), to accomodate the possibility of not having targets,
% ...and July 2001, to accomodate GUI
%%%%%%%%%

% Possible nargin errors...
if (nargin>=2 | nargin<=4)
   if (nargin == 2)
      quiet=0;
      targmode=0;
   elseif (nargin==3)
      if (isstr(datarg))
         if (strcmp(datarg, 'quiet'))
            quiet = 1;
            targmode=0;
         else
            errordlg('Input argument q should be the STRING quiet','ERROR in gtm_hit2');
         end
      else
         quiet=0;
         targmode = 1;
         errordlg('Input argument q should be a STRING','ERROR in gtm_hit2');
      end
   elseif (nargin == 4)
      targmode = 1;
      if (strcmp(q, 'quiet'))
         quiet = 1;
      else
         errordlg('Input argument q should be the STRING quiet','ERROR in gtm_hit2');
      end
   else
      errordlg('The number of input arguments should be between 2 and 4','ERROR in gtm_hit2');
   end
end

% Initialisation of variables...

noLatVar=size(X,1);
if targmode, count=zeros(size(X,1),2); end  % count of observations belonging to each class, assigned to each node 
countot=zeros(size(X,1),1);  % total number of observations assigned to each node
clusmemb=zeros(1,length(mod));  % node membership of each observation

for i=1:length(mod)
   for j=1:noLatVar
      if mod(i,:)==X(j,:)
         clusmemb(i)=j;
         if ~targmode
            countot(j)=countot(j)+1;
         else
            if datarg(i)==1 
               count(j,1)=count(j,1)+1;
            else
               count(j,2)=count(j,2)+1;
      		end
    	 end
   	  end
  	end
end

if targmode, countot=(sum(count'))'; end
empties=find(countot==0);  % Which are the nodes with no observation assigned?
noempts=find(countot~=0);  % ...and which are not empty?

% Proportion of observations falling on each unit
hitprop=zeros(sqrt(noLatVar),sqrt(noLatVar));
hitpfake=zeros(sqrt(noLatVar),sqrt(noLatVar));

hitprop(noempts)=countot(noempts)./max(countot);
hitprop(empties)=0;
hitpfake(noempts)=1;

rati1=zeros(1,length(X)); % Initialise the class proportion in each unit (leave it as zeros if we 
                          % do not have class membership information)
if targmode
   class1pr=sum(datarg)/length(datarg); % Prior of class 1
   rati1(noempts)=(count(noempts,1)./class1pr)./(count(noempts,1)./class1pr+count(noempts,2)./(1-class1pr));
   rati1(empties)=1; rati0=1-rati1;
   if ~quiet
      figure  % It only makes sense to plot the evenly-sized units map when we have class information
      r1=reshape(rati1,sqrt(noLatVar),sqrt(noLatVar));
      gtm_plane('rect',r1,hitpfake);
      colormap('gray') % In presence of a GUI, this could be considered an option
   end
end

% Class-membership ratio map with node size scaled according to number of assigned
% observations. If no targets available, homogeneously shaded display 
if ~quiet
   %figure
   r1=reshape(rati1,sqrt(noLatVar),sqrt(noLatVar));
   gtm_plane('rect',r1,hitprop);
   colormap('gray')
end

% Organize output
if targmode, hit.count=count; hit.hitpfake=hitpfake; end
hit.countot=countot;
hit.empties=empties;
hit.noempts=noempts;
hit.clusmemb=clusmemb;
hit.hitprop=hitprop;
hit.rati1=rati1;
