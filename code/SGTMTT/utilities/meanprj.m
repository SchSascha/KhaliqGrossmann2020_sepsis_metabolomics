function u_mean = meanprj(resp,targets)
% Calculate the mean projection

error(nargchk(1, 2, nargin, 'struct'))

if nargout>1
    error('the number of output args is wrong');
end



[K N] = size(resp);

U = rctg(K);
u_mean = resp'*U;

if nargout == 0
    if nargin==1 targets=ones(N,1); end;
    maxtarg = max(targets);
    markedvector = {'ko','b+','r<','c>','g^','mv','ys','kd'};
    for targ=1:maxtarg
        plot(u_mean(targ==targets,1),u_mean(targ==targets,2),markedvector{targ});
        hold on;
    end
   % plot(U(:,1),U(:,2),'ro');
    axis([-1.1 1.1 -1.1 1.1]);
    hold off;
    clear u_mean;
   % grid on;
end
