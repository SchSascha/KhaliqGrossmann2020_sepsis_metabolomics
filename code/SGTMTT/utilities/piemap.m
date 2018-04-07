function piemap(noLatVar, clusm, targets)
% Display a cluster membership map describing each node by a pie chart.
% This visualization is useful when it is desired incorporate class
% information "post-training"
%
%   PIEMAP(K,CLUSM,TARGETS)
%
% K is the number of latent points (it must be a squared integer), CLUSM is
% a vector which each of its component represent the cluster "assigned" by
% the mode projection procedure (It can be calculated by PMD function).
% TARGETS is a vector which represents the corresponding class for each
% data point. TARGETS accepts only integer numbers, starting with 1 (for
% class 1) and so on.
%
% Author: Iván Olier,       2009

error(nargchk(2, 3, nargin));
maxtarg = max(targets);
sqrLV = sqrt(noLatVar);

for k=1:noLatVar
    ctotxcluster(k) = sum(clusm==k);       
    targtxclus = targets(clusm==k);
    for t1=1:maxtarg
        cntxclust{k}(t1) = sum(targtxclus==t1);
    end
end

ctotxcluster = ctotxcluster/max(ctotxcluster);


figure('Name','Cluster Membership Map');
k=0;
for i=1:sqrLV
    for j=1:sqrLV
        k=k+1;
        x=cntxclust{k};
        %x = x(x~=0);
        xsum = sum(x);
        if xsum~=0
            if xsum > 1+sqrt(eps), x = x/xsum; end
            subplot(1,2,1),
            piec(x,i,j,ctotxcluster(k));
            title('MAP with "unit size" scale');
            subplot(1,2,2)
            piec(x,i,j,1);
           title('MAP without "unit size" scale');
        end
    end
end

for map=1:2
    subplot(1,2,map);
    axis('ij');
    axis('equal');
    axis([.5 sqrLV+.5 .5 sqrLV+.5]);
    set(gca,'XaxisLocation','Top');
    grid on;
    set(gca,'Box','On')
end

function piec(x,px,py,spi)
theta0 = -pi/2;
maxpts = 100;

for i=1:length(x)
  n = max(1,ceil(maxpts*x(i)));
  r = [0;ones(n+1,1);0];
  theta = theta0 + [0;x(i)*(0:n)'/n;0]*2*pi;
  [xx,yy] = pol2cart(theta,0.5*spi*r);
  theta0 = max(theta);
  patch('XData',xx+px,'YData',yy+py,'CData',i*ones(size(xx)),...
     'FaceColor','Flat' );        
end
