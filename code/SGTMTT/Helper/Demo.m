% define three sample classes in 1-D and map them to 4-D
stream = RandStream.getGlobalStream;
reset(stream);
N=50; % samples per class 
T=30; % time points        
K=9 ; % number of states
D=zeros(N*3,T);S=D;R=D;
for(k=1:N);[vSignal,vRand,vState,mP] = FiniteAutomateShort(T);         D(k,:)     = vSignal;   S(k,:)     = vState; R(k,:)     = vRand; end
for(k=1:N);[vSignal,vRand,vState,mP] = FiniteAutomateLoopForward(T,1); D(N+k,:)   = vSignal;   S(N+k,:)   = vState; R(N+k,:)   = vRand; end
for(k=1:N);[vSignal,vRand,vState,mP] = FiniteAutomateLoopForward(T);   D(2*N+k,:) = vSignal;   S(2*N+k,:) = vState; R(2*N+k,:) = vRand; end
L=[ones(N,1);2*ones(N,1);3*ones(N,1)];
X=[]
for(k=1:6)
    X(:,:,k) = 0.001*rand(N*3,T);
end
X(:,:,3) = D+(0.01*rand(N*3,T));

% learn the supervised gtmtt
[vModels] = GTMTT_Supervised(X,L,K,4);

% plot the learned metric parameter
pause
E=get(vModels{1}.oModel1,'ometric')
bar(E.Par)

% plot the map and the path
for(c=1:3)
    Resp = get(vModels{c}.oModel1,'r');
    mod  = pmd(Resp{1});
    iTDim = size(X,2);
    %... now including paths
    membmap(K,mod,[1:iTDim]);
end
