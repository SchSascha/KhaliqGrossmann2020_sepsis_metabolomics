%%
% Plot data upto 5 classes
% vUnique =unique(SDL.labels);Marker = [{'r'},{'g'},{'b'},{'y'},{'k'}]; hold; for(k=1:numel(vUnique)) I=find(SDL.labels == vUnique(k)); plot(SDL.data(I,:,2)','b-s','MarkerFaceColor',Marker{k});end
%
%%
function mAccs = TestModelOnFrameShifted(SDL,vFrame,varargin)
iRepeats = 1;
if(nargin==3)
    iRepeats = varargin{1};
end
mAccs = [];
for(i=1:iRepeats)
    sFilename = sprintf('frame_model_nose_800_94_repeat_%d.mat',i);
    T    = size(SDL.data,2); % equal times assumed
    N    = size(SDL.data,1);
    perm = randperm(N);
    TrainIndices = perm(1:ceil(N*0.75));
    TestIndices  = perm(find(~ismember(perm,TrainIndices)));
    SDLTrain = SDL;
    SDLTrain.data   = SDLTrain.data(TrainIndices,:,:);
    SDLTrain.labels = SDLTrain.labels(TrainIndices,:);
    SDLTest.data    = SDL.data(TestIndices,:,:);
    SDLTest.labels  = SDL.labels(TestIndices,:);
    [oModelResults.Models,vUniqueLabels] = GTMTT_Supervised(SDLTrain.data,SDLTrain.labels,9,4);
    Model = oModelResults.Models;
    save(sFilename,'Model','vUniqueLabels','perm');
    vAcc = zeros(T,1); % T shifts
    for(k=1:T) % this is to shift the frame over the time range   
        CurrentFrame = k+vFrame;
        vValidRange  = find(CurrentFrame>0 & CurrentFrame <= T);
        if(~isempty(vValidRange))
            CurrentFrame = CurrentFrame(vValidRange);            
            vAcc(k) = TestModelOnFrame(SDLTrain,SDLTest,oModelResults,CurrentFrame,vUniqueLabels);
        else
            vAcc(k) = 0;
        end
    end
    mAccs = [mAccs;vAcc];
    save(sFilename,'Model','vUniqueLabels','vAcc','vFrame','SDL','perm');
end