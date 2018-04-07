function mConfusion = ConfusionMatrix(vClusters,vExpectedCluster,vUniqueClasses)
    vUniqueClasses       = [vUniqueClasses; numel(vUniqueClasses)+1]; % add reject class
    iNumberOfClasses     = length(vUniqueClasses);    
    mTemp                = bsxfun(@eq,vClusters,vUniqueClasses');
    mConfusion           = zeros(iNumberOfClasses-1,iNumberOfClasses);
    for(k=1:iNumberOfClasses-1) % rows - expected
        mConfusion(k,:) = sum(mTemp(find(vExpectedCluster==vUniqueClasses(k)),:)',2); % find entries of current class and how they are mapped to the observed labels
    end
end