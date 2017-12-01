function llr_g = patientGhmmGeneLLR(patientGhmm, x, varargin)
%PATIENTGHMMGENELLR returns a 1 by G vector of LLR according to each gene.
%   LLR_G = PATIENTGHMMGENELLR(PATIENTGHMM, X) returns a 1 by G vector of
%   log-likelihood-ratio according to each gene. PATIENTGHMM{m} is a cell
%   array of the two models, and X{q}(g,t) is the gene expression of
%   patient q, gene g and time t. LLR_G(1,g) is the log-likelihood-ratio
%   according to gene g. 

nPatient = length(x);
nGene = size(x{1}, 1);

llr_g = zeros(1, nGene);
loglik = zeros(2, nGene);
for q = 1 : nPatient
    for m = 1 : 2
        tr=patientGhmm{m}.tr; mu=patientGhmm{m}.mu; sigma=patientGhmm{m}.sigma;
        nState = size(patientGhmm{m}.mu,2);
%         estS = mghmmViterbi(x{q}, tr, mu, sigma, varargin{:});
%         for s = 1 : size(patientGhmm{m}.mu, 2) - 1
%             xs = x{q}(:, (estS == s + 1));
%             T = size(xs, 2);
%             loglik(m,:) = loglik(m,:) + sum(log(normpdf(xs,...
%                 repmat(mu(:,s+1),1,T), repmat(sigma(:,s+1),1,T))), 2)';
%         end
        [emilik, emilik3d] = mghmmEmiLik(x{q}, mu, sigma, varargin{:});
        T = size(emilik3d, 3);
        pStates = hmmDecodeLik(tr, emilik);
        % if seq likelihood is 0, skip this seq, following HMMER
        if size(pStates,1) == 0, continue; end %seq likelihood zero
        pStates = permute(repmat(pStates(1:nState, 1:T),[1,1,nGene]), [3,1,2]);
        %loglik(m,:) = loglik(m,:) + sum(sum(pStates.*emilik3d,2),3)';
        %loglik(m,:) = loglik(m,:) + sum(sum(pStates.*log(emilik3d),2),3)';
        loglik(m,:) = loglik(m,:) + sum(log(sum(pStates.*emilik3d,2)),3)';%GxSxT
    end
    llr_g = llr_g + loglik(2,:) - loglik(1,:);
end
