function [nllk] = nll_en_is(x,trl,stiPar,config_pool,condIdx)
%--Entropy model with both inference noise and sensory noise--
% This script uses Monte-Carlo simualtions to estimate response distribution
% and computes negative log-likelihood of the model.
% trl: structure containing data and config label
% required field: trl.estD, trl.estC, trl.config, trl.target_x
% stiPar: the structure containing stimulus locations

if ~exist('condIdx','var') || isempty(condIdx)
    condIdx = ':';
end

simntrl = 10000; %sample size of MonteCarlo simulation per trial

ux = stiPar.ux;  
uy = stiPar.uy;  
sig_s = stiPar.sig_s;   
nr = stiPar.nr;      
ncat = stiPar.ncat;

%%
trl.estD = trl.estD(condIdx);
trl.estC = trl.estC(condIdx);
trl.config = trl.config(condIdx);
trl.target_x = trl.target_x(condIdx);
trl.target_y = trl.target_y(condIdx);

%%
k     = sort([-Inf x(1:3) Inf]); %boundaries for 4-point confidence report
sig_x = 10^x(5); %sensory noise
alpha = 10^x(4); %Diricelet inference noise
lapse = x(6);

%%
pMat = cell(max(config_pool),1);
for config = config_pool
    
    temp_s   = trl.target_x(trl.config==config)';
    ntrl     = length(temp_s);
    ux_vec    = ux(config,:);
    pMat{config} = nan(ncat, nr, ntrl);
    
    sim_xMat = repmat(temp_s,simntrl,1) + randn(simntrl,ntrl)*sig_x; %measurement matrix: simntrl X ntrial
    sim_xMat = repmat(sim_xMat,1,1,ncat); %measurement matrix: simntrl X ntrial X ncat
    cMatx     = kron(ux_vec',ones(simntrl*ntrl,1));
    cMatx     = reshape(cMatx,simntrl,ntrl,ncat);
    %sim_lh   = 1./sqrt(2*pi*sig_s^2) .* exp(-(sim_xMat - cMatx).^2 ./ (2*sig_s^2)); %likelihood: category X ntrial
    
    %Run this part for 2 dimensions
    temp_sy   = trl.target_y(trl.config==config)';
    uy_vec    = uy(config,:);
    sim_yMat = repmat(temp_sy,simntrl,1) + randn(simntrl,ntrl)*sig_x; %measurement matrix: simntrl X ntrial
    sim_yMat = repmat(sim_yMat,1,1,ncat); %measurement matrix: simntrl X ntrial X ncat
    cMaty    = kron(uy_vec',ones(simntrl*ntrl,1));
    cMaty    = reshape(cMaty,simntrl,ntrl,ncat);
    
    sig_p = sqrt(sig_s^2 + sig_x^2);
    sim_lh   = 1./(2*pi*sig_p*sig_p) .* exp(-((sim_xMat - cMatx).^2 +(sim_yMat - cMaty).^2) ./ (2*sig_p^2)); %likelihood: category X ntrial
    
    %Inject Dirichlet decision noise
    sim_post = bsxfun(@rdivide,sim_lh,sum(sim_lh,3));
    sim_post = gamrnd(sim_post*alpha,1);
    sim_post = bsxfun(@rdivide,sim_post,sum(sim_post,3));
    [~,sim_d] = max(sim_post,[],3); %category decision
    sim_post(sim_post==0) = eps; %deal with numerical issue
    sim_c_en   = 1-(-sum(sim_post.*log2(sim_post),3) / log2(ncat)); %ENTROPY confidence

    %Get probability map given a set of parameter value
    [~,sim_r]  = histc(sim_c_en,k,1); %categorical confidence report based on k
    sim_r = max(min(sim_r,4),1); %prevent some numerical issue (if ever happens)
    
    for tr = 1:ntrl
        for d = 1:ncat
            for r = 1:nr
                temp = sum(sim_d(:,tr)==d & sim_r(:,tr)==r)/simntrl;
                pMat{config}(d,r,tr) = (1-lapse)*temp + lapse*(1/(ncat*nr));
            end
        end
    end
    pMat{config} = max(pMat{config},eps);
end

%%
nllk = 0;
for config = config_pool
    [ncat,nr,ntrl] = size(pMat{config});
    temp_D = trl.estD(trl.config==config);
    temp_C = trl.estC(trl.config==config);
    ind = sub2ind([ncat,nr,ntrl],temp_D,temp_C,(1:ntrl)');
    llk = sum(log(pMat{config}(ind)));
    nllk = nllk-llk;
end

