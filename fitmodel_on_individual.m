function output = fitmodel_on_individual(exp, subject, model_type, optimize)

%Inputs:
%exp: expeirment to fit (1, 2 or 3)
%subject: subject ID (1-13 in Exp1; 1-11 in Exp2 and Exp3)
%model_type: 1: Max; 2: Difference 3: Entropy model
%optimzie: 0 to use patternsearhc, 1 to use BADS for optimization

%Outputs are saved in a structure:
%output.par: best-fit parameters
%output.nll: negative log likelihood
%output.AIC: AIC score
%output.x0: initial values of the free parameters
%output.exitflag: exitflag of the optimization function

if ~exist('optimize', 'var')
    optimize = 0;
end

%% Load data to fit
fprintf('fitting exp%i sub%i model%i ----- \n',exp, subject, model_type);
filename_data = sprintf('data/EXP%i/S%i_log.mat',exp,subject);
load(filename_data,'trl');

%Which configurations to fit. We fit all four configurations in each experiment. 
%Experiment 1 and 3: configuration 1~4; Experiment 2: configuration 5~8
switch exp 
    case {1,3}
        configpool = 1:4;
    case 2
        configpool = 5:8; 
end

%% Intiate experiment parameters
% x, y location of the center of each category. 8 configures are tested (1-4 in Experiment 1 and 3; 5-8 in Expeirment 2)
% Category locations of 1th to 8th configures correspond to 1th to 8th row in stiPar.ux and stiPar.uy respectively
stiPar.ux  = [-96,0,96; -128,0,128; -96,-59,96; -96,59,96; -64,0,64; -51,0,51; -64,-64,64; -64,64,64]; %horizontal location (relative to screen center)
stiPar.uy  = [0,0,0; 0,0,0; 0,0,0; 0,0,0; 37,-74,37; 30,-59,30; 37,0,37; 37,0,37]; %vertical location (relative to screen center)
stiPar.sig_s   = 64; %standard deviation of stimulus distribution (same for all three categories)
stiPar.nr      = 4;  %number of levels of confidence reports allowed
stiPar.ncat    = 3;  %number of categories

%% Start optimization

options    = psoptimset('Display','iter','MeshAccelerator','on'); %Options of optimization

switch model_type
    case 1
        % MAX model with sensory noise and inference noise
        %[b1, b2, b3, log10(sig_x), log10(alpha), lapse_rate]
        %sig_x and alpha can have infinite boundaries [-Inf Inf], here using 
        %[-10 10] ([10^-10 10^10] actually) as demo
        
        lb  = [0 0 0 -20 -20 0];
        ub  = [1 1 1 20 20 1];
        plb = [0.2 0.5 .75 1 .5  0];
        pub = [0.5 0.75 1  1.5 1.5 .2];
        objFunc = @(x)nll_map_is(x,trl,stiPar,configpool);
        
    case 2
        % DIFF model with sensory noise and inference noise
        %[b1, b2, b3, log10(sig_x), log10(alpha), lapse_rate]
        
        lb  = [0 0 0 -20 -20 0];
        ub  = [1 1 1 20 20 1];
        plb = [0   0.2 .5 1 .5  0];
        pub = [0.2 .5  1  1.5 1.5 .2];
        objFunc = @(x)nll_diff_is(x,trl,stiPar,configpool);
        
    case 3
        % Entropy model with sensory noise and inference noise
        %[b1, b2, b3, log10(sig_x), log10(alpha), lapse_rate]
        
        lb  = [0 0 0 -20 -20  0];
        ub  = [1 1 1  20 20  1];
        plb = [0 0.2 .7 1 .5  0];
        pub = [0.2 0.7 1 1.5 1.5 .2];
        objFunc = @(x)nll_en_is(x,trl,stiPar,configpool);
end

randseed = rand(1,length(plb));
x0 = plb + randseed.*(pub-plb);
x0(1:3) = sort(x0(1:3));

switch optimize
    case 0
        %Use patternsearch for optimization. An optimization function comes with MATLAB global optimization toolbox
        A = -[-1 1 0 0 0 0; 0 -1 1 0 0 0];
        b = [0,0];
        [estX,fval,exitflag] = patternsearch(objFunc,x0,A,b,[],[],lb,ub,[],options);
    case 1
        %Use BADS for optimization; Download BADS here: https://github.com/lacerbi/bads
        nonbcon = @(x) x(:,1)>x(:,2) | x(:,2)>x(:,3);
        [estX,fval,exitflag] = bads(objFunc,x0,lb,ub,plb,pub,nonbcon,options);
end

AIC = 2*(length(x0) + fval);

output.par = estX;
output.AIC = AIC;
output.nll = fval;
output.exitflag = exitflag;
output.x0 = x0;

%% Print some results in the command Window
fprintf('Done fitting Experiment %i Subject %i Model %i \n', exp, subject, model_type);
fprintf('AIC score %.1f \n', AIC);

