algoParam.isInvivo      = true;     % true for using initial guesses for in vivo imaging
algoParam.isParallel    = false;     % true: using parfor parallel processing; false: no parfor
algoParam.DEBUG         = false;    % true: debug mode to display some info
algoParam.isNormData    = false;     % true: normalise the data by a global constant so that absolute tolerance will be used for fitting, here we use false because we did normalisation outside the function, see below
% fitting option
algoParam.maxIter       = 200;      % maximum number of iterations
algoParam.fcnTol        = 1e-4;     % fitting tolerance, this valuse is for normalised data
algoParam.stepTol       = 1e-4;     % step tolerance, this valuse is for normalised data
% residual option
algoParam.numMagn       = 1;        % 0: complex fitting
algoParam.isWeighted    = true;     % true: using magnitude signal weighting
algoParam.weightMethod  = 'quadratic_1stEcho';  % '1stEcho': weighted by 1st echo magnitude; 'quadratic_1stEcho': weighted by squared 1st echo magnitude
% T1 model
algoParam.isExchange    = 1;        % BM model
algoParam.isEPG         = 1;        % Using EPG-X for signal simulation
algoParam.npulse        = 50;       % number of pulses to reach steady-state
if ~isfield(input,'MRvendor')
    algoParam.rfphase       = 50;       % RF phase for EPG, degree
else
    if strcmp(input.MRvendor,'siemens')||strcmp(input.MRvendor,'Siemens')||strcmp(input.MRvendor,'SIEMENS')
        algoParam.rfphase       = 50;       % RF phase for EPG, degree
    else
        algoParam.rfphase       = 150;       % phase increment for Philips
    end
end

algoParam.isT1mw        = false;    % true: fitting myelin water T1, false: use fixed value
algoParam.T1mw          = 234e-3;   % define fixed myelin water T1 value
% No DIMWI
algoParam.DIMWI.isVic       = false;    % false: no extra DWI info for DIMWI
algoParam.DIMWI.isR2sEW     = false;    % false: no extra DWI info for DIMWI
algoParam.DIMWI.isFreqMW    = false;    % false: no extra DWI info for DIMWI
algoParam.DIMWI.isFreqIW    = false;    % false: no extra DWI info for DIMWI
% initial guess
%algoParam.advancedStarting = 'robust';   % initial guesses for multi-comp S0
algoParam.advancedStarting = 'default';  % initial guesses for multi-comp S0