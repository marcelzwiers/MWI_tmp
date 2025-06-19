kappa_mw                = 0.36; % Jung, NI., myelin water density
kappa_iew               = 0.86; % Jung, NI., intra-/extra-axonal water density
fixed_params.B0     	= sepia_header{end}.B0;    % field strength, in tesla
fixed_params.rho_mw    	= kappa_mw/kappa_iew; % relative myelin water density
fixed_params.E      	= 0.02; % exchange effect in signal phase, in ppm
fixed_params.x_i      	= -0.1; % myelin isotropic susceptibility, in ppm
fixed_params.x_a      	= -0.1; % myelin anisotropic susceptibility, in ppm
fixed_params.B0dir      = sepia_header{end}.B0_dir;
fixed_params.t1_mw      = 234e-3;


fitting = [];
fitting.Nepoch              = 4000;
fitting.initialLearnRate    = 0.001;    %    start from 0.001 % 0.01 I could not get MWF betweeen 0 and 6...
fitting.decayRate           = 0;
fitting.convergenceValue    = -inf% 1e-8;     %   -inf in this case it keeps running
fitting.tol                 = 1e-8;
fitting.display             = false;
fitting.lossFunction        = 'l1';
fitting.start               = 'prior';   
fitting.patience            = 20;       %   it waits n iterations after getting to a local minimum


fitting.DIMWI.isFitIWF      = 1;
fitting.DIMWI.isFitFreqMW   = 1;
fitting.DIMWI.isFitFreqIW   = 1;
fitting.DIMWI.isFitR2sEW    = 1;
fitting.isFitExchange       = 1;
fitting.isEPG               = 0;

