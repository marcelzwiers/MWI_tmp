
function func_MCR_AfterCoregistration_gpu(input,output)
% % if I am to make this into a function I need the following input:
% input.acq_str{countflip} %input acq_str of MGRE data
%
% input.derivative_SEPIA_dir % directory where the following datasets are stored
% %     magn_fn         = [gre_basename '_part-mag_MEGRE_space-withinGRE.nii.gz'];
% %     phase_fn        = [gre_basename '_MEGRE_space-withinGRE_part-phase_unwrapped.nii.gz'];
% %     totalField_fn   = [gre_basename '_MEGRE_space-withinGRE_fieldmap.nii.gz'];
% %     % mask_fn         = [gre_basename '_MEGRE_space-withinGRE_mask_refine.nii.gz'];
% % %     mask_fn         = [gre_basename '_MEGRE_space-withinGRE_mask_localfield.nii.gz'];
% input.derivative_FSL_dir % directory where B1 map saved in the MEGRE space is stored
% input.B1scaleFactor scaling factor for B1 map
% input.subj_label
% input.run_label
% input.MRvendor "Siemens"(deafault) or "Philips" selects the rf spoiling phase
% input.CorrectionFactorVFA  scale factor for VFA acquisition (should be some length as number of flip angles) 
% input.Configfile % config file for MCR-MWI fitting parameters

if ~isfield(output,'MPPCAdenoise')
    output.MPPCAdenoise = 0;
end

% example script to process MCR-MWI on a local computer with multi-threads
% The MCR-MWI processing scripts inside the subject folders were setted up
% slightly differently because of the parallelisation on HPC vs local PC

%% load GRE data
% for kp = protocol

% prot_SEPIA_dir  = fullfile(derivative_SEPIA_dir, kp{1});
% load GRE data
img             = [];
unwrappedPhase  = [];
mask            = [];
totalField      = [];
sepia_header    = [];
nFA             = (length(input.acq_str));
fa              = zeros(1,nFA);
for countflip= 1:length(input.acq_str)

    seq_SEPIA_dir = fullfile(input.derivative_SEPIA_dir,input.acq_str{countflip});

    % general GRE basename
    gre_basename    = [input.subj_label '_' input.acq_str{countflip} '_' input.run_label];

    % magnitude nifti image filename
    magn_fn         = [gre_basename '_part-mag_MEGRE_space-withinGRE.nii.gz'];
    phase_fn        = [gre_basename '_MEGRE_space-withinGRE_part-phase_unwrapped.nii.gz'];
    totalField_fn   = [gre_basename '_MEGRE_space-withinGRE_fieldmap.nii.gz'];
    % mask_fn         = [gre_basename '_MEGRE_space-withinGRE_mask_refine.nii.gz'];
    mask_fn         = [gre_basename '_MEGRE_space-withinGRE_mask_localfield.nii.gz'];

    sepia_header_fn = [gre_basename '_header.mat'];

    nii                     = load_untouch_nii(fullfile(input.derivative_SEPIA_dir, magn_fn));
    img                     = cat(5,img,nii.img);
    nii             = load_untouch_nii(fullfile(seq_SEPIA_dir, phase_fn));
    unwrappedPhase = cat(5,unwrappedPhase,nii.img);
    sepia_header{countflip}	= load(fullfile(input.derivative_SEPIA_dir, sepia_header_fn));
    mask            = cat(5,mask, load_nii_img_only(fullfile(seq_SEPIA_dir, mask_fn)));

    totalField = cat(4,totalField, load_nii_img_only(fullfile(seq_SEPIA_dir, totalField_fn)));

    fa(countflip)   	= sepia_header{countflip}.FA;
    tr              = sepia_header{countflip}.TR; % note that here there is an assumption that all protocols have the same TR
    te              = sepia_header{countflip}.TE;

end
dims =size(img);
% B1 info
% true_flip_angle_fn      = [subj_label '_' sess_label '_acq-famp_run-1_TB1TFL_space-withinGRE-' kp{1} '.nii.gz'];
true_flip_angle_json	= [input.subj_label '_acq-famp_run-1_TB1TFL.json'];
true_flip_angle_fn      = [input.subj_label '_acq-famp_run-1_TB1TFLProtocolSpace.nii.gz'];
true_flip_angle         = load_nii_img_only( fullfile(input.derivative_FSL_dir, true_flip_angle_fn));

% b1_header               = jsondecode( fileread( fullfile( converted_b1_dir,true_flip_angle_json)));
%
% b1                      = true_flip_angle / 10 / b1_header.FlipAngle;
b1                      = true_flip_angle / input.B1scaleFactor;

figure
Orthoview2(sum(mask,5),[],[],'tight' )
mask= prod(mask,5);
% mask_fn         = [subj_label '_' sess_label '_acq-' kp{1} '_' run_label '_brain_mask.nii.gz'];

clear true_flip_angle

%% obtain the initial estimation of the initial B1 phase
img = img .* exp(1i*unwrappedPhase);

mask_nonnan = min(min( ~isnan(img), [], 5), [],4);

mask = and(mask,mask_nonnan);
mask = prod(mask,4);

%%


if ~isfield(input,'Configfile')

    input.Configfile = [];
    [fixed_params] =get_default_GPUfixParam(input,sepia_header)
    [fitting] =     get_default_GPUfitParam(input)
else

    % load the config file
    if isfile(input.Configfile)
        run(input.Configfile);
    else
        error('Config file not found');
    end
end;





%% denoising part
if ~exist('MPPCAdenoise')
    MPPCAdenoise=0;
end

if output.MPPCAdenoise == 1

    [denoised,~,~] = denoise(reshape(img,[dims(1:3) prod(dims(4:5))]),[5 5 5],mask);
    img = reshape(denoised,dims);
    clear denoised

    PreProcessing ='MPPCAdenoising';
else
    if output.MPPCAdenoise == 2

        [denoised,~,~] = denoise(reshape(img,[dims(1:3) prod(dims(4:5))]),[3 3 3],mask);
        img = reshape(denoised,dims);
        clear denoised

        PreProcessing ='MPPCAdenoising3';
    else
        PreProcessing ='no_preprocessing';
    end;
end;
%%

pini = squeeze(unwrappedPhase(:,:,:,1,:)) - 2*pi*totalField .* sepia_header{end}.TE(1);

Debug = 0;
if Debug==1
% figure
for k=1:nFA
    figure(100)
    tiledlayout(ceil(sqrt(nFA)), ceil(sqrt(nFA)))
    nexttile(1)
    nexttile
    Orthoview2(pini(:,:,:,k),[],[],'tight')
    title(['Initial phase map for FA = ', num2str(fa(k))])
end;
end;


pini = polyfit3D_NthOrder( double(mean(pini(:,:,:,1:(end-1)), 4)), prod(mask,4), 6);

clear unwrappedPhase


% % data normalisation
mask_tmp = mask>0;

if isfield(input,'CorrectionFactorVFA')
    nFA = length(input.CorrectionFactorVFA);
    for countflip= 1:nFA
        img(:,:,:,:,countflip) = img(:,:,:,:,countflip) * (input.CorrectionFactorVFA(countflip));
    end
end

[scaleFactor, ~] = mwi_image_normalisation(img, mask_tmp);
clear mask_tmp


imgParam.output_dir = fullfile(output.derivative_MWI_dir, 'MCR', PreProcessing, ['using_', num2str(nFA),'_flipangle'],['quadraticW']);


%%
objGPU                     = gpuMCRMWI(te,tr,fa,fixed_params);

% slices = round(3/4*dims(3));
slices = 1:dims(3);

% fitting.display  = 1;

extraData = [];
extraData.freqBKG   = single(squeeze(totalField(:,:,slices,:,:)) / (gpuGREMWI.gyro*fixed_params.B0)); % in ppm
extraData.pini      = single(pini(:,:,slices));
extraData.b1        = single(b1(:,:,slices));
DIMWI=0;
if DIMWI==1
    extraData.IWF       = single(iwf(:,:,slices));
    extraData.theta     = single(theta(:,:,slices,:));
    extraData.ff        = single(ff(:,:,slices,:)./sum(ff(:,:,slices,:),4));
else
    % extraData.IWF       = [];
    % extraData.theta     = [];
    % extraData.ff        = [];

end

% keyboard
[out_askadam_mcr]    = objGPU.estimate(img(:,:,slices,:,:), mask(:,:,slices,:), extraData, fitting);




%% export to NIFTI
gre_basename    = [input.subj_label '_' output.acq_str '_' input.run_label];

imgParam.output_dir = fullfile(output.derivative_MWI_dir, 'MCR', PreProcessing,[ 'using_',num2str(nFA),'_flipangle'],['quadraticW']);
imgParam.output_filename    = [gre_basename '_MEGRE_MWI-MCR_',num2str(nFA),'FA'];
% imgParam.output_filename    = [gre_basename '_MEGRE_MWI-MCRgpuBMC_',num2str(nFA),'FA_'];
% imgParam.output_filename    = [gre_basename '_MEGRE_MWI-MCRgpuBMC_LR0001_',num2str(nFA),'FA_'];
% imgParam.output_filename    = [gre_basename '_MEGRE_MWI-MCRgpuBMC_AmpcorrLR001_',num2str(nFA),'FA_'];

% mask = [];
estimates = []; resnorm = []; iterations=[]; exitflag = [];
S0_MW = []; S0_IW=[]; S0_EW=[]; R2s_MW = []; R2s_IW = []; R2s_EW = []; Freq_MW=[]; Freq_IW=[];
T1_IEW = []; kiewm = []; Freq_BKG = []; pini = [];
%%
save_nii_quick(nii,out_askadam_mcr.final.MWF*100,	fullfile(imgParam.output_dir, [imgParam.output_filename '_MWFmap.nii.gz']));
save_nii_quick(nii,out_askadam_mcr.final.MWF.*out_askadam_mcr.final.S0,           	fullfile(imgParam.output_dir, [imgParam.output_filename '_M0map-myelinwater.nii.gz']));
save_nii_quick(nii,(1-out_askadam_mcr.final.MWF).*out_askadam_mcr.final.S0,           fullfile(imgParam.output_dir, [imgParam.output_filename '_M0map-freewater.nii.gz']));
save_nii_quick(nii,out_askadam_mcr.final.R2sMW,            fullfile(imgParam.output_dir, [imgParam.output_filename '_R2starmap-myelinwater.nii.gz']));
save_nii_quick(nii,out_askadam_mcr.final.R2sIW,            fullfile(imgParam.output_dir, [imgParam.output_filename '_R2starmap-intraaxonal.nii.gz']));
save_nii_quick(nii,1./out_askadam_mcr.final.R1IEW,            fullfile(imgParam.output_dir, [imgParam.output_filename '_T1map-freewater.nii.gz']));
save_nii_quick(nii,out_askadam_mcr.final.R1IEW,            fullfile(imgParam.output_dir, [imgParam.output_filename '_R1map-freewater.nii.gz']));
save_nii_quick(nii,out_askadam_mcr.final.kIEWM,             fullfile(imgParam.output_dir, [imgParam.output_filename '_exchangerate-freewatertomyelinwater.nii.gz']));
save_nii_quick(nii,out_askadam_mcr.final.dfreqBKG,          fullfile(imgParam.output_dir, [imgParam.output_filename '_Frequencymap-background.nii.gz']));
save_nii_quick(nii,extraData.pini +out_askadam_mcr.final.dpini,              fullfile(imgParam.output_dir, [imgParam.output_filename '_Initialphase.nii.gz']));
% save_nii_quick(nii,exitflag,          fullfile(imgParam.output_dir, [imgParam.output_filename '_exitflag.nii.gz']));
% save_nii_quick(nii,iterations,        fullfile(imgParam.output_dir, [imgParam.output_filename '_numiterations.nii.gz']));
save_nii_quick(nii,mask,       fullfile(imgParam.output_dir, [imgParam.output_filename '_mask_fittedvoxel.nii.gz']));
% MWF = S0_MW./(S0_IW+S0_EW+S0_MW)*100;

% %%
% for kz = 1:dims(3)
%     delete(fullfile(imgParamCell{kz}.output_dir,imgParamCell{kz}.output_filename )) ;
% end



function [fixed_params] = get_default_GPUfixParam(input,sepia_header)
kappa_mw                = 0.36; % Jung, NI., myelin water density
kappa_iew               = 0.86; % Jung, NI., intra-/extra-axonal water density
fixed_params.B0     	= sepia_header{end}.B0;    % field strength, in tesla
fixed_params.rho_mw    	= kappa_mw/kappa_iew; % relative myelin water density
fixed_params.E      	= 0.02; % exchange effect in signal phase, in ppm
fixed_params.x_i      	= -0.1; % myelin isotropic susceptibility, in ppm
fixed_params.x_a      	= -0.1; % myelin anisotropic susceptibility, in ppm
fixed_params.B0dir      = sepia_header{end}.B0_dir;
fixed_params.t1_mw      = 234e-3;

function [fitting] = get_default_GPUfitParam(input)
fitting = [];
fitting.Nepoch              = 4000;
fitting.initialLearnRate    = 0.001;    %    start from 0.001 % 0.01 I could not get MWF betweeen 0 and 6...
fitting.decayRate           = 0;
fitting.convergenceValue    = -inf% 1e-8;     %   -inf in this case it keeps running
fitting.tol                 = 1e-8;
fitting.display             = false;
fitting.lossFunction        = 'l2';
fitting.start               = 'prior';   
fitting.patience            = 20;       %   it waits n iterations after getting to a local minimum


fitting.DIMWI.isFitIWF      = 1;
fitting.DIMWI.isFitFreqMW   = 1;
fitting.DIMWI.isFitFreqIW   = 1;
fitting.DIMWI.isFitR2sEW    = 1;
fitting.isFitExchange       = 1;
fitting.isEPG               = 0;

