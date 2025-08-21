function MCR_fitwrite_gpu(input, output)
% % if I am to make this into a function I need the following input:
% input.acq_str{countflip} %input acq_str of MGRE data
%
% input.derivative_SEPIA_dir % directory where the following datasets are stored
%     magn_fn         = [gre_basename '_part-mag_MEGRE_space-withinGRE.nii.gz'];
%     phase_fn        = [gre_basename '_MEGRE_space-withinGRE_part-phase_unwrapped.nii.gz'];
%     totalField_fn   = [gre_basename '_MEGRE_space-withinGRE_fieldmap.nii.gz'];
%     % mask_fn         = [gre_basename '_MEGRE_space-withinGRE_mask_refine.nii.gz'];
%     mask_fn         = [gre_basename '_MEGRE_space-withinGRE_mask_localfield.nii.gz'];
% input.derivative_FSL_dir % directory where B1 map saved in the MEGRE space is stored
% input.B1scaleFactor scaling factor for B1 map
% input.subj_label
% input.run_label
% input.MRvendor "Siemens"(deafault) or "Philips" selects the rf spoiling phase
% input.CorrectionFactorVFA  scale factor for VFA acquisition (should be some length as number of flip angles) 
% input.Configfile % config file for MCR-MWI fitting parameters

if ~isfield(output, 'MPPCAdenoise')
    output.MPPCAdenoise = 0;
end

%% load GRE data
img             = [];
unwrappedPhase  = [];
mask            = [];
totalField      = [];
nFA             = length(input.acq_str);
fa              = zeros(1, nFA);
for countflip = 1:length(input.acq_str)

    seq_SEPIA_dir = fullfile(input.derivative_SEPIA_dir, input.acq_str{countflip});

    % general GRE basename
    gre_basename    = [input.subj_label '_' input.acq_str{countflip} '_' input.run_label];

    % magnitude nifti image filename
    magn_fn         = [gre_basename '_part-mag_MEGRE_space-withinGRE.nii.gz'];
    phase_fn        = [gre_basename '_MEGRE_space-withinGRE_part-phase_unwrapped.nii.gz'];
    totalField_fn   = [gre_basename '_MEGRE_space-withinGRE_fieldmap.nii.gz'];
    mask_fn         = [gre_basename '_MEGRE_space-withinGRE_mask_localfield.nii.gz'];
    sepia_header_fn = [gre_basename '_header.mat'];
    img             = cat(5, img, spm_read_vols(spm_vol(fullfile(input.derivative_SEPIA_dir, magn_fn))));
    unwrappedPhase  = cat(5, unwrappedPhase, spm_read_vols(spm_vol(fullfile(seq_SEPIA_dir, phase_fn))));
    sepia_header(countflip)	= load(fullfile(input.derivative_SEPIA_dir, sepia_header_fn));
    mask            = cat(5, mask, spm_read_vols(spm_vol(fullfile(seq_SEPIA_dir, mask_fn))));
    totalField      = cat(4, totalField, spm_read_vols(spm_vol(fullfile(seq_SEPIA_dir, totalField_fn))));
    fa(countflip)   = sepia_header(countflip).FA;
    tr              = sepia_header(countflip).TR; % note that here there is an assumption that all protocols have the same TR
    te              = sepia_header(countflip).TE;

end
dims = size(img);
% B1 info
hdr                     = spm_vol(fullfile(input.derivative_FSL_dir, [input.subj_label '_acq-famp_run-1_TB1TFLProtocolSpace.nii.gz']));
true_flip_angle         = spm_read_vols(hdr);
b1                      = true_flip_angle / input.B1scaleFactor;
clear true_flip_angle

figure
Orthoview2(sum(mask, 5), [], [], 'tight')

%% obtain the initial estimation of the initial B1 phase
img = img .* exp(1i*unwrappedPhase);

mask = all(mask, 5) & all(~isnan(img), [4 5]);  % true where mask is true along 5th dim and no NaNs along 4th/5th dims

if ~isfield(input, 'Configfile')
    input.Configfile = [];
    fixed_params     = get_default_GPUfixParam(sepia_header);
    fitting          = get_default_GPUfitParam();
else
    % load the config file
    if isfile(input.Configfile)
        run(input.Configfile);
    else
        error('Config file not found');
    end
end

%% denoising part
if output.MPPCAdenoise == 1
    denoised = denoise(reshape(img, [dims(1:3) prod(dims(4:5))]), [5 5 5], mask);
    img = reshape(denoised, dims);
    clear denoised
    PreProcessing = 'MPPCAdenoising';
else
    if output.MPPCAdenoise == 2
        denoised = denoise(reshape(img, [dims(1:3) prod(dims(4:5))]), [3 3 3], mask);
        img = reshape(denoised, dims);
        clear denoised
        PreProcessing ='MPPCAdenoising3';
    else
        PreProcessing ='no_preprocessing';
    end
end

pini = squeeze(unwrappedPhase(:,:,:,1,:)) - 2*pi*totalField .* sepia_header(end).TE(1);

Debug = 0;
if Debug==1
    % figure
    for k = 1:nFA
        figure(100)
        tiledlayout(ceil(sqrt(nFA)), ceil(sqrt(nFA)))
        nexttile(1)
        nexttile
        Orthoview2(pini(:,:,:,k), [], [], 'tight')
        title(['Initial phase map for FA = ' num2str(fa(k))])
    end
end

pini = polyfit3D_NthOrder(double(mean(pini(:,:,:,1:(end-1)), 4)), prod(mask,4), 6);

clear unwrappedPhase

% data normalisation
if isfield(input,'CorrectionFactorVFA')
    nFA = length(input.CorrectionFactorVFA);
    for countflip= 1:nFA
        img(:,:,:,:,countflip) = img(:,:,:,:,countflip) * (input.CorrectionFactorVFA(countflip));
    end
end

slices = 1:dims(3);
extraData = [];
extraData.freqBKG   = single(squeeze(totalField(:,:,slices,:,:)) / (gpuGREMWI.gyro*fixed_params.B0)); % in ppm
extraData.pini      = single(pini(:,:,slices));
extraData.b1        = single(b1(:,:,slices));
DIMWI = 0;
if DIMWI == 1
    extraData.IWF   = single(iwf(:,:,slices));
    extraData.theta = single(theta(:,:,slices,:));
    extraData.ff    = single(ff(:,:,slices,:)./sum(ff(:,:,slices,:),4));
end

objGPU          = gpuMCRMWI(te,tr,fa,fixed_params);
out_askadam_mcr = objGPU.estimate(img(:,:,slices,:,:), mask(:,:,slices,:), extraData, fitting);


%% export to NIFTI
gre_basename             = [input.subj_label '_' output.acq_str '_' input.run_label];
imgParam.output_dir      = fullfile(output.derivative_MWI_dir, 'MCR', PreProcessing, ['using_' num2str(nFA) '_flipangle'], 'quadraticW');
imgParam.output_filename = [gre_basename '_MEGRE_MWI-MCR_' num2str(nFA) 'FA'];

hdr.dim = dims(1:3);
spm_write_vol_gz(hdr, out_askadam_mcr.final.MWF * 100,	                         fullfile(imgParam.output_dir, [imgParam.output_filename '_MWFmap.nii.gz']));
spm_write_vol_gz(hdr, out_askadam_mcr.final.MWF .* out_askadam_mcr.final.S0,     fullfile(imgParam.output_dir, [imgParam.output_filename '_M0map-myelinwater.nii.gz']));
spm_write_vol_gz(hdr, (1-out_askadam_mcr.final.MWF) .* out_askadam_mcr.final.S0, fullfile(imgParam.output_dir, [imgParam.output_filename '_M0map-freewater.nii.gz']));
spm_write_vol_gz(hdr, out_askadam_mcr.final.R2sMW,                               fullfile(imgParam.output_dir, [imgParam.output_filename '_R2starmap-myelinwater.nii.gz']));
spm_write_vol_gz(hdr, out_askadam_mcr.final.R2sIW,                               fullfile(imgParam.output_dir, [imgParam.output_filename '_R2starmap-intraaxonal.nii.gz']));
spm_write_vol_gz(hdr, 1 ./ out_askadam_mcr.final.R1IEW,                          fullfile(imgParam.output_dir, [imgParam.output_filename '_T1map-freewater.nii.gz']));
spm_write_vol_gz(hdr, out_askadam_mcr.final.R1IEW,                               fullfile(imgParam.output_dir, [imgParam.output_filename '_R1map-freewater.nii.gz']));
spm_write_vol_gz(hdr, out_askadam_mcr.final.kIEWM,                               fullfile(imgParam.output_dir, [imgParam.output_filename '_exchangerate-freewatertomyelinwater.nii.gz']));
% spm_write_vol_gz(hdr, out_askadam_mcr.final.dfreqBKG,                          fullfile(imgParam.output_dir, [imgParam.output_filename '_Frequencymap-background.nii.gz']));  % = 4D, spm_write_vol can only handle a maximum of 3 dimensions
spm_write_vol_gz(hdr, extraData.pini + out_askadam_mcr.final.dpini,              fullfile(imgParam.output_dir, [imgParam.output_filename '_Initialphase.nii.gz']));
spm_write_vol_gz(hdr, mask,                                                      fullfile(imgParam.output_dir, [imgParam.output_filename '_mask_fittedvoxel.nii.gz']));


function fixed_params = get_default_GPUfixParam(sepia_header)
kappa_mw                = 0.36; % Jung, NI., myelin water density
kappa_iew               = 0.86; % Jung, NI., intra-/extra-axonal water density
fixed_params.B0     	= sepia_header(end).B0;    % field strength, in tesla
fixed_params.rho_mw    	= kappa_mw/kappa_iew; % relative myelin water density
fixed_params.E      	= 0.02; % exchange effect in signal phase, in ppm
fixed_params.x_i      	= -0.1; % myelin isotropic susceptibility, in ppm
fixed_params.x_a      	= -0.1; % myelin anisotropic susceptibility, in ppm
fixed_params.B0dir      = sepia_header(end).B0_dir;
fixed_params.t1_mw      = 234e-3;


function fitting = get_default_GPUfitParam()
fitting = [];
fitting.Nepoch              = 4000;
fitting.initialLearnRate    = 0.001;    %    start from 0.001 % 0.01 I could not get MWF betweeen 0 and 6...
fitting.decayRate           = 0;
fitting.convergenceValue    = -inf;     % 1e-8;     %   -inf in this case it keeps running
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
