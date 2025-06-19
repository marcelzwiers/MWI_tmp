clear;

%% Subject info and directories
% subj	= '002';
run     = '1';
% ses     = 'mri01';

subj_label  = ['sub-' subj];
sess_label  = ['ses-' ses];
run_label   = ['run-' run];

% load path
addpath('/project/3015069.05/bids/code/');
load_module_despot1;
load_module_epg_epgx;
load_module_mwi;
load_module_sepia;

% call diretcories
subject_directory_master

subj_script_dir = fullfile(code_dir, subj_label, sess_label);

% GRE Protocols and flip angles
% protocol    = { 'TR38NTE12', 'TR50NTE15', 'TR55NTE13' };
% protocol    = { 'TR38NTE12'};
flipAngle   = { 'FA5', 'FA10', 'FA20', 'FA50', 'FA70' };

%% load GRE data
for kp = protocol
    
    prot_ANTs_dir   = fullfile(ANTs_gre_within_dir, kp{1});
    prot_SEPIA_dir  = fullfile(derivative_SEPIA_dir, kp{1});

    % load GRE data
    counter         = 0;
    img             = [];
    unwrappedPhase  = [];
    totalField      = [];
    sepia_header    = [];
    fa              = zeros(1,length(flipAngle));
    for kf = flipAngle
        counter = counter + 1;
        seq_SEPIA_dir   = fullfile(prot_SEPIA_dir, kf{1});
        
        % general GRE basename
        seq             = [kp{1} kf{1}];
        acq_label       = ['acq-' seq];
        gre_basename    = [subj_label '_' sess_label '_' acq_label '_' run_label];
        
        % magnitude nifti image filename
        magn_fn         = [gre_basename '_part-mag_MEGRE_space-withinGRE.nii.gz'];
        phase_fn        = [gre_basename '_MEGRE_space-withinGRE_EchoCombine-OW_part-phase_unwrapped.nii.gz'];
        totalField_fn   = [gre_basename '_MEGRE_space-withinGRE_EchoCombine-OW_fieldmap.nii.gz'];
        
        sepia_header_fn = [gre_basename '_MEGRE_sepia_header.mat'];
        
        nii                     = load_untouch_nii(fullfile(ANTs_gre_apply_within_dir, magn_fn));
        img                     = cat(5,img,nii.img);
        nii = load_untouch_nii(fullfile(seq_SEPIA_dir, phase_fn));
        unwrappedPhase = cat(5,unwrappedPhase,nii.img);
        sepia_header{counter}	= load(fullfile(prot_SEPIA_dir, sepia_header_fn));
        
        totalField = cat(4,totalField, load_nii_img_only(fullfile(seq_SEPIA_dir, totalField_fn)));
        
        fa(counter)   	= sepia_header{counter}.FA;
        tr              = sepia_header{counter}.TR;

    end
    
    % B1 info
    true_flip_angle_fn      = [subj_label '_' sess_label '_acq-famp_run-1_TB1TFL_space-withinGRE-' kp{1} '.nii.gz'];
    true_flip_angle_json	= [subj_label '_' sess_label '_acq-famp_run-1_TB1TFL.json'];
    
    true_flip_angle         = load_nii_img_only( fullfile( ANTs_b1_apply_gre_dir, true_flip_angle_fn));
    b1_header               = jsondecode( fileread( fullfile( converted_b1_dir,true_flip_angle_json)));
    
    b1                      = true_flip_angle / 10 / b1_header.FlipAngle;
    
    mask_fn         = [subj_label '_' sess_label '_acq-' kp{1} '_' run_label '_brain_mask.nii.gz'];
    mask_filename   = fullfile(prot_SEPIA_dir, mask_fn);
    mask            = load_nii_img_only(mask_filename);
    
end

clear true_flip_angle

%%
img = img .* exp(1i*unwrappedPhase);

mask_nonnan = min(min( ~isnan(img), [], 5), [],4);

mask = and(mask,mask_nonnan);

pini = unwrappedPhase - 2*pi*permute(totalField,[1 2 3 5 4]) .* permute(sepia_header{end}.TE(:),[2 3 4 1 5]);
pini = polyfit3D_NthOrder( double(mean(pini(:,:,:,1,1:3), 5)), mask, 6);

clear unwrappedPhase

%%
% setup algorithm parameters
algoParam.isInvivo      = true;
algoParam.isParallel    = false;
algoParam.DEBUG         = false;
algoParam.isNormData    = false; % normalised outside the function
% fitting option
algoParam.maxIter       = 200;
algoParam.fcnTol        = 1e-4;
algoParam.stepTol       = 1e-4;
% residual option
algoParam.numMagn       = 0;
algoParam.isWeighted    = true;
algoParam.weightMethod  = '1stEcho';
% T1 model
algoParam.isExchange    = 1;
algoParam.isEPG         = 1;
algoParam.npulse        = 50;
algoParam.rfphase       = 50;
algoParam.isT1mw        = false;
algoParam.T1mw          = 234e-3;
% No DIMWI
algoParam.DIMWI.isVic       = false;
algoParam.DIMWI.isR2sEW     = false;
algoParam.DIMWI.isFreqMW    = false;
algoParam.DIMWI.isFreqIW    = false;
% initial guess
algoParam.advancedStarting = 'robust';

mask_tmp = mask>0;
[scaleFactor, ~] = mwi_image_normalisation(img, mask_tmp);
clear mask_tmp

% fixed parameters
kappa_mw    = 0.36; % jung
kappa_iew   = 0.86; % jung
imgParam.b0          = sepia_header{end}.B0;
imgParam.b0dir       = sepia_header{end}.B0_dir;
imgParam.rho_mw      = kappa_mw/kappa_iew;
imgParam.E           = 0.02;
imgParam.x_i         = -0.1;
imgParam.x_a         = -0.1;

imgParam.te         = sepia_header{end}.TE;
imgParam.tr         = tr;
imgParam.fa         = fa;
imgParam.img        = img(:,:,slice,:,:)/scaleFactor;
imgParam.mask       = mask(:,:,slice);
imgParam.fieldmap   = totalField(:,:,slice,:);
imgParam.pini       = pini(:,:,slice,:);
imgParam.b1map      = b1(:,:,slice);

clear img mask totalField pini b1

% output
% general GRE basename
seq             = [protocol{1}];
acq_label       = ['acq-' seq];
gre_basename    = [subj_label '_' sess_label '_' acq_label '_' run_label];

imgParam.output_dir         = fullfile(derivative_MWI_dir, 'MCR', 'no_preprocessing', 'using_5_flipangle',['tmp_' seq]);
imgParam.output_filename    = [gre_basename '_MEGRE_MWI-MCR_5FA_slice-' num2str(slice)];

fitres_mcr_ow = mwi_3cx_2R1R2s_dimwi(algoParam,imgParam);

fitres_mcr_ow.estimates(:,:,:,1:3) = fitres_mcr_ow.estimates(:,:,:,1:3) * scaleFactor;
fitres_mcr_ow.S0_IW                = fitres_mcr_ow.S0_IW * scaleFactor;
fitres_mcr_ow.S0_EW                = fitres_mcr_ow.S0_EW * scaleFactor;
fitres_mcr_ow.S0_MW                = fitres_mcr_ow.S0_MW * scaleFactor;
save(fullfile(imgParam.output_dir,imgParam.output_filename),'fitres_mcr_ow','scaleFactor','-append')

% %% export to NIFTI
% save_nii_quick(nii,ComputeMWF(fitres_mcrdimwi_ow)*100,	fullfile(imgParam.output_dir, [imgParam.output_filename '_MWFmap.nii.gz']));
% save_nii_quick(nii,fitres_mcrdimwi_ow.S0_MW,           	fullfile(imgParam.output_dir, [imgParam.output_filename '_M0map-myelinwater.nii.gz']));
% save_nii_quick(nii,fitres_mcrdimwi_ow.S0_IEW,           fullfile(imgParam.output_dir, [imgParam.output_filename '_M0map-freewater.nii.gz']));
% save_nii_quick(nii,fitres_mcrdimwi_ow.R2s_MW,            fullfile(imgParam.output_dir, [imgParam.output_filename '_R2starmap-myelinwater.nii.gz']));
% save_nii_quick(nii,fitres_mcrdimwi_ow.R2s_IW,            fullfile(imgParam.output_dir, [imgParam.output_filename '_R2starmap-intraaxonal.nii.gz']));
% save_nii_quick(nii,fitres_mcrdimwi_ow.T1_IEW,            fullfile(imgParam.output_dir, [imgParam.output_filename '_T1map-freewater.nii.gz']));
% save_nii_quick(nii,fitres_mcrdimwi_ow.kiewm,             fullfile(imgParam.output_dir, [imgParam.output_filename '_exchangerate-freewatertomyelinwater.nii.gz']));
% save_nii_quick(nii,fitres_mcrdimwi_ow.Freq_BKG,          fullfile(imgParam.output_dir, [imgParam.output_filename '_Frequencymap-background.nii.gz']));
% save_nii_quick(nii,fitres_mcrdimwi_ow.pini,              fullfile(imgParam.output_dir, [imgParam.output_filename '_Initialphase.nii.gz']));
% save_nii_quick(nii,fitres_mcrdimwi_ow.exitflag,          fullfile(imgParam.output_dir, [imgParam.output_filename '_exitflag.nii.gz']));
% save_nii_quick(nii,fitres_mcrdimwi_ow.iterations,        fullfile(imgParam.output_dir, [imgParam.output_filename '_numiterations.nii.gz']));
% save_nii_quick(nii,fitres_mcrdimwi_ow.mask_fitted,       fullfile(imgParam.output_dir, [imgParam.output_filename '_mask_fittedvoxel.nii.gz']));
