load_module_mwi;
clear;
%% Subject info and directories
subj	= '003';
run     = '1';
ses     = 'mri02';

subj_label  = ['sub-' subj];
sess_label  = ['ses-' ses];
run_label   = ['run-' run];

% load path
addpath('/project/3015069.05/bids/code/');
load_module_mwi;

% call diretcories
subject_directory_master

subj_script_dir = fullfile(code_dir, subj_label, sess_label);

% GRE Protocols and flip angles
% protocol    = { 'TR38NTE12', 'TR50NTE15', 'TR55NTE13' };
protocol    = { 'TR38NTE12'};
flipAngle   = { 'FA5', 'FA10', 'FA20', 'FA50', 'FA70' };

seq             = [protocol{1}];
acq_label       = ['acq-' seq];
gre_basename    = [subj_label '_' sess_label '_' acq_label '_' run_label];


tmp_dir         = fullfile(derivative_MWI_dir, 'MCR', 'no_preprocessing', 'using_5_flipangle',['tmp_' seq]);

prot_SEPIA_dir  = fullfile(derivative_SEPIA_dir, protocol{1});
mask_fn         = [subj_label '_' sess_label '_acq-' protocol{1} '_' run_label '_brain_mask.nii.gz'];
mask_filename   = fullfile(prot_SEPIA_dir, mask_fn);

nii             = load_untouch_nii(mask_filename);

output_dir = fullfile(derivative_MWI_dir, 'MCR', 'no_preprocessing', 'using_5_flipangle');
output_filename = [gre_basename '_MEGRE_MWI-MCR_5FA'];
    
%%

mask = [];
estimates = []; resnorm = []; iterations=[]; exitflag = [];
S0_MW = []; S0_IW=[]; S0_EW=[]; R2s_MW = []; R2s_IW = []; R2s_EW = []; Freq_MW=[]; Freq_IW=[];
T1_IEW = []; kiewm = []; Freq_BKG = []; pini = [];
for kz = 1:160
    tmp_output_fn = [gre_basename '_MEGRE_MWI-MCR_5FA_slice-' num2str(kz)];
    load( fullfile(tmp_dir, tmp_output_fn))
    mask        = cat(3,mask,fitres_mcr_ow.mask_fitted);
    estimates   = cat(3,estimates,fitres_mcr_ow.estimates);
    resnorm     = cat(3,resnorm,fitres_mcr_ow.resnorm);
    iterations  = cat(3,iterations,fitres_mcr_ow.iterations);
    exitflag    = cat(3,exitflag,fitres_mcr_ow.exitflag);
    S0_MW       = cat(3,S0_MW,fitres_mcr_ow.S0_MW);
    S0_IW       = cat(3,S0_IW,fitres_mcr_ow.S0_IW);
    S0_EW       = cat(3,S0_EW,fitres_mcr_ow.S0_EW);
    R2s_MW      = cat(3,R2s_MW,fitres_mcr_ow.R2s_MW);
    R2s_IW      = cat(3,R2s_IW,fitres_mcr_ow.R2s_IW);
    R2s_EW      = cat(3,R2s_EW,fitres_mcr_ow.R2s_EW);
    Freq_MW     = cat(3,Freq_MW,fitres_mcr_ow.Freq_MW);
    Freq_IW     = cat(3,Freq_IW,fitres_mcr_ow.Freq_IW);
    T1_IEW      = cat(3,T1_IEW,fitres_mcr_ow.T1_IEW);
    kiewm       = cat(3,kiewm,fitres_mcr_ow.kiewm);
    Freq_BKG    = cat(3,Freq_BKG,fitres_mcr_ow.Freq_BKG);
    pini        = cat(3,pini,fitres_mcr_ow.pini);
end
%%
fitres_mcr_ow = [];
fitres_mcr_ow.mask_fitted     = mask;
fitres_mcr_ow.estimates = estimates;
fitres_mcr_ow.resnorm  = resnorm;
fitres_mcr_ow.iterations = iterations;
fitres_mcr_ow.exitflag = exitflag;
fitres_mcr_ow.S0_MW    = S0_MW;
fitres_mcr_ow.S0_IW    = S0_IW;
fitres_mcr_ow.S0_EW    = S0_EW;
fitres_mcr_ow.R2s_MW   = R2s_MW;
fitres_mcr_ow.R2s_IW   = R2s_IW;
fitres_mcr_ow.R2s_EW   = R2s_EW;
fitres_mcr_ow.Freq_MW  = Freq_MW;
fitres_mcr_ow.Freq_IW  = Freq_IW;
fitres_mcr_ow.T1_IEW   = T1_IEW;
fitres_mcr_ow.kiewm    = kiewm;
fitres_mcr_ow.Freq_BKG = Freq_BKG;
fitres_mcr_ow.pini     = pini;
%%
save( fullfile(output_dir, output_filename),'fitres_mcr_ow');
%% export to NIFTI
save_nii_quick(nii,ComputeMWF(fitres_mcr_ow)*100,	fullfile(output_dir, [output_filename '_MWFmap.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.S0_MW,           	fullfile(output_dir, [output_filename '_M0map-myelinwater.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.S0_IW,             fullfile(output_dir, [output_filename '_M0map-intraaxonalwater.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.S0_EW,             fullfile(output_dir, [output_filename '_M0map-extraaxonalwater.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.R2s_MW,            fullfile(output_dir, [output_filename '_R2starmap-myelinwater.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.R2s_IW,            fullfile(output_dir, [output_filename '_R2starmap-intraaxonalwater.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.R2s_EW,            fullfile(output_dir, [output_filename '_R2starmap-extraaxonalwater.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.Freq_MW,           fullfile(output_dir, [output_filename '_Frequencymap-myelinwater.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.Freq_IW,           fullfile(output_dir, [output_filename '_Frequencymap-intraaxonalwater.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.T1_IEW,            fullfile(output_dir, [output_filename '_T1map-freewater.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.kiewm,             fullfile(output_dir, [output_filename '_exchangerate-freewatertomyelinwater.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.Freq_BKG,          fullfile(output_dir, [output_filename '_Frequencymap-background.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.pini,              fullfile(output_dir, [output_filename '_Initialphase.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.exitflag,          fullfile(output_dir, [output_filename '_exitflag.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.iterations,        fullfile(output_dir, [output_filename '_numiterations.nii.gz']));
save_nii_quick(nii,fitres_mcr_ow.mask_fitted,       fullfile(output_dir, [output_filename '_mask_fittedvoxel.nii.gz']));

