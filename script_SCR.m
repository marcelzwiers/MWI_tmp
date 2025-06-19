% example script to process MCR-MWI on a local computer with multi-threads
% The MCR-MWI processing scripts inside the subject folders were setted up
% slightly differently because of the parallelisation on HPC vs local PC

%% load GRE data
% for kp = protocol

% prot_SEPIA_dir  = fullfile(derivative_SEPIA_dir, kp{1});
% load GRE data
S0              = [];
R2s             = [];
Chi             = [];
mask            = [];
totalField      = [];
sepia_header    = [];
fa              = zeros(1,length(prot{1}.flip));
for countflip= 1:length(prot{1}.flip)
    
    seq_SEPIA_dir = fullfile(derivative_SEPIA_dir,prot{1}.acq_str{countflip});
    
    % general GRE basename
    gre_basename    = [subj_label '_' prot{1}.acq_str{countflip} '_' run_label];
    
    % magnitude nifti image filename
    S0_fn           = [gre_basename '_MEGRE_space-withinGRE_S0map.nii.gz'];
    R2s_fn          = [gre_basename '_MEGRE_space-withinGRE_R2starmap.nii.gz'];
    Chi_fn          = [gre_basename '_MEGRE_space-withinGRE_Chimap.nii.gz'];
    mask_fn         = [gre_basename '_MEGRE_space-withinGRE_mask_localfield.nii.gz'];
    
    sepia_header_fn = [gre_basename '_header.mat'];
    
    nii             = load_untouch_nii(fullfile(seq_SEPIA_dir, S0_fn));
    S0              = cat(5, S0, nii.img);
    nii             = load_untouch_nii(fullfile(seq_SEPIA_dir, R2s_fn));
    R2s             = cat(5, R2s, nii.img);
    nii             = load_untouch_nii(fullfile(seq_SEPIA_dir, Chi_fn));
    Chi             = cat(5, Chi, nii.img);
    sepia_header{countflip}	= load(fullfile(derivative_SEPIA_dir, sepia_header_fn));
    mask            = cat(5, mask, load_nii_img_only(fullfile(seq_SEPIA_dir, mask_fn)));
    
    fa(countflip)   = sepia_header{countflip}.FA;
    tr              = sepia_header{countflip}.TR; % note that here there is an assumption that all protocols have the same TR
end

% B1 info
% true_flip_angle_fn      = [subj_label '_' sess_label '_acq-famp_run-1_TB1TFL_space-withinGRE-' kp{1} '.nii.gz'];
true_flip_angle_json	= [subj_label '_acq-famp_run-1_TB1TFL.json'];
true_flip_angle_fn      = [subj_label '_acq-famp_run-1_TB1TFLProtocolSpace.nii.gz'];
true_flip_angle         = load_nii_img_only( fullfile(derivative_FSL_dir, true_flip_angle_fn));
b1_header               = jsondecode( fileread( fullfile( converted_b1_dir,true_flip_angle_json)));

b1                      = true_flip_angle / 10 / b1_header.FlipAngle;
figure
Orthoview2(sum(mask,5), [], [], 'tight')
mask = prod(mask,5);
% mask_fn         = [subj_label '_' sess_label '_acq-' kp{1} '_' run_label '_brain_mask.nii.gz'];

sts = mkdir(derivative_R1R2s_dir);
gre_basename    = [subj_label '_' prot{1}.rec '_' run_label];

R2smean = sum(S0.^2 .* R2s,5)./sum(S0.^2 ,5);

save_nii_quick(nii,R2smean.*mask, fullfile(derivative_R1R2s_dir,[ gre_basename, '_MEGRE_space-withinGRE_R2starmap.nii.gz']));

Chimean = sum(S0.^2 .* Chi,5)./sum(S0.^2 ,5);

save_nii_quick(nii,Chimean.*mask, fullfile(derivative_R1R2s_dir,[ gre_basename, '_MEGRE_space-withinGRE_Chimap.nii.gz']));


[T1, M0] = despot1_mapping(double(squeeze(S0(:,:,:,1:2))), fa(1:2), tr, mask, b1);

% [T1_noB1, M0_noB1] = despot1_mapping(double(squeeze(S0(:,:,:,1:2))),fa(1:3),tr);
R1 = double(1./double(T1).*mask) * 1000;
R1(isnan(R1)) = 0;
R1(isinf(R1)) = 0;

save_nii_quick(nii,R1, fullfile(derivative_R1R2s_dir, [gre_basename '_MEGRE_space-withinGRE_R1map.nii.gz']));
save_nii_quick(nii,M0.*mask, fullfile(derivative_R1R2s_dir, [gre_basename '_MEGRE_space-withinGRE_M0map.nii.gz']));


% derivative_R1R2s_dir
% %%
% save_nii_quick(nii,S0_MW./(S0_IW+S0_EW+S0_MW)*100,	fullfile(imgParam.output_dir, [imgParam.output_filename '_MWFmap.nii.gz']));
