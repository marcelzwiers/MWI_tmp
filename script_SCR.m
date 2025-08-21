%% load GRE data
S0           = [];
R2s          = [];
Chi          = [];
mask         = [];
fa           = zeros(1, length(prot.flip));
for countflip = 1:length(prot.flip)
    
    seq_SEPIA_dir = fullfile(derivative_SEPIA_dir, prot.acq_str{countflip});
    
    % general GRE basename
    gre_basename    = [subj_label '_' prot.acq_str{countflip} '_' run_label];
    
    % magnitude nifti image filename
    S0_fn           = [gre_basename '_MEGRE_space-withinGRE_S0map.nii.gz'];
    R2s_fn          = [gre_basename '_MEGRE_space-withinGRE_R2starmap.nii.gz'];
    Chi_fn          = [gre_basename '_MEGRE_space-withinGRE_Chimap.nii.gz'];
    mask_fn         = [gre_basename '_MEGRE_space-withinGRE_mask_localfield.nii.gz'];
    sepia_header_fn = [gre_basename '_header.mat'];
    
    hdr             = spm_vol(fullfile(seq_SEPIA_dir, S0_fn));
    S0              = cat(5, S0, spm_read_vols(hdr));
    R2s             = cat(5, R2s, spm_read_vols(spm_vol(fullfile(seq_SEPIA_dir, R2s_fn))));
    Chi             = cat(5, Chi, spm_read_vols(spm_vol(fullfile(seq_SEPIA_dir, Chi_fn))));
    sepia_header{countflip}	= load(fullfile(derivative_SEPIA_dir, sepia_header_fn));
    mask            = cat(5, mask, spm_read_vols(spm_vol(fullfile(seq_SEPIA_dir, mask_fn))));
    
    fa(countflip)   = sepia_header{countflip}.FA;
    tr              = sepia_header{countflip}.TR; % note that here there is an assumption that all protocols have the same TR
end

% B1 info
true_flip_angle_json	= [subj_label '_acq-famp_run-1_TB1TFL.json'];
true_flip_angle_fn      = [subj_label '_acq-famp_run-1_TB1TFLProtocolSpace.nii.gz'];
true_flip_angle         = load_nii_img_only(fullfile(derivative_FSL_dir, true_flip_angle_fn));
b1_header               = jsondecode(fileread(fullfile(converted_b1_dir, true_flip_angle_json)));
b1                      = true_flip_angle / 10 / b1_header.FlipAngle;
figure
Orthoview2(sum(mask,5), [], [], 'tight')
mask = all(mask, 5);

sts = mkdir(derivative_R1R2s_dir);
gre_basename = [subj_label '_' prot.rec '_' run_label];

R2smean = sum(S0.^2 .* R2s,5)./sum(S0.^2 ,5);
Chimean = sum(S0.^2 .* Chi,5)./sum(S0.^2 ,5);

spm_write_vol_gz(hdr, R2smean.*mask, fullfile(derivative_R1R2s_dir, [gre_basename '_MEGRE_space-withinGRE_R2starmap.nii.gz']));
spm_write_vol_gz(hdr, Chimean.*mask, fullfile(derivative_R1R2s_dir, [gre_basename '_MEGRE_space-withinGRE_Chimap.nii.gz']));

[T1, M0] = despot1_mapping(double(squeeze(S0(:,:,:,1:2))), fa(1:2), tr, mask, b1);

R1 = (mask ./ T1) * 1000;
R1(~isfinite(R1)) = 0;          % set NaN and Inf to 0

spm_write_vol_gz(hdr, R1,       fullfile(derivative_R1R2s_dir, [gre_basename '_MEGRE_space-withinGRE_R1map.nii.gz']));
spm_write_vol_gz(hdr, M0.*mask, fullfile(derivative_R1R2s_dir, [gre_basename '_MEGRE_space-withinGRE_M0map.nii.gz']));
