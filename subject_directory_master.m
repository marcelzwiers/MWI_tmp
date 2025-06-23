%% major directories
derivative_dir  = fullfile(bids_dir, 'derivatives');

%% bidscoiner dicm2niix output directories
% bids_ses_dir        = fullfile(bids_dir, subj_label, sess_label);
bids_ses_dir        = fullfile(bids_dir, subj_label);
converted_dwi_dir   = fullfile(bids_ses_dir, 'dwi');
converted_anat_dir  = fullfile(bids_ses_dir, 'anat');
converted_b1_dir    = fullfile(bids_ses_dir, 'fmap');

derivative_ANTs_dir     = fullfile(derivative_dir, 'ANTs',      subj_label);
derivative_FSL_dir      = fullfile(derivative_dir, 'FSL',       subj_label);
derivative_MRI_SYNTHSEG_dir = fullfile(derivative_dir, 'MRI_SYNTHSEG',       subj_label);
derivative_MP2RAGE_dir  = fullfile(derivative_dir, 'MP2RAGE', 	subj_label);
derivative_ROMEO_dir    = fullfile(derivative_dir, 'ROMEO',     subj_label);
derivative_DESPOT1_dir  = fullfile(derivative_dir, 'DESPOT1',	subj_label);
derivative_SEPIA_dir    = fullfile(derivative_dir, 'SEPIA',     subj_label);
derivative_MWI_dir      = fullfile(derivative_dir, 'MWI',       subj_label);
derivative_R2star_dir  	= fullfile(derivative_dir, 'R2star',   	subj_label);
derivative_R1R2s_dir  	= fullfile(derivative_dir, 'SimultaneousR1R2star',   	subj_label);

%% ANTs related directories
ANTs_gre_within_dir            = fullfile(derivative_ANTs_dir, 'gre_within_protocol_transformation');
ANTs_gre_apply_within_dir      = fullfile(derivative_ANTs_dir, 'gre_within_protocol_space');

FSL_gre_within_dir             = fullfile(derivative_FSL_dir, 'gre_within_protocol_transformation');
FSL_gre_apply_within_dir       = fullfile(derivative_FSL_dir, 'gre_within_protocol_space');

% ANTs_gre_across_dir             = fullfile(derivative_ANTs_dir, 'gre_across_protocols_transformation');
% ANTs_gre_apply_registration_dir = fullfile(derivative_ANTs_dir, 'gre_across_protocols_space');

ANTs_b1_b12gre_dir              = fullfile(derivative_ANTs_dir, 'b1_2_gre_transformation');
ANTs_b1_apply_gre_dir           = fullfile(derivative_ANTs_dir, 'b1_2_gre_space');
prot_ANTs_dir = ANTs_b1_apply_gre_dir;
% ANTs_b1_b12mp2rage_dir          = fullfile(derivative_ANTs_dir, 'b1_2_mp2rage_transformation');
% ANTs_b1_apply_mp2rage_dir       = fullfile(derivative_ANTs_dir, 'b1_2_mp2rage_space');
