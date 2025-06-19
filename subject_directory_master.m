
%% major directories
proj_dir        = '/project/3055010.04/RunningProjects/MyelinWaterImaging/';
bids_dir        = fullfile(proj_dir, 'bidsSiemensVariantsNew');
derivative_dir  = fullfile(bids_dir, 'derivatives');
code_dir        = fullfile(bids_dir, 'code');

%% bidscoiner dicm2niix output directories
% bids_ses_dir        = fullfile(bids_dir, subj_label, sess_label);
bids_ses_dir        = fullfile(bids_dir, subj_label);
converted_dwi_dir   = fullfile(bids_ses_dir, 'dwi');
converted_anat_dir  = fullfile(bids_ses_dir, 'anat');
% converted_b1_dir    = fullfile(derivative_dir, 'SIEMENS', subj_label, sess_label, 'fmap');
converted_b1_dir    = fullfile(bids_ses_dir, 'fmap');

%% derivative top directories (for futher sub structure)
% if isvarname('sess_label')
%     derivative_ANTs_dir     = fullfile(derivative_dir, 'ANTs',      subj_label, sess_label);
%     derivative_SIEMENS_dir  = fullfile(derivative_dir, 'SIEMENS',   subj_label, sess_label);
%     derivative_FSL_dir      = fullfile(derivative_dir, 'FSL',       subj_label, sess_label);
% 
%     derivative_MP2RAGE_dir  = fullfile(derivative_dir, 'MP2RAGE', 	subj_label, sess_label);
%     % derivative_smt_dir      = fullfile(derivative_dir, 'SMT',       subj_label, sess_label);
%     derivative_ROMEO_dir    = fullfile(derivative_dir, 'ROMEO',     subj_label, sess_label);
%     derivative_DESPOT1_dir  = fullfile(derivative_dir, 'DESPOT1',	subj_label, sess_label);
%     derivative_SEPIA_dir    = fullfile(derivative_dir, 'SEPIA',     subj_label, sess_label);
%     derivative_MWI_dir      = fullfile(derivative_dir, 'MWI',       subj_label, sess_label);
%     derivative_R2star_dir  	= fullfile(derivative_dir, 'R2star',   	subj_label, sess_label);
%     derivative_R1R2s_dir  	= fullfile(derivative_dir, 'SimultaneousR1R2star',   	subj_label, sess_label);
% 
% else
    derivative_ANTs_dir     = fullfile(derivative_dir, 'ANTs',      subj_label);
%    derivative_SIEMENS_dir  = fullfile(derivative_dir, 'SIEMENS',   subj_label);
    derivative_SIEMENS_dir  = bids_ses_dir;
    derivative_FSL_dir      = fullfile(derivative_dir, 'FSL',       subj_label);
    derivative_MRI_SYNTHSEG_dir = fullfile(derivative_dir, 'MRI_SYNTHSEG',       subj_label);

    derivative_MP2RAGE_dir  = fullfile(derivative_dir, 'MP2RAGE', 	subj_label);
    % derivative_smt_dir      = fullfile(derivative_dir, 'SMT',       subj_label);
    derivative_ROMEO_dir    = fullfile(derivative_dir, 'ROMEO',     subj_label);
    derivative_DESPOT1_dir  = fullfile(derivative_dir, 'DESPOT1',	subj_label);
    derivative_SEPIA_dir    = fullfile(derivative_dir, 'SEPIA',     subj_label);
    derivative_MWI_dir      = fullfile(derivative_dir, 'MWI',       subj_label);
    derivative_R2star_dir  	= fullfile(derivative_dir, 'R2star',   	subj_label);
    derivative_R1R2s_dir  	= fullfile(derivative_dir, 'SimultaneousR1R2star',   	subj_label);
% end


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


%% %% FSL related directories
% fsl_anat_dir = fullfile(derivative_FSL_dir, 'fsl_anat');
%
% %% Diffusion related
% % Diffusion - DIPY output directories
% dipy_dir = fullfile(derivative_dir, 'dipy', subj_label, sess_label);
%
% % Diffusion - MRtrix3 output directories
% mrtrix3_dir = fullfile(derivative_dir, 'mrtrix3', subj_label, sess_label);
%
% % Diffusion - FDT output directories
% FDT_dir         = fullfile(derivative_FSL_dir, 'FDT');
% DTI_dir         = fullfile(FDT_dir, 'DTI');
% bedpostx_dir    = fullfile(FDT_dir, 'BEDPOSTX');
%
% %  Diffusion - SMT output directories
% % mcmicro_dir = fullfile(derivative_smt_dir, 'mcmicro');

%% tool directories
% tool_dir            = '~/Tools/';
% smt_tool            = fullfile(tool_dir, 'Diffusion', 'smt-0.4-linux-x64', 'bin');
% proj_python3_dir    = fullfile(code_dir, 'python3-env');
