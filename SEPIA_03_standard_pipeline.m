%% Subject info and directories

% General algorithm parameters
algorParam = struct();
algorParam.general.isBET = 0;
algorParam.general.isRefineBrainMask = 0;
algorParam.general.isInvert = 0;
% Total field recovery algorithm parameters
algorParam.unwrap.echoCombMethod = 'Optimum weights';
algorParam.unwrap.unwrapMethod = 'SEGUE';
algorParam.unwrap.isEddyCorrect = 0;
algorParam.unwrap.isSaveUnwrappedEcho = 1;
algorParam.unwrap.excludeMaskThreshold = 0.4;
algorParam.unwrap.excludeMethod = 'Weighting map';
% Background field removal algorithm parameters
algorParam.bfr.refine_method = 'Spherical harmonic';
algorParam.bfr.refine_order = 2;
algorParam.bfr.erode_radius = 1;
algorParam.bfr.method = 'PDF';
algorParam.bfr.tol = 0.001;
algorParam.bfr.iteration = 200;
algorParam.bfr.padSize = 40;
% QSM algorithm parameters
algorParam.qsm.reference_tissue = 'CSF';
algorParam.qsm.method = 'MRI Suscep. Calc.';
algorParam.qsm.solver = 'Iterative Tikhonov';
algorParam.qsm.lambda = 0.05;
algorParam.qsm.tolerance = 0.03;

clear input output_basename mask_filename
for flip = 1:length(prot.flip)
    seq_SEPIA_dir = fullfile(derivative_SEPIA_dir,prot.acq_str{flip});
    mkdir(seq_SEPIA_dir)
    % general GRE basename
    gre_basename    = [subj_label '_' prot.acq_str{flip} '_' run_label];

    % magnitude nifti image filename
    phase_fn        = [gre_basename '_part-phase_MEGRE_space-withinGRE.nii.gz'];
    magn_fn         = [gre_basename '_part-mag_MEGRE_space-withinGRE.nii.gz'];
    mask_fn         = [gre_basename '_mask_MEGRE_space-withinGRE.nii.gz'];
    sepia_header_fn = [gre_basename '_header.mat'];
    output_prefix   = [gre_basename '_MEGRE_space-withinGRE'];

    % Input/Output filenames
    input{flip}(1).name   = fullfile(derivative_SEPIA_dir, phase_fn);
    input{flip}(2).name   = fullfile(derivative_SEPIA_dir, magn_fn);
    input{flip}(3).name   = '';
    input{flip}(4).name   = fullfile(derivative_SEPIA_dir, sepia_header_fn);
    output_basename{flip} = fullfile(seq_SEPIA_dir, output_prefix);
    mask_filename{flip}   = fullfile(derivative_SEPIA_dir, mask_fn);
end
if ~isdeployed
    qsubcellfun(@sepiaIO, input, output_basename, mask_filename, repmat({algorParam},size(input)), 'memreq', 4*1024^3, 'timreq', 60*60)
else
    cellfun(@sepiaIO, input, output_basename, mask_filename, repmat({algorParam},size(input)), 'UniformOutput', false)
end

algorParamR2star = struct();
algorParamR2star.general.isBET = 0;
algorParamR2star.general.isInvert = 0;
algorParamR2star.general.isRefineBrainMask = 0;

% R2* algorithm parameters
algorParamR2star.r2s.method = 'Trapezoidal';
algorParamR2star.r2s.s0mode = '1st echo';

for flip = 1:length(prot.flip)
    seq_SEPIA_dir   = fullfile(derivative_SEPIA_dir,prot.acq_str{flip});
    gre_basename    = [subj_label '_' prot.acq_str{flip} '_' run_label];
    magn_fn         = [gre_basename '_part-mag_MEGRE_space-withinGRE.nii.gz'];
    mask_fn         = [gre_basename '_mask_MEGRE_space-withinGRE.nii.gz'];
    sepia_header_fn = [gre_basename '_header.mat'];
    output_prefix   = [gre_basename '_MEGRE_space-withinGRE'];

    input{flip}(1).name   = fullfile(derivative_SEPIA_dir, phase_fn);
    input{flip}(2).name   = fullfile(derivative_SEPIA_dir, magn_fn);
    input{flip}(3).name   = '';
    input{flip}(4).name   = fullfile(derivative_SEPIA_dir, sepia_header_fn);
    output_basename{flip} = fullfile(seq_SEPIA_dir, output_prefix);
    mask_filename{flip}   = fullfile(derivative_SEPIA_dir, mask_fn);
end
if ~isdeployed
    qsubcellfun(@sepiaIO, input, output_basename, mask_filename, repmat({algorParamR2star},size(input)), 'memreq', 4*1024^3, 'timreq', 60*60)
else
    cellfun(@sepiaIO, input, output_basename, mask_filename, repmat({algorParamR2star},size(input)), 'UniformOutput', false)
end
