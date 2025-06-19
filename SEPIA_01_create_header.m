%% Subject info and directories
subj	= '001';
run     = '1';
ses     = 'mri01';

subj_label  = ['sub-' subj];
sess_label  = ['ses-' ses];
sess_label  = [];
run_label   = ['run-' run];

% load path
% addpath('.');           	% 'code' directory
addpath(genpath('sepia-1.0.0/'));	% https://github.com/kschan0214/sepia/releases/tag/v1.0.0 
addpath(genpath('utils/'));	% https://github.com/kschan0214/sepia/releases/tag/v1.0.0 

% call diretcories
subject_directory_master
mkdir(derivative_SEPIA_dir)

subj_script_dir = fullfile(code_dir, subj_label, sess_label);

% GRE Protocols and flip angles
protocol    = { 'fl3d',};
flipAngle   = { 'FA10', 'FA20', 'FA50', 'FA70' };

%% create SEPIA header
for kp = protocol
    
     
    prot_SEPIA_dir  = fullfile(derivative_SEPIA_dir, kp{1});
    mkdir(prot_SEPIA_dir)
    
    for kf = flipAngle
        
        % general GRE basename
        seq             = [kp{1} kf{1}];
        acq_label       = ['acq-' seq];
        gre_basename    = [subj_label '_' sess_label '_' acq_label '_' run_label];
        gre_basename    = [subj_label '_' acq_label '_' run_label];
        sepia_header_fn = [gre_basename '_MEGRE_sepia_header.mat'];
        
        % magnitude nifti image filename
        magn_fn = [gre_basename '_echo-1_part-mag_MEGRE.nii.gz'];
        
        % all json filenames
        list        = dir(fullfile(converted_anat_dir, [gre_basename '*_part-mag_MEGRE.json']));
        json_list   = []; for kl = 1:length(list); json_list{kl} = fullfile(list(kl).folder, list(kl).name); end
        
        % fsl's flirt transformation matrix
        fsl_mat_fn	= [gre_basename '_echo-1_part-mag_MEGRE_to_reference_0GenericAffine.fsl'];
        
        % load nifti for header
        nii         = load_untouch_nii(fullfile(converted_anat_dir, magn_fn));
        % load json for header
        json_str    = fileread(json_list{end}); 
        data        = jsondecode(json_str);
        % load fsl transformation matrix
        fsl_mat     = GetFSLRotationMatrix(fullfile(prot_ANTs_dir, fsl_mat_fn));
        
        % get image size
        matrixSize  = nii.hdr.dime.dim(2:4);
        voxelSize   = nii.hdr.dime.pixdim(2:4);
        % get B0 direction
        B0_dir_native = get_B0_dir_from_nifti(nii);
        % rotate B0 direction after motion correction
        B0_dir      = fsl_mat * B0_dir_native(:);
        
        % get other relevant parameters
        FA          = data.FlipAngle;
        TR          = data.RepetitionTime;
        B0          = data.MagneticFieldStrength;
        CF          = data.ImagingFrequency * 1e6;
        TE          = readTE_dcm2niix_JSON(json_list);
        delta_TE    = TE(2) - TE(1);
        
        save(fullfile(prot_SEPIA_dir,sepia_header_fn),'voxelSize','matrixSize','CF','delta_TE',...
        'TE','B0_dir','B0_dir_native','B0','FA','TR');
        
    end
end
