function ProcessingPipelineModular(prot, subj_label, run_label, bids_ses_dir, derivative_FSL_dir, derivative_SEPIA_dir, derivative_MRI_SYNTHSEG_dir)

% ProcessingPipeline

% 0 move data to derivatives folder derivatives/sub-x002/Prot1/Flip1/
% magnitude and phase
% 1 mcflirt first echo of each acquisition magnitude data
% 2 convert each echo  of each protocol to real and imaginary
% 2 could be skipped in the assumption of high SNR
% 3 apply motion correction
% 3b one could compute first a field map for each measurment and remove the main field, but it does not seem necessary

% 6 Coregister B1 map to S0

% Preparing data

for count_echo = 1:length(prot.echo)
    prot.echo_str{count_echo} = ['echo-' num2str(count_echo)];
end
for count_flip = 1:length(prot.flip)
    prot.flip_str{count_flip} = ['FA' num2str(prot.flip(count_flip))];
end
for count_flip = 1:length(prot.flip)
    prot.acq_str{count_flip} = [prot.rec 'FA' num2str(prot.flip(count_flip))];
end

bids_anat_dir = [bids_ses_dir '/anat/'];

[sts, msg] = mkdir(derivative_FSL_dir);
[sts, msg] = mkdir(derivative_SEPIA_dir);

for count_flip = 1:length(prot.flip)

    Dataset{count_flip}.rec    = prot.rec;
    Dataset{count_flip}.flip   = ['flip-' num2str(prot.flip(count_flip))];
    Dataset{count_flip}.outdir = fullfile(derivative_FSL_dir, prot.acq_str{count_flip});

    [sts, msg] = mkdir(Dataset{count_flip}.outdir);
    gre_basename               = [subj_label '_' prot.acq_str{count_flip} '_' run_label];
    list_json{count_flip}      = dir(fullfile(bids_anat_dir, [gre_basename '*_part-mag_MEGRE.json']));
    list_mag{count_flip}       = dir(fullfile(bids_anat_dir, [gre_basename '*_part-mag_MEGRE.nii.gz']));
    list_phase{count_flip}     = dir(fullfile(bids_anat_dir, [gre_basename '*_part-phase_MEGRE.nii.gz']));

    for count_echo = 1:length(prot.echo)
        % copies nifti and json phase
        run_command(['cp ' fullfile(bids_anat_dir, list_json{count_flip}(count_echo).name)  ' ' fullfile(Dataset{count_flip}.outdir, list_json{count_flip}(count_echo).name)], true);
        run_command(['cp ' fullfile(bids_anat_dir, list_mag{count_flip}(count_echo).name)   ' ' fullfile(Dataset{count_flip}.outdir, list_mag{count_flip}(count_echo).name)], true);
        run_command(['cp ' fullfile(bids_anat_dir, list_phase{count_flip}(count_echo).name) ' ' fullfile(Dataset{count_flip}.outdir, list_phase{count_flip}(count_echo).name)], true);
    end
    save_sepia_header(Dataset{count_flip}.outdir, [], fullfile(derivative_SEPIA_dir, gre_basename))
end

%%
Dataset{count_flip+1}.outdir = fullfile(derivative_FSL_dir, 'B1map');
[sts, msg] = mkdir(Dataset{count_flip+1}.outdir);

% Creating the reference volumes by doing a fast T1 map

% Extracting TR
json_str = fileread(fullfile(list_json{1}(1).folder, list_json{1}(1).name));
data     = jsondecode(json_str);
TR       = data.RepetitionTime;
S0       = load_untouch_nii(fullfile(list_mag{1}(1).folder, list_mag{1}(1).name));
dims     = size (S0.img);
tempS0   = zeros([dims length(prot.flip)]);

% Computes a T1 maps and S0 map
for count_flip = 1:length(prot.flip)
    gre_basename = [subj_label '_' prot.acq_str{count_flip} '_' run_label '_echo-1_part-mag_MEGRE.nii.gz'];
    S0  = load_untouch_nii(fullfile(list_mag{count_flip}(1).folder, gre_basename));
    tempS0 (:,:,:,count_flip) = S0.img;
    S0_target{count_flip} = ((fullfile(Dataset{count_flip}.outdir, 'EstimateFromT1mappping')));
end
[T1, M0] = despot1_mapping(tempS0, prot.flip, TR, []);

GRESignal = @(FlipAngle, TRep, T1) sind(FlipAngle) .* (1-exp(-TRep./T1)) ./ (1-(exp(-TRep./T1)) .* cosd(FlipAngle));

% Saves a T1 weighted images
for count_flip = 1:length(prot.flip)
    temp = M0 .* GRESignal(prot.flip(count_flip), TR, T1);
    temp(isinf(temp)) = 0;
    temp(isnan(temp)) = 0;
    save_nii_quick(S0, temp, S0_target{count_flip});
end

S0_target{count_flip+1} = fullfile(Dataset{count_flip+1}.outdir, 'EstimateFromT1mappping');
save_nii_quick(S0, M0, S0_target{count_flip+1});
clear tempS0

%% Making a fast M0 and T1 mapping

%% 3 Motion correction of all the first echo the similar volume

%  RefProt = fullfile(Dataset{RefFlip_index}.outdir,list_mag{RefFlip_index}(1).name);

for count_flip = 1:length(prot.flip)
    RefProt = S0_target{count_flip};
    gre_basename = [subj_label '_' prot.acq_str{count_flip} '_' run_label '_echo-1_part-mag_MEGRE.nii.gz'];
    S0 = fullfile(Dataset{count_flip}.outdir, gre_basename);

    extension = fullfile(Dataset{count_flip}.outdir, [prot.flip_str{count_flip} 'TransformToProtocolSpace.mat']);
    a = findstr(gre_basename, '.');

    S0_output = fullfile(Dataset{count_flip}.outdir, [gre_basename(1:(a(1)-1)) 'ProtocaolSpace.nii.gz']);
    run_command(['flirt -interp sinc -dof 6 -in ' S0 ' -ref ' RefProt ' -out ' S0_output ' -omat ' extension]);
end

RefProt       = S0_target{count_flip+1};
AnatB1        = fullfile(bids_ses_dir, 'fmap', [subj_label '_acq-anat_' run_label '_TB1TFL.nii.gz']);
B1map         = fullfile(bids_ses_dir, 'fmap', [subj_label '_acq-famp_' run_label '_TB1TFL.nii.gz']);
AnatB1_target = fullfile(derivative_FSL_dir, [subj_label '_acq-anat_' run_label '_TB1TFLProtocolSpace.nii.gz']);
B1map_target  = fullfile(derivative_FSL_dir, [subj_label '_acq-famp_' run_label '_TB1TFLProtocolSpace.nii.gz']);
extension     = fullfile(Dataset{count_flip+1}.outdir, 'TransformToProtocolSpace.mat');

% TODO: Fix "No image files match: ./sub-004/fmap/sub-004_acq-anat_run-1_TB1TFL"
run_command(['flirt -cost mutualinfo -in ' AnatB1 ' -ref ' RefProt ' -out ' AnatB1_target ' -omat ' extension]);
run_command(['flirt -in ' B1map ' -ref ' RefProt ' -applyxfm -init ' extension ' -out ' B1map_target]);


%% Applies co-registration to all the separate datasets

clear input
count = 0;
for count_flip = 1:length(prot.flip)
    RefProt = S0_target{count_flip};

    transform2common       = fullfile(Dataset{count_flip}.outdir, [prot.flip_str{count_flip} 'TransformToProtocolSpace.mat']);
    gre_basename           = [subj_label '_' prot.acq_str{count_flip} '_' run_label];
    list_mag{count_flip}   = dir(fullfile(Dataset{count_flip}.outdir, [gre_basename '*_part-mag_MEGRE.nii.gz']));
    list_phase{count_flip} = dir(fullfile(Dataset{count_flip}.outdir, [gre_basename '*_part-phase_MEGRE.nii.gz']));

    for count_meas = 1:length(list_phase{count_flip})
        count = count + 1;
        input{count}.MagName   = fullfile(list_mag{count_flip}(count_meas).folder, list_mag{count_flip}(count_meas).name);
        input{count}.PhaseName = fullfile(list_phase{count_flip}(count_meas).folder, list_phase{count_flip}(count_meas).name);
        input{count}.RefVol    = RefProt;
        input{count}.MotionMat = transform2common;
        
        temp = list_phase{count_flip}(count_meas).name;
        a = findstr(temp, '.');
        list_phase{count_flip}(count_meas).nameout = [temp(1:(a(1)-1)) 'ProtocolSpace' temp((a(1)):end)];
        input{count}.PhaseNameOutput = fullfile(derivative_FSL_dir, list_phase{count_flip}(count_meas).nameout);

        temp = list_mag{count_flip}(count_meas).name;
        a = findstr(temp, '.');
        list_mag{count_flip}(count_meas).nameout = [temp(1:(a(1)-1)) 'ProtocolSpace' temp((a(1)):end)];
        input{count}.MagNameOutput = fullfile(derivative_FSL_dir, list_mag{count_flip}(count_meas).nameout);
        % applyxfm4D_MagnPhase( input{count})
    end
end

if ~isdeployed
    res = qsubcellfun(@applyxfm4D_MagnPhase, input, 'memreq', 3*1024^3, 'timreq', 10*60, 'stack', 5);
else
    res = cellfun(@applyxfm4D_MagnPhase, input, 'UniformOutput', false);
end

%% Create a SEPIA folder with 4D data and a header file
for count_flip = 1:length(prot.flip)
    gre_mag      = [subj_label '_' prot.acq_str{count_flip} '_' run_label '_part-mag_MEGRE_space-withinGRE.nii.gz '];
    gre_phase    = [subj_label '_' prot.acq_str{count_flip} '_' run_label '_part-phase_MEGRE_space-withinGRE.nii.gz '];
    comand_mag   = ['fslmerge -t ' fullfile(derivative_SEPIA_dir, gre_mag)];
    comand_phase = ['fslmerge -t ' fullfile(derivative_SEPIA_dir, gre_phase)];
    for count_meas = 1:length(list_phase{count_flip})
        gre_basename = [subj_label '_' prot.acq_str{count_flip} '_' run_label '_echo-' num2str(count_meas) '_part-'];
        gre_mag      = [gre_basename 'mag_MEGREProtocolSpace.nii.gz '];
        gre_phase    = [gre_basename 'phase_MEGREProtocolSpace.nii.gz '];
        comand_mag   = [comand_mag   ' ' fullfile(derivative_FSL_dir, gre_mag)];
        comand_phase = [comand_phase ' ' fullfile(derivative_FSL_dir, gre_phase)];
    end
    run_command(comand_mag);
    run_command(comand_phase);
end

%% Get a brain mask for all datasets based on MRI syntseg
% 
% mri_synthseg --i ${converted_b1_dir}${in_vol} --o  ${ANTs_b1_b12gre_dir}${in_vol}_brain --cpu --robust

for count_flip = 1:length(prot.flip)
    gre_basename   = [subj_label '_' prot.acq_str{count_flip} '_' run_label '_echo-1_part-'];
    gre            = [gre_basename 'mag_MEGREProtocolSpace.nii.gz '];
    gre_mag        = [gre_basename 'mag_MEGREProtocolSpace_1mm.nii.gz '];
    gre_seg        = [gre_basename 'mag_MEGREProtocolSpace_1mmseg.nii.gz '];
    command{count_flip} = ['mri_synthseg --i ' fullfile(derivative_FSL_dir, gre) ' --o ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_seg) ' --cpu --robust --resample ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mag)];
end
fprintf('$ %s\n$ [..] %d times more\n', command{1}, length(command) - 1);
if ~isdeployed
    [status, output] = qsubcellfun(@run_command, command, 'memreq', 15*1024^3, 'timreq', 20*60);
else
    [status, output] = cellfun(@run_command, command, 'UniformOutput', false);
end

for count_flip = 1:length(prot.flip)
    gre_basename   = [subj_label '_' prot.acq_str{count_flip} '_' run_label '_echo-1_part-'];
    gre            = [gre_basename 'mag_MEGREProtocolSpace.nii.gz '];
    gre_mag        = [gre_basename 'mag_MEGREProtocolSpace_1mm.nii.gz '];
    gre_seg        = [gre_basename 'mag_MEGREProtocolSpace_1mmseg.nii.gz '];
    gre_mask       = [gre_basename 'mag_MEGREProtocolSpace_1mmmask.nii.gz '];
    gre_mag_ori    = [gre_basename 'mag_MEGREProtocolSpace_OriginalResampled.nii.gz '];
    gre_mask_ori   = [gre_basename 'mag_MEGREProtocolSpace_mask.nii.gz '];
    gre_mask_sepia = [subj_label '_' prot.acq_str{count_flip} '_' run_label '_mask_MEGRE_space-withinGRE.nii.gz '];
    % Check for errors
    if status{count_flip} ~= 0
        error('Command failed with status %d\nCommand:\n%s\nOutput:\n%s', status{count_flip}, command{count_flip}, output{count_flip});
    elseif ~isempty(output{count_flip})
        fprintf('%s\n', output{count_flip});
    end
    if ~isfile(strtrim(fullfile(derivative_MRI_SYNTHSEG_dir, gre_seg)))
        error(['No output file: ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_seg)])
    end

    % From here onwards it is just to get the mask on the SEPIA folder
    run_command(['fslmaths ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_seg) ' -bin ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mask)]);

    extension = fullfile(derivative_MRI_SYNTHSEG_dir, 'Transform');
    run_command(['flirt -in ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mag) ' -ref ' fullfile(derivative_FSL_dir, gre) ' -out ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mag_ori) ' -omat ' extension '.mat']);
    run_command(['flirt -interp nearestneighbour -in ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mask) ' -ref ' fullfile(derivative_FSL_dir, gre) ' -applyxfm -init ' extension '.mat -out ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mask_ori)]);
    run_command(['rm ' extension '.mat'], true);
    run_command(['cp ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mask_ori) fullfile(derivative_SEPIA_dir, gre_mask_sepia)], true);
    run_command(['rm ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mag_ori)], true);
    run_command(['rm ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mask_ori)], true);
end
