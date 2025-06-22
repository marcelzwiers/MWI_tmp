function ProcessingPipelineModular(subj_label, run_label, bids_ses_dir, derivative_FSL_dir, derivative_SEPIA_dir, derivative_MRI_SYNTHSEG_dir)

% ProcessingPipeline
debugmode = 1;
% 0 move data to derivatives folder derivatives/sub-x002/Prot1/Flip1/
% magnitude and phase
% 1 mcflirt first echo of each acquisition magnitude data
% 2 convert each echo  of each protocol to real and imaginary
% 2 could be skipped in the assumption of high SNR
% 3 apply motion correction
% 3b one could compute first a field map for each measurment and remove the main field, but it does not seem necessary

% 6 Coregister B1 map to S0

skipDataPreparation = 0;
% preparing data

count_prot = 1;
for count_echo = 1:length(prot{count_prot}.echo)
    prot{count_prot}.echo_str{count_echo} = ['echo-' num2str(count_echo)];
end
for count_flip = 1:length(prot{count_prot}.flip)
    prot{count_prot}.flip_str{count_flip} = ['FA' num2str(prot{count_prot}.flip(count_flip))];
end
for count_flip = 1:length(prot{count_prot}.flip)
    prot{count_prot}.acq_str{count_flip} = [prot{count_prot}.rec 'FA' num2str(prot{count_prot}.flip(count_flip))];
end

bids_anat_dir = [bids_ses_dir '/anat/'];

[sts, msg] = mkdir(derivative_FSL_dir);
[sts, msg] = mkdir(derivative_SEPIA_dir);
% count = 0;
% for count_prot = [1 4]
count_prot = 1;
for count_flip = 1:length(prot{1}.flip)
    % count = count + 1;

    Dataset{count_flip}.rec    = prot{count_prot}.rec;
    Dataset{count_flip}.flip   = ['flip-' num2str(prot{count_prot}.flip(count_flip))];
    Dataset{count_flip}.outdir = fullfile(derivative_FSL_dir, prot{count_prot}.acq_str{count_flip});

    [sts, msg] = mkdir(Dataset{count_flip}.outdir);
    gre_basename               = [subj_label '_' prot{count_prot}.acq_str{count_flip} '_' run_label];
    list_json{count_flip}      = dir(fullfile(bids_anat_dir, [gre_basename '*_part-mag_MEGRE.json']));
    list_mag{count_flip}       = dir(fullfile(bids_anat_dir, [gre_basename '*_part-mag_MEGRE.nii.gz']));
    list_phase{count_flip}     = dir(fullfile(bids_anat_dir, [gre_basename '*_part-phase_MEGRE.nii.gz']));

    if skipDataPreparation == 0
        for count_echo = 1:length(prot{count_prot}.echo)
            % copies nifti and json phase
            unix(['cp ' fullfile(bids_anat_dir, list_json{count_flip}(count_echo).name)  ' ' fullfile(Dataset{count_flip}.outdir, list_json{count_flip}(count_echo).name)]);
            unix(['cp ' fullfile(bids_anat_dir, list_mag{count_flip}(count_echo).name)   ' ' fullfile(Dataset{count_flip}.outdir, list_mag{count_flip}(count_echo).name)]);
            unix(['cp ' fullfile(bids_anat_dir, list_phase{count_flip}(count_echo).name) ' ' fullfile(Dataset{count_flip}.outdir, list_phase{count_flip}(count_echo).name)]);
        end
    end
    save_sepia_header(Dataset{count_flip}.outdir, [], fullfile(derivative_SEPIA_dir, gre_basename))
end

%%
Dataset{count_flip+1}.outdir = fullfile(derivative_FSL_dir, 'B1map');
[sts, msg] = mkdir(Dataset{count_flip+1}.outdir);

% Creating the reference volumes by doing a fast T1 map

% extracting TR
json_str = fileread(fullfile(list_json{1}(1).folder,list_json{1}(1).name));
data = jsondecode(json_str);
TR   = data.RepetitionTime;
S0   = load_untouch_nii(fullfile(list_mag{1}(1).folder,list_mag{1}(1).name));
dims = size (S0.img);
tempS0 = zeros([dims length(prot{count_prot}.flip)]);

% computes a T1 maps and S0 map
for count_flip = 1:length(prot{count_prot}.flip)
    gre_basename = [subj_label '_' prot{count_prot}.acq_str{count_flip} '_' run_label '_echo-1_part-mag_MEGRE.nii.gz'];
    S0  = load_untouch_nii(fullfile(list_mag{count_flip}(1).folder, gre_basename));
    tempS0 (:,:,:,count_flip) = S0.img;
    S0_target{count_flip} = ((fullfile(Dataset{count_flip}.outdir, 'EstimateFromT1mappping')));
end
[T1, M0] = despot1_mapping(tempS0, prot{count_prot}.flip, TR, []);

GRESignal = @(FlipAngle, TRep, T1) sind(FlipAngle) .* (1-exp(-TRep./T1)) ./ (1-(exp(-TRep./T1)) .* cosd(FlipAngle));

% saves a T1 weighted images
for count_flip = 1:length(prot{count_prot}.flip)
    temp = M0 .* GRESignal(prot{count_prot}.flip(count_flip), TR, T1);
    temp(isinf(temp)) = 0;
    temp(isnan(temp)) = 0;
    save_nii_quick(S0, temp, S0_target{count_flip});
end

S0_target{count_flip+1} = fullfile(Dataset{count_flip+1}.outdir, 'EstimateFromT1mappping');
save_nii_quick(S0, M0, S0_target{count_flip+1});
clear tempS0

%% making a fast M0 and T1 mapping

%% 3 Motion correction of all the first echo the similar volume

%  RefProt = fullfile(Dataset{RefFlip_index}.outdir,list_mag{RefFlip_index}(1).name);

for count_flip = 1:length(prot{1}.flip)
    RefProt = S0_target{count_flip};
    gre_basename = [subj_label '_' prot{count_prot}.acq_str{count_flip} '_' run_label '_echo-1_part-mag_MEGRE.nii.gz'];
    S0 = fullfile(Dataset{count_flip}.outdir, gre_basename);

    extension = fullfile(Dataset{count_flip}.outdir, [prot{count_prot}.flip_str{count_flip} 'TransformToProtocolSpace.mat']);
    a = findstr(gre_basename, '.');

    S0_output = fullfile(Dataset{count_flip}.outdir, [gre_basename(1:(a(1)-1)) 'ProtocaolSpace.nii.gz']);
    unix(['flirt -interp sinc -dof 6 -in ' S0 ' -ref ' RefProt ' -out ' S0_output ' -omat ' extension]);
end

RefProt       = S0_target{count_flip+1};
AnatB1        = fullfile(derivative_SIEMENS_dir, 'fmap', [subj_label '_acq-anat_' run_label '_TB1TFL.nii.gz']);
B1map         = fullfile(derivative_SIEMENS_dir, 'fmap', [subj_label '_acq-famp_' run_label '_TB1TFL.nii.gz']);
AnatB1_target = fullfile(derivative_FSL_dir, [subj_label '_acq-anat_' run_label '_TB1TFLProtocolSpace.nii.gz']);
B1map_target  = fullfile(derivative_FSL_dir, [subj_label '_acq-famp_' run_label '_TB1TFLProtocolSpace.nii.gz']);
extension     = fullfile(Dataset{count_flip+1}.outdir, 'TransformToProtocolSpace.mat');

unix(['flirt -cost mutualinfo -in ' AnatB1 ' -ref ' RefProt ' -out ' AnatB1_target ' -omat ' extension]);
unix(['flirt -in ' B1map ' -ref ' RefProt ' -applyxfm -init ' extension ' -out ' B1map_target]);


%% applies co-registration to all the separate datasets
clear input
count = 0;
for count_prot = [1]
    for count_flip = 1:length(prot{1}.flip)
        RefProt = S0_target{count_flip};
 
        transform2common       = fullfile(Dataset{count_flip}.outdir, [prot{count_prot}.flip_str{count_flip} 'TransformToProtocolSpace.mat']);
        gre_basename           = [subj_label '_' prot{count_prot}.acq_str{count_flip} '_' run_label];
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
end
res = qsubcellfun(@applyxfm4D_MagnPhase, input, 'memreq', 3 * 1024^3, 'timreq', 600);


%% Create a SEPIA folder with 4D data and a header file
count_prot = 1;
for count_flip = 1:length(prot{1}.flip)
    gre_mag      = [subj_label '_' prot{count_prot}.acq_str{count_flip} '_' run_label '_part-mag_MEGRE_space-withinGRE.nii.gz '];
    gre_phase    = [subj_label '_' prot{count_prot}.acq_str{count_flip} '_' run_label '_part-phase_MEGRE_space-withinGRE.nii.gz '];
    comand_mag   = ['fslmerge -t ' fullfile(derivative_SEPIA_dir, gre_mag)];
    comand_phase = ['fslmerge -t ' fullfile(derivative_SEPIA_dir, gre_phase)];
    for count_meas = 1:length(list_phase{count_flip})
        gre_basename = [subj_label '_' prot{count_prot}.acq_str{count_flip} '_' run_label '_echo-' num2str(count_meas) '_part-'];
        gre_mag      = [gre_basename 'mag_MEGREProtocolSpace.nii.gz '];
        gre_phase    = [gre_basename 'phase_MEGREProtocolSpace.nii.gz '];
        comand_mag   = [comand_mag   ' ' fullfile(derivative_FSL_dir, gre_mag)];
        comand_phase = [comand_phase ' ' fullfile(derivative_FSL_dir, gre_phase)];
    end
    unix(comand_mag);
    unix(comand_phase);
end

%% get a brain mask for all datasets based on MRI syntseg
% 
% mri_synthseg --i ${converted_b1_dir}${in_vol} --o  ${ANTs_b1_b12gre_dir}${in_vol}_brain --cpu --robust

count_prot = 1;
for count_flip = 1:length(prot{1}.flip)
    gre_basename    = [subj_label '_' prot{count_prot}.acq_str{count_flip} '_' run_label '_echo-1_part-'];
    gre             = [gre_basename 'mag_MEGREProtocolSpace.nii.gz '];
    gre_mag         = [gre_basename 'mag_MEGREProtocolSpace_1mm.nii.gz '];
    gre_seg         = [gre_basename 'mag_MEGREProtocolSpace_1mmseg.nii.gz '];
    gre_mask        = [gre_basename 'mag_MEGREProtocolSpace_1mmmask.nii.gz '];
    gre_mag_ori     = [gre_basename 'mag_MEGREProtocolSpace_OriginalResampled.nii.gz '];
    gre_mask_ori    = [gre_basename 'mag_MEGREProtocolSpace_mask.nii.gz '];
    gre_mask_sepia  = [subj_label '_' prot{count_prot}.acq_str{count_flip} '_' run_label '_mask_MEGRE_space-withinGRE.nii.gz '];

    command    = ['mri_synthseg --i ' fullfile(derivative_FSL_dir, gre) ' --o ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_seg) ' --cpu --robust --resample ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mag)];
    disp(['Running: ' command])
    [sts, out] = unix(command);
    if sts ~= 0
        error(['Failed to execute: ' command '\nError %i: ' out], sts)
    elseif ~isfile(strtrim(fullfile(derivative_MRI_SYNTHSEG_dir, gre_seg)))
        error(['No output file: ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_seg)])
    end

    % from here onwards it is just to get the mask on the SEPIA folder
    unix(['fslmaths ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_seg) ' -bin ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mask)]);

    extension = fullfile(derivative_MRI_SYNTHSEG_dir, 'Transform');
    unix(['flirt -in ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mag) ' -ref ' fullfile(derivative_FSL_dir, gre) ' -out ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mag_ori) ' -omat ' extension '.mat']);
    unix(['flirt -interp nearestneighbour -in ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mask) ' -ref ' fullfile(derivative_FSL_dir, gre) ' -applyxfm -init ' extension '.mat -out ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mask_ori)]);
    unix(['rm ' extension '.mat']);
    unix(['cp ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mask_ori) fullfile(derivative_SEPIA_dir, gre_mask_sepia)]);
    unix(['rm ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mag_ori)]);
    unix(['rm ' fullfile(derivative_MRI_SYNTHSEG_dir, gre_mask_ori)]);
end
