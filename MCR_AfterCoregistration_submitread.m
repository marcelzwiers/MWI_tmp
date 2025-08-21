function MCR_AfterCoregistration_submitread(input,output,task)
% MCR_AfterCoregistration_submitread(input,output)
% input.acq_str{countflip} %input acq_str of MGRE data
%
% input.derivative_SEPIA_dir % directory where the following datasets are stored
%     magn_fn         = [gre_basename '_part-mag_MEGRE_space-withinGRE.nii.gz'];
%     phase_fn        = [gre_basename '_MEGRE_space-withinGRE_part-phase_unwrapped.nii.gz'];
%     totalField_fn   = [gre_basename '_MEGRE_space-withinGRE_fieldmap.nii.gz'];
%     mask_fn         = [gre_basename '_MEGRE_space-withinGRE_mask_localfield.nii.gz'];
% input.derivative_FSL_dir % directory where B1 map saved in the MEGRE space is stored
% input.B1scaleFactor scaling factor for B1 map it is not compulsory, if not present uses json file field
% input.subj_label
% input.run_label
% input.MRvendor "Siemens"(deafault) or "Philips" selects the rf spoiling phase
% input.CorrectionFactorVFA  scale factor for VFA acquisition (should be some length as number of flip angles)
% input.Configfile % config file for MCR-MWI fitting parameters
% task.Submit_Job               % default = 1
% task.ReSubmit_MissingJobs     % default = 0
% task.Read_JobResults          % default = 0
if ~isfield(task, 'Submit_Job')
    task.Submit_Job = 1;
end
if ~isfield(task, 'ReSubmit_MissingJobs')
    task.ReSubmit_MissingJobs = 0;
end
if ~isfield(task, 'Read_JobResults')
    task.Read_JobResults = 0;
end
if ~isfield(output, 'MPPCAdenoise')
    output.MPPCAdenoise = 0;
end

% set up and run MCR-MWI fitting
if ~isfield(input,'Configfile')
    input.Configfile = [];
    algoParam = get_default_algoParam(input);
else
    % load the config file
    if exist(input.Configfile,'file')
        run(input.Configfile);
    else
        error('Config file not found');
    end
end

% example script to process MCR-MWI on a local computer with multi-threads
% The MCR-MWI processing scripts inside the subject folders were setted up
% slightly differently because of the parallelisation on HPC vs local PC

%% load GRE data

img             = [];
unwrappedPhase  = [];
mask            = [];
totalField      = [];
sepia_header    = [];
nFA             = (length(input.acq_str));
fa              = zeros(1,nFA);
for countflip = 1:length(input.acq_str)

    seq_SEPIA_dir = fullfile(input.derivative_SEPIA_dir, input.acq_str{countflip});

    % general GRE basename
    gre_basename    = [input.subj_label '_' input.acq_str{countflip} '_' input.run_label];

    % magnitude nifti image filename
    magn_fn         = [gre_basename '_part-mag_MEGRE_space-withinGRE.nii.gz'];
    phase_fn        = [gre_basename '_MEGRE_space-withinGRE_part-phase_unwrapped.nii.gz'];
    totalField_fn   = [gre_basename '_MEGRE_space-withinGRE_fieldmap.nii.gz'];
    mask_fn         = [gre_basename '_MEGRE_space-withinGRE_mask_localfield.nii.gz'];
    sepia_header_fn = [gre_basename '_header.mat'];
    nii             = load_untouch_nii(fullfile(input.derivative_SEPIA_dir, magn_fn));
    img             = cat(5, img,nii.img);
    nii             = load_untouch_nii(fullfile(seq_SEPIA_dir, phase_fn));
    unwrappedPhase  = cat(5, unwrappedPhase,nii.img);
    sepia_header{countflip}	= load(fullfile(input.derivative_SEPIA_dir, sepia_header_fn));
    mask            = cat(5, mask, load_nii_img_only(fullfile(seq_SEPIA_dir, mask_fn)));

    totalField      = cat(4, totalField, load_nii_img_only(fullfile(seq_SEPIA_dir, totalField_fn)));

    fa(countflip)   = sepia_header{countflip}.FA;
    tr              = sepia_header{countflip}.TR; % note that here there is an assumption that all protocols have the same TR

end

dims = size(img);
% B1 info
true_flip_angle_json    = [input.subj_label '_acq-famp_run-1_TB1TFL.json'];
true_flip_angle_fn      = [input.subj_label '_acq-famp_run-1_TB1TFLProtocolSpace.nii.gz'];
true_flip_angle         = load_nii_img_only(fullfile(input.derivative_FSL_dir, true_flip_angle_fn));

if isfield(input,'B1scaleFactor')
    b1                  = true_flip_angle / input.B1scaleFactor;
else
    b1_header           = jsondecode(fileread(fullfile(converted_b1_dir, true_flip_angle_json)));
    b1                  = true_flip_angle / 10 / b1_header.FlipAngle;
end

clear true_flip_angle

img = img .* exp(1i*unwrappedPhase);
if isfield(input,'CorrectionFactorVFA')
    nFA = length(input.CorrectionFactorVFA);
    for countflip= 1:nFA
        img(:,:,:,:,countflip) = img(:,:,:,:,countflip) * (input.CorrectionFactorVFA(countflip));
    end
end

mask = all(mask, 5) & all(~isnan(img), [4 5]);  % true where mask is true along 5th dim and no NaNs along 4th/5th dims

%% denoising part
if output.MPPCAdenoise == 1
    disp(['Denoising:' fullfile(input.derivative_SEPIA_dir, magn_fn)])
    if  or(task.ReSubmit_MissingJobs,task.Submit_Job)
        [denoised,~,~] = denoise(reshape(img,[dims(1:3) prod(dims(4:5))]), [5 5 5], mask);
        img = reshape(denoised,dims);
        clear denoised
    end
    PreProcessing = 'MPPCAdenoising';
else
    if output.MPPCAdenoise == 2
        disp(['Denoising:' fullfile(input.derivative_SEPIA_dir, magn_fn)])
        if  or(task.ReSubmit_MissingJobs,task.Submit_Job)
            [denoised,~,~] = denoise(reshape(img,[dims(1:3) prod(dims(4:5))]), [3 3 3], mask);
            img = reshape(denoised,dims);
            clear denoised
        end
        PreProcessing = 'MPPCAdenoising3';
    else
        disp('No Denoising')
        PreProcessing = 'no_preprocessing';
    end
end

%% obtain the initial estimation of the initial B1 phase
pini = squeeze(unwrappedPhase(:,:,:,1,:)) - 2*pi*totalField .* sepia_header{end}.TE(1);
Debug = 0;
if Debug==1
    % figure
    for k = 1:nFA
        figure(100)
        tiledlayout(ceil(sqrt(nFA)), ceil(sqrt(nFA)))
        nexttile(1)
        nexttile
        Orthoview2(pini(:,:,:,k),[],[],'tight')
        title(['Initial phase map for FA = ' num2str(fa(k))])
    end
end
pini = polyfit3D_NthOrder( double(mean(pini(:,:,:,1:(end-1)), 4)), mask, 6);

clear unwrappedPhase

% data normalisation
[scaleFactor, ~] = mwi_image_normalisation(img, mask);

% fixed parameters
kappa_mw    = 0.36; % jung
kappa_iew   = 0.86; % jung

imgParam.b0          = sepia_header{end}.B0;
imgParam.b0dir       = sepia_header{end}.B0_dir;
imgParam.rho_mw      = kappa_mw/kappa_iew;
imgParam.E           = 0.02;
imgParam.x_i         = -0.1;
imgParam.x_a         = -0.1;
imgParam.output_dir  = fullfile(output.derivative_MWI_dir, 'MCR', PreProcessing, ['using_' num2str(nFA) '_flipangle'], 'quadraticW');
gre_basename = [input.subj_label '_' output.acq_str '_' input.run_label];
imgParam.output_filename = [gre_basename '_MEGRE_MWI-MCR_' num2str(nFA) 'FA'];

% acquisition parameters
imgParam.te         = sepia_header{end}.TE;
imgParam.tr         = tr;
imgParam.fa         = fa;

if ~input.orthofit
    for slice = 1:dims(3)
        algoParamCell{slice} = algoParam;
        imgParamCell{slice} = imgParam;
        imgParamCell{slice}.output_filename = [gre_basename '_MEGRE_MWI-MCR_' num2str(nFA) 'FA_slice' num2str(slice)];
        if task.ReSubmit_MissingJobs || task.Submit_Job
            imgParamCell{slice}.img        = img(:,:,slice,:,:)/scaleFactor;
            imgParamCell{slice}.mask       = mask(:,:,slice);
            imgParamCell{slice}.fieldmap   = totalField(:,:,slice,:);
            imgParamCell{slice}.pini       = pini(:,:,slice,:);
            imgParamCell{slice}.b1map      = b1(:,:,slice);
            imgParamCell{slice}.autoSave   = 'true';
        end
    end
else
    slice = 1;
    ds_factor = 1;
    x = 1:ds_factor:dims(1);
    y = 1:ds_factor:dims(2);
    z = 1:ds_factor:dims(3);
    algoParamCell{slice} = algoParam;
    imgParamCell{slice}  = imgParam;
    for fa = 1:dims(5)
        for echo = 1:dims(4)
            imgParamCell{slice}.img(:,:,1,echo,fa) = OrthoSlice(img(x,y,z,echo,fa)/scaleFactor.*mask(x,y,z), [], 'tight');
        end
        imgParamCell{slice}.fieldmap(:,:,1,fa)   = OrthoSlice(totalField(x,y,z,fa).*mask(x,y,z), [], 'tight');
        imgParamCell{slice}.b1map(:,:,1,fa)      = OrthoSlice(b1(x,y,z).*mask(x,y,z), [] , 'tight');
    end
    imgParamCell{slice}.pini(:,:,1) = OrthoSlice(pini(x,y,z).*mask(x,y,z), [], 'tight');
    imgParamCell{slice}.mask        = OrthoSlice(mask(x,y,z).*mask(x,y,z), [], 'tight');
    imgParamCell{slice}.b1map       = OrthoSlice(b1(x,y,z).*mask(x,y,z), [], 'tight');
    imgParamCell{slice}.autoSave    = 'true';
    dims = size(imgParamCell{slice}.img);
    nii = make_nii(imgParamCell{slice}.img, [1 1 1]); % creates a fake nifti with the right shape
end

if task.Submit_Job
    for slice = 1:dims(3)
        if ~isdeployed
            jobid{slice} = qsubfeval(@mwi_3cx_2R1R2s_dimwi, algoParamCell{slice}, imgParamCell{slice}, 'memreq' , 1e10, 'timreq', 6*3600);
        else
            jobid{slice} = mwi_3cx_2R1R2s_dimwi(algoParamCell{slice}, imgParamCell{slice});
        end
    end
end

if task.ReSubmit_MissingJobs
    for slice = 1:dims(3)
        a = dir(fullfile(imgParamCell{slice}.output_dir, [imgParamCell{slice}.output_filename '.mat']));
        if isempty(a)
            disp(['Slice ' num2str(slice) ' was not present: Resubmitting job'])
            if ~isdeployed
                jobid{slice} = qsubfeval(@mwi_3cx_2R1R2s_dimwi, algoParamCell{slice}, imgParamCell{slice}, 'memreq', 1e10 , 'timreq', 6*3600);
            else
                jobid{slice} = mwi_3cx_2R1R2s_dimwi(algoParamCell{slice}, imgParamCell{slice});
            end
        end
    end
end
if and(task.ReSubmit_MissingJobs,task.Read_JobResults)
    for kz = 1:dims(3)
        a = dir(fullfile(imgParamCell{kz}.output_dir, [imgParamCell{kz}.output_filename '.mat']));
        while isempty(a)
            a = dir(fullfile(imgParamCell{kz}.output_dir, [imgParamCell{kz}.output_filename '.mat']));
            T = timer('TimerFcn', @(~,~) disp(['Slice ' num2str(kz)' is not yet ready']), 'StartDelay', 60);
            start(T)
            wait(T)
        end
        disp(['Slice ' num2str(kz) ' is ready'])
    end
end


MWF = @(x) x.S0_MW ./ (x.S0_MW+x.S0_EW+x.S0_IW);

%%
if task.Read_JobResults
    % read the job results

    try % normaly this is run inside the function that submitted the jobs,
        for kz = 1:dims(3)
            fitRes(kz) = qsubget(jobid{kz});
        end
        Alljobsread = 1;
    catch % but it could also have been run as a script and load the files
        clear fitRes
        Alljobsread = 1;
        for kz = 1:dims(3)
            try
                a =load(fullfile(imgParamCell{kz}.output_dir,[ imgParamCell{kz}.output_filename, '.mat'])) ;
                fitRes(kz)=a.fitRes;
            catch
                Alljobsread = 0;

                fitRes(kz)=fitRes(kz-1);
                disp([' Did not exist: ' fullfile(imgParamCell{kz}.output_dir, imgParamCell{kz}.output_filename)])
            end
        end
    end
    fitRes = ConcatStructure(fitRes,3);

    fitRes.MWF = MWF(fitRes);
    save_nii_quick(nii,fitRes.MWF*100,	          fullfile(imgParam.output_dir, [imgParam.output_filename '_MWFmap.nii.gz']));
    save_nii_quick(nii,fitRes.S0_MW,           	  fullfile(imgParam.output_dir, [imgParam.output_filename '_M0map-myelinwater.nii.gz']));
    save_nii_quick(nii,fitRes.S0_IW+fitRes.S0_EW, fullfile(imgParam.output_dir, [imgParam.output_filename '_M0map-freewater.nii.gz']));
    save_nii_quick(nii,fitRes.R2s_MW,             fullfile(imgParam.output_dir, [imgParam.output_filename '_R2starmap-myelinwater.nii.gz']));
    save_nii_quick(nii,fitRes.R2s_IW,             fullfile(imgParam.output_dir, [imgParam.output_filename '_R2starmap-intraaxonal.nii.gz']));
    save_nii_quick(nii,fitRes.T1_IEW,             fullfile(imgParam.output_dir, [imgParam.output_filename '_T1map-freewater.nii.gz']));
    save_nii_quick(nii,1./fitRes.T1_IEW,          fullfile(imgParam.output_dir, [imgParam.output_filename '_R1map-freewater.nii.gz']));
    save_nii_quick(nii,fitRes.kiewm,              fullfile(imgParam.output_dir, [imgParam.output_filename '_exchangerate-freewatertomyelinwater.nii.gz']));
    save_nii_quick(nii,fitRes.Freq_BKG,           fullfile(imgParam.output_dir, [imgParam.output_filename '_Frequencymap-background.nii.gz']));
    save_nii_quick(nii,fitRes.pini,               fullfile(imgParam.output_dir, [imgParam.output_filename '_Initialphase.nii.gz']));
    save_nii_quick(nii,fitRes.exitflag,           fullfile(imgParam.output_dir, [imgParam.output_filename '_exitflag.nii.gz']));
    save_nii_quick(nii,fitRes.iterations,         fullfile(imgParam.output_dir, [imgParam.output_filename '_numiterations.nii.gz']));
    save_nii_quick(nii,fitRes.mask_fitted,        fullfile(imgParam.output_dir, [imgParam.output_filename '_mask_fittedvoxel.nii.gz']));
else
    Alljobsread = 0;
end
if Alljobsread == 1
    for kz = 1:dims(3)
        delete(fullfile(imgParamCell{kz}.output_dir, imgParamCell{kz}.output_filename));
    end
end


function algoParam = get_default_algoParam(input)
    % %% set up and run MCR-MWI fitting
    % % setup algorithm parameters
    algoParam.isInvivo      = true;     % true for using initial guesses for in vivo imaging
    algoParam.isParallel    = false;     % true: using parfor parallel processing; false: no parfor
    algoParam.DEBUG         = false;    % true: debug mode to display some info
    algoParam.isNormData    = false;     % true: normalise the data by a global constant so that absolute tolerance will be used for fitting, here we use false because we did normalisation outside the function, see below
    % fitting option
    algoParam.maxIter       = 200;      % maximum number of iterations
    algoParam.fcnTol        = 1e-4;     % fitting tolerance, this valuse is for normalised data
    algoParam.stepTol       = 1e-4;     % step tolerance, this valuse is for normalised data
    % residual option
    algoParam.numMagn       = 0;        % 0: complex fitting
    algoParam.isWeighted    = true;     % true: using magnitude signal weighting
    algoParam.weightMethod  = 'quadratic_1stEcho';  % '1stEcho': weighted by 1st echo magnitude; 'quadratic_1stEcho': weighted by squared 1st echo magnitude
    % T1 model
    algoParam.isExchange    = 1;        % BM model
    algoParam.isEPG         = 1;        % Using EPG-X for signal simulation
    algoParam.npulse        = 50;       % number of pulses to reach steady-state
    if ~isfield(input,'MRvendor')
        algoParam.rfphase       = 50;       % RF phase for EPG, degree
    else
        if strcmpi(input.MRvendor, 'siemens')
            algoParam.rfphase       = 50;       % RF phase for EPG, degree
        else
            algoParam.rfphase       = 150;      % phase increment for Philips
        end
    end

    algoParam.isT1mw        = false;    % true: fitting myelin water T1, false: use fixed value
    algoParam.T1mw          = 234e-3;   % define fixed myelin water T1 value
    % No DIMWI
    algoParam.DIMWI.isVic       = false;    % false: no extra DWI info for DIMWI
    algoParam.DIMWI.isR2sEW     = false;    % false: no extra DWI info for DIMWI
    algoParam.DIMWI.isFreqMW    = false;    % false: no extra DWI info for DIMWI
    algoParam.DIMWI.isFreqIW    = false;    % false: no extra DWI info for DIMWI
    % initial guess
    algoParam.advancedStarting = 'default';  % initial guesses for multi-comp S0


function mosaic = OrthoSlice(volume, xyz, showim)
    % mosaic = OrthoSlice(volume, xyz, showim)
    %
    % Extract orthogonal slices from a 3D volume.
    %
    % Inputs:
    %   volume - 3D volume
    %   xyz    - slice positions (default: center of volume)
    %   showim - display type: 'normal' (default) or 'tight'

    % Defaults
    if nargin < 3 || isempty(showim)
        showim = 'normal';
    end
    if strcmp(showim, 'tight')
        [x,y,z] = ind2sub(size(volume), find(volume));
        volume  = volume(min(x):max(x), min(y):max(y), min(z):max(z));
    end
    dims = size(volume);
    if nargin < 2 || isempty(xyz)
        xyz = round(dims/2);
    end

    mosaic = zeros([max(dims(2:3)) 2*dims(1) + dims(2)]);
    temp1  = zeros([size(mosaic,1), dims(2)]);
    temp2  = zeros([size(mosaic,1), dims(1)]);
    temp3  = zeros([size(mosaic,1), dims(1)]);
    if ismember(showim, {'normal','tight'})
        temp1a = permute(volume(xyz(1),:,:), [3,2,1]);
        temp2a = permute(volume(:,xyz(2),:), [3,1,2]);
        temp3a = permute(volume(:,:,xyz(3)), [2,1,3]);
    end

    temp1(round((size(temp1,1) - size(temp1a,1))/2) + (1:size(temp1a,1)), :) = temp1a;
    temp2(round((size(temp2,1) - size(temp2a,1))/2) + (1:size(temp2a,1)), :) = temp2a;
    temp3(round((size(temp3,1) - size(temp3a,1))/2) + (1:size(temp3a,1)), :) = temp3a;
    mosaic = cat(2, temp1, temp2, temp3);
