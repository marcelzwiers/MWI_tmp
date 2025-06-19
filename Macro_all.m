cd '/project/3055010.04/RunningProjects/MyelinWaterImaging/bidsSiemensVariantsNew/code'

preprocessing = 1;
SepiaPrep     = 1;
fittingMCR    = 1;
writingMCR    = 1;
subjcell{2}   = 'fl3dmcrG2x2';
subjcell{1}   = 'fl3dmcrG2x2new';

for subjn = 1:length(subjcell)
    fprintf('Processing: %s (%d/%d)\n', subjcell{subjn}, subjn, length(subjcell));

    subj         = subjcell{subjn};
    subj_label   = ['sub-' subj];
    subject_directory_master
    sub          = subj_label;
    ses          = 'ses-mri01';
    acqname      = 'GRE';
    run_label    = 'run-1';

    prot{1}.rec  = ['acq-' acqname];
    prot{1}.flip = [10 20 50 70];
    prot{1}.echo = 1:12;

    for k = 1:length(prot{1}.flip)
        acq_label{k} = [prot{1}.rec 'FA' num2str(prot{1}.flip(k))];
    end

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

    addpath('/home/common/matlab/sepia/sepia_1.2.2.5/');    % https://github.com/kschan0214/sepi
    sepia_addpath
    addpath(code_dir);                                      % 'code' directory
    addpath(genpath(fullfile(code_dir,'despot1')));         % https://github.com/kschan0214/despot1
    addpath(genpath(fullfile(code_dir,'EPG-X'))); 	        % KCL extended phase graphs
    addpath(genpath(fullfile(code_dir,'utils')));
    addpath(fullfile(code_dir,'MP-PCA-Denoising'));
    addpath('/home/common/matlab/fieldtrip/qsub')

    if preprocessing
        ProcessingPipelineModular
    end

    if SepiaPrep
        SEPIA_03_standard_pipeline
        script_SCR
    end

    % some solvers have similar names in various packages and we have to make sure the MWI is the one that gets used
    % rmpath(genpath(fullfile(code_dir,'despot1')));
    warning off
    rmpath(genpath('/home/common/matlab/sepia/sepia_1.2.2.5/'));
    rmpath(genpath(fullfile(code_dir,'mwi')));
    warning on
    addpath(genpath(fullfile(code_dir,'mwi')));

    if fittingMCR || writingMCR
        % Check that MWI toolbox is at the top of the path it has a solver that has conflicts with SEPIA and despot
        % MPPCAdenoise = 1; %this actually has a positive effect on maps
        if fittingMCR
            task.Submit_Job           = 1;
            task.ReSubmit_MissingJobs = 0;
            task.Read_JobResults      = 0; % ideally one would have this one as also 1, but that just takes too much time
        end
        if writingMCR
            task.Submit_Job           = 0;
            task.ReSubmit_MissingJobs = 0;
            task.Read_JobResults      = 0; % only do this if enough slices have successfully been processed
        end

        input                      = struct();
        input.derivative_SEPIA_dir = derivative_SEPIA_dir;
        input.derivative_FSL_dir   = derivative_FSL_dir;
        input.acq_str              = prot{1}.acq_str;
        input.B1scaleFactor        = 800;               % directory where json b1 information is present alternatively it can be the scaling factor of B1 field
        input.subj_label           = subj_label;
        input.run_label            = run_label;
        output.acq_str             = prot{1}.rec;
        output.derivative_MWI_dir  = derivative_MWI_dir; % main output directory

        output.MPPCAdenoise        = 0;
        func_MCR_AfterCoregistration_qsubfeval_submitread(input, output, task)

        output.MPPCAdenoise        = 1;
        func_MCR_AfterCoregistration_qsubfeval_submitread(input, output, task)
    end
end
disp("Finished processing")

return
%% LOAD DATA FOR FIGURES

Preprocessingcell{1} = 'MPPCAdenoising';
Preprocessingcell{2} = 'no_preprocessing';
MWI_dir = '/project/3055010.04/RunningProjects/MyelinWaterImaging/bidsSiemensVariantsNew/derivatives/MWI/';
SCR_dir = '/project/3055010.04/RunningProjects/MyelinWaterImaging/bidsSiemensVariantsNew/derivatives/SimultaneousR1R2star/';
MWIthreshold = 5;
for subjn = 1:length(subjcell)
    subj       = subjcell{subjn};
    subj_label = ['sub-' subj];
    for prepn = 1:length(Preprocessingcell)
        prep_label                     = ['MCR/' Preprocessingcell{prepn} '/using_4_flipangle/quadraticW'];
        MWI_data(subjn,prepn)          = load_untouch_nii(fullfile(MWI_dir,subj_label,prep_label,[subj_label,'_acq-GRE_run-1_MEGRE_MWI-MCR_4FA_MWFmap.nii.gz']));
        R1mapMCR_data(subjn,prepn)     = load_untouch_nii(fullfile(MWI_dir,subj_label,prep_label,[subj_label,'_acq-GRE_run-1_MEGRE_MWI-MCR_4FA_R1map-freewater.nii.gz']));
        ExchangeRate_data(subjn,prepn) = load_untouch_nii(fullfile(MWI_dir,subj_label,prep_label,[subj_label,'_acq-GRE_run-1_MEGRE_MWI-MCR_4FA_exchangerate-freewatertomyelinwater.nii.gz']));
        R2sIW_MCR_data(subjn,prepn)    = load_untouch_nii(fullfile(MWI_dir,subj_label,prep_label,[subj_label,'_acq-GRE_run-1_MEGRE_MWI-MCR_4FA_R2starmap-intraaxonal.nii.gz']));
        R2sMW_MCR_data(subjn,prepn)    = load_untouch_nii(fullfile(MWI_dir,subj_label,prep_label,[subj_label,'_acq-GRE_run-1_MEGRE_MWI-MCR_4FA_R2starmap-myelinwater.nii.gz']));
        M0_MCR_data(subjn,prepn)       = load_untouch_nii(fullfile(MWI_dir,subj_label,prep_label,[subj_label,'_acq-GRE_run-1_MEGRE_MWI-MCR_4FA_M0map-freewater.nii.gz']));

        mask_data(subjn,prepn) = load_untouch_nii(fullfile(MWI_dir,subj_label,prep_label,[subj_label,'_acq-GRE_run-1_MEGRE_MWI-MCR_4FA_mask_fittedvoxel.nii.gz']));
        MWI_data(subjn,prepn).img((mask_data(subjn,prepn).img)==0)=0;
        R2sIW_MCR_data(subjn,prepn).img((mask_data(subjn,prepn).img)==0)=0;
        R2sMW_MCR_data(subjn,prepn).img((mask_data(subjn,prepn).img)==0)=0;
        ExchangeRate_data(subjn,prepn).img((mask_data(subjn,prepn).img)==0)=0;
        R2sMW_MCR_data(subjn,prepn).img((MWI_data(subjn,prepn).img)<MWIthreshold)=0;
        ExchangeRate_data(subjn,prepn).img((MWI_data(subjn,prepn).img)<MWIthreshold)=0;
        R1mapMCR_data(subjn,prepn).img((mask_data(subjn,prepn).img)==0)=0;

    end
        R1mapSCR_data(subjn)  = load_untouch_nii(fullfile(SCR_dir,subj_label,[subj_label,'_acq-GRE_run-1_MEGRE_space-withinGRE_R1map.nii.gz']));
        R2smapSCR_data(subjn) = load_untouch_nii(fullfile(SCR_dir,subj_label,[subj_label,'_acq-GRE_run-1_MEGRE_space-withinGRE_R2starmap.nii.gz']));
        M0mapSCR_data(subjn)  = load_untouch_nii(fullfile(SCR_dir,subj_label,[subj_label,'_acq-GRE_run-1_MEGRE_space-withinGRE_M0map.nii.gz']));
end

%% MAKE FIGURES

Preprocessingcell    = [];
Preprocessingcell{1} = 'MPPCAdenoising';

figureJ(1)
close(1)
figureJ(1)
for subjn = 1:length(subjcell)
    for prepn = 1:length(Preprocessingcell)
        nexttile()
        Orthoview2(MWI_data(subjn,prepn).img.*mask_data(subjn,prepn).img,[],[0 20],'tight');
        title(['Prococol - ' ,subjcell{subjn}, '  ', Preprocessingcell{prepn}],'Interpreter','none');
        colorbar('south', 'Color', [1 1 1]);
    end
end
set(gcf,'position',[1,1,810,1070])

figureJ(2)
close(2)
figureJ(2)
for subjn = 1:length(subjcell)
    for prepn = 1:length(Preprocessingcell)
        nexttile()
        Orthoview2(R1mapMCR_data(subjn,prepn).img.*mask_data(subjn,prepn).img,[],[0.2 1.600],'tight');
        title(['Prococol - ' ,subjcell{subjn}, '  ', Preprocessingcell{prepn}],'Interpreter','none');
        colorbar('south', 'Color', [1 1 1]);
    end
end
set(gcf,'position',[1,1,810,1070])

figureJ(3)
close(3)
figureJ(3)
for subjn = 1:length(subjcell)
    for prepn = 1:length(Preprocessingcell)
        nexttile()
        Orthoview2(ExchangeRate_data(subjn,prepn).img.*mask_data(subjn,prepn).img,[],[0 5],'tight');
        title(['Prococol - ' ,subjcell{subjn}, '  ', Preprocessingcell{prepn}],'Interpreter','none')
        colorbar('south', 'Color', [1 1 1]);
    end
end
set(gcf,'position',[1,1,810,1070])

figureJ(4)
close(4)
figureJ(4)
for subjn = 1:length(subjcell)
    for prepn = 1:length(Preprocessingcell)
        nexttile()
        Orthoview2(R2sMW_MCR_data(subjn,prepn).img.*mask_data(subjn,prepn).img,[],[50 150],'tight');
        title(['Prococol - ' ,subjcell{subjn}, '  ', Preprocessingcell{prepn}],'Interpreter','none')
        colorbar('south', 'Color', [1 1 1]);
    end
end
set(gcf,'position',[1,1,810,1070])

figureJ(5)
close(5)
figureJ(5)
for subjn = 1:length(subjcell)
    for prepn = 1:length(Preprocessingcell)
        nexttile()
        Orthoview2(R2sIW_MCR_data(subjn,prepn).img.*mask_data(subjn,prepn).img,[],[0 50],'tight');
        title(['Prococol - ' ,subjcell{subjn}, '  ', Preprocessingcell{prepn}],'Interpreter','none')
        colorbar('south', 'Color', [1 1 1]);
    end
end
set(gcf,'position',[1,1,810,1070])

figureJ(6)
close(6)
figureJ(6)
for subjn = 1:length(subjcell)
    for prepn = 1:length(Preprocessingcell)
        nexttile()
        Orthoview2(M0_MCR_data(subjn,prepn).img.*mask_data(subjn,prepn).img,[],[],'tight');
        title(['Prococol - ' ,subjcell{subjn}, '  ', Preprocessingcell{prepn}],'Interpreter','none')
        colorbar('south', 'Color', [1 1 1]);
    end
end
set(gcf,'position',[1,1,810,1070])


%%
figureJ(12)
close(12)
figureJ(12)
for subjn = 1:length(subjcell)
    % for prepn = 1:length(Preprocessingcell)
        nexttile()
        Orthoview2(R1mapSCR_data(subjn).img.*mask_data(subjn,prepn).img,[],[0 1600],'tight')
        title(['Prococol - ' ,subjcell{subjn}],'Interpreter','none')
        colorbar('south', 'Color', [1 1 1]);
    % end
end
set(gcf,'position',[1,1,810,1070])


figureJ(13)
close(13)
figureJ(13)
for subjn = 1:length(subjcell)
    % for prepn = 1:length(Preprocessingcell)
        nexttile()
        R2smapSCR_data(subjn).img(isnan(R2smapSCR_data(subjn).img))=0;
        R2smapSCR_data(subjn).img(isinf(R2smapSCR_data(subjn).img))=0;
        
        Orthoview2(R2smapSCR_data(subjn).img.*mask_data(subjn,prepn).img,[],[0 50],'tight')
        title(['Prococol - ' ,subjcell{subjn}, '  '],'Interpreter','none')
        colorbar('south', 'Color', [1 1 1]);
    % end
end
set(gcf,'position',[1,1,810,1070])

figureJ(14)
close(14)
figureJ(14)
for subjn = 1:length(subjcell)
    % for prepn = 1:length(Preprocessingcell)
        nexttile()
        Orthoview2(M0mapSCR_data(subjn).img.*mask_data(subjn,prepn).img,[],[],'tight')
        title(['Prococol - ' ,subjcell{subjn} ],'Interpreter','none')
    % end
end
set(gcf,'position',[1,1,810,1070])

%%
RawData

%% this in only for visalization purposes
R2smean(isnan(R2smean)) = 0;
MWF(isnan(MWF))         = 0;
Chimean(isnan(Chimean)) = 0;
Nsubplots               = 6;

widthheigth = [1./round(sqrt(Nsubplots)), 1./(ceil(Nsubplots/round(sqrt(Nsubplots))))];
count = 0;
for l = 1:(ceil(Nsubplots/round(sqrt(Nsubplots))))
    for k = 1:round(sqrt(Nsubplots))
        count = count + 1;
        pos{count} = [(k-1)*widthheigth(1), 1-l*widthheigth(2), widthheigth];
    end
end

% pos{1} = [0   0.5 0.5 0.5];
% pos{2} = [0.5 0.5 0.5 0.5];
% pos{3} = [0   0   0.5 0.5];
% pos{4} = [0.5 0   0.5 0.5];
position = [74 70 100];

figure(1);

subplot('position', pos{1})
Orthoview2(R1.*mask, [position], [0.2 2]*1000, 'tight');
title('R1 [KHz]'); colorbar('south', 'Color', [1 1 1]);

subplot('position', pos{2})
Orthoview2(M0.*mask, [position], [], 'tight');
title('M0 [au]'); colorbar('south', 'Color', [1 1 1]);

subplot('position', pos{3})
Orthoview2(R2smean.*mask, [position], [0.2 50], 'tight');
title('R2s [Hz]'); colorbar('south', 'Color', [1 1 1]);

subplot('position', pos{4})
Orthoview2(Chimean.*mask, [position], [-0.15 0.45], 'tight');
title('Chi [ppm]'); colorbar('south', 'Color', [1 1 1]);

subplot('position', pos{5})
Orthoview2(MWF.*mask, [position], [0 15], 'tight');
title('MWF [%]'); colorbar('south', 'Color', [1 1 1]);


Nslices  = 15;
index    = find(mask~=0);
[x,y,z ] = ind2sub(size(mask),index);
xrange   = min(x):max(x);
yrange   = min(y):max(y);
zrange   = round(min(z):(max(z)-min(z))/(Nslices+1):max(z));
zrange   = zrange(2:end-1);

figure(2)
set(gcf, 'Color', [1,1,1]);

subplot('position', pos{1})
imab(R1(xrange,yrange,zrange)/1000, [0.0 1.4]);
ylabel('R1 [Kz]'); colorbar('eastoutside', 'Color', [1 1 1]*0);

subplot('position', pos{2})
imab(M0(xrange, yrange, zrange), [6000 12000]);
ylabel('M0 [au]'); colorbar('eastoutside', 'Color', [1 1 1]*0);

subplot('position', pos{3})
imab(R2smean(xrange, yrange, zrange), [0.2 50]);
ylabel('R2s [Hz]'); colorbar('eastoutside','Color', [1 1 1]*0);

subplot('position', pos{4})
imab(Chimean(xrange,yrange,zrange),[-0.15 0.45]);
ylabel('Chi [ppm]'); colorbar('eastoutside','Color',[1 1 1]*0);

subplot('position', pos{5})
imab(MWF(xrange, yrange, zrange), [0 12]);
ylabel('MWF[%]'); colorbar('eastoutside', 'Color', [1 1 1]*0);

subplot('position', pos{6})
imab(1./T1_IEW(xrange, yrange, zrange), [0 1.2]);
ylabel('R1 myelin free MCR[Hz]'); colorbar('eastoutside', 'Color', [1 1 1]*0);

colormap("gray")
