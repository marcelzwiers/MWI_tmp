function CreateFiguresForReport(bidsdir,subj_label,acq_label)



% Preprocessingcell{1} = 'MPPCAdenoising';
Preprocessingcell{1} = 'no_preprocessing';
MWI_dir = fullfile(bidsdir,'derivatives/MWI');
SCR_dir = fullfile(bidsdir,'derivatives/SimultaneousR1R2star');
MWIthreshold = 5;
subjn = 1;
prepn =1;
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

        R1mapSCR_data(subjn)  = load_untouch_nii(fullfile(SCR_dir,subj_label,[subj_label,'_acq-GRE_run-1_MEGRE_space-withinGRE_R1map.nii.gz']));
        R2smapSCR_data(subjn) = load_untouch_nii(fullfile(SCR_dir,subj_label,[subj_label,'_acq-GRE_run-1_MEGRE_space-withinGRE_R2starmap.nii.gz']));
        M0mapSCR_data(subjn)  = load_untouch_nii(fullfile(SCR_dir,subj_label,[subj_label,'_acq-GRE_run-1_MEGRE_space-withinGRE_M0map.nii.gz']));


pos =[1,1,810,1070]
pos =[6    ,    -339,810,1070]

% figureJ(1)
% close(1)
% figureJ(1)
% for subjn = 1:length(subjcell)
%     for prepn = 1:length(Preprocessingcell)
%         nexttile()
%         Orthoview2(MWI_data(subjn,prepn).img.*mask_data(subjn,prepn).img,[],[0 20],'tight');
%         title(['MWI - Prococol - ' ,subjcell{subjn}, '  ', Preprocessingcell{prepn}],'Interpreter','none');
%         colorbar('south', 'Color', [1 1 1]);
%     end
% 
% end
% set(gcf,'position',pos)
% 
% figureJ(2)
% close(2)
% figureJ(2)
% for subjn = 1:length(subjcell)
%     for prepn = 1:length(Preprocessingcell)
%         nexttile()
%         Orthoview2(R1mapMCR_data(subjn,prepn).img.*mask_data(subjn,prepn).img,[],[0.2 1.600],'tight');
%         title(['R1map  - Prococol - ' ,subjcell{subjn}, '  ', Preprocessingcell{prepn}],'Interpreter','none');
%         colorbar('south', 'Color', [1 1 1]);
%     end
% end
% set(gcf,'position',pos)
% 
% figureJ(3)
% close(3)
% figureJ(3)
% for subjn = 1:length(subjcell)
%     for prepn = 1:length(Preprocessingcell)
%         nexttile()
%         Orthoview2(ExchangeRate_data(subjn,prepn).img.*mask_data(subjn,prepn).img,[],[0 5],'tight');
%         title(['ExchangeRate - Prococol - ' ,subjcell{subjn}, '  ', Preprocessingcell{prepn}],'Interpreter','none')
%         colorbar('south', 'Color', [1 1 1]);
%     end
% end
% set(gcf,'position',pos)
% 
% figureJ(4)
% close(4)
% figureJ(4)
% for subjn = 1:length(subjcell)
%     for prepn = 1:length(Preprocessingcell)
%         nexttile()
%         Orthoview2(R2sMW_MCR_data(subjn,prepn).img.*mask_data(subjn,prepn).img,[],[50 150],'tight');
%         title(['R2sMW - Prococol - ' ,subjcell{subjn}, '  ', Preprocessingcell{prepn}],'Interpreter','none')
%         colorbar('south', 'Color', [1 1 1]);
%     end
% end
% set(gcf,'position',pos)
% 
% figureJ(5)
% close(5)
% figureJ(5)
% for subjn = 1:length(subjcell)
%     for prepn = 1:length(Preprocessingcell)
%         nexttile()
%         Orthoview2(R2sIW_MCR_data(subjn,prepn).img.*mask_data(subjn,prepn).img,[],[0 50],'tight');
%         title(['R2sIW  - Prococol - ' ,subjcell{subjn}, '  ', Preprocessingcell{prepn}],'Interpreter','none')
%         colorbar('south', 'Color', [1 1 1]);
%     end
% end
% set(gcf,'position',pos)
% 
% figureJ(6)
% close(6)
% figureJ(6)
% for subjn = 1:length(subjcell)
%     for prepn = 1:length(Preprocessingcell)
%         nexttile()
%         Orthoview2(M0_MCR_data(subjn,prepn).img.*mask_data(subjn,prepn).img,[],[],'tight');
%         title(['M0 - Prococol - ' ,subjcell{subjn}, '  ', Preprocessingcell{prepn}],'Interpreter','none')
%         colorbar('south', 'Color', [1 1 1]);
%     end
% end
% set(gcf,'position',pos)
% 
% 
% %%
% 
% figureJ(12)
% close(12)
% figureJ(12)
% for subjn = 1:length(subjcell)
%     % for prepn = 1:length(Preprocessingcell)
%         nexttile()
%         Orthoview2(R1mapSCR_data(subjn).img.*mask_data(subjn,prepn).img,[],[0 1600],'tight')
%         title(['Prococol - ' ,subjcell{subjn}],'Interpreter','none')
%         colorbar('south', 'Color', [1 1 1]);
%     % end
% end
% set(gcf,'position',pos)
% 
% 
% figureJ(13)
% close(13)
% figureJ(13)
% for subjn = 1:length(subjcell)
%     % for prepn = 1:length(Preprocessingcell)
%         nexttile()
%         R2smapSCR_data(subjn).img(isnan(R2smapSCR_data(subjn).img))=0;
%         R2smapSCR_data(subjn).img(isinf(R2smapSCR_data(subjn).img))=0;
% 
%         Orthoview2(R2smapSCR_data(subjn).img.*mask_data(subjn,prepn).img,[],[0 50],'tight')
%         title(['Prococol - ' ,subjcell{subjn}, '  '],'Interpreter','none')
%         colorbar('south', 'Color', [1 1 1]);
%     % end
% end
% set(gcf,'position',pos)
% 
% figureJ(14)
% close(14)
% figureJ(14)
% for subjn = 1:length(subjcell)
%     % for prepn = 1:length(Preprocessingcell)
%         nexttile()
%         Orthoview2(M0mapSCR_data(subjn).img.*mask_data(subjn,prepn).img,[],[],'tight')
%         title(['Prococol - ' ,subjcell{subjn} ],'Interpreter','none')
%     % end
% end
% set(gcf,'position',pos)
% 
% %%
% RawData

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
