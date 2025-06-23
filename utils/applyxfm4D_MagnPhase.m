function res = applyxfm4D_MagnPhase(input);


% input.MagName - should have _mag in the end of the name, prior to .nii
% input.PhaseName - should have _phase in the end of the name, prior to .nii
% input.MotionMat - directory where motion matrices computed by mcflirt are stored
% input.RefVol - referenceVolume
%
%         MagnData = [Dataset{count}.outdir ,'Data_',Dataset{count}.echo{count_echo},'_mag.nii.gz']
%         PhaseData = [Dataset{count}.outdir,'Data_',Dataset{count}.echo{count_echo},'_phase.nii.gz']
%         RealData = [Dataset{count}.outdir ,'Data_',Dataset{count}.echo{count_echo},'_real.nii.gz']
%         ImagData = [Dataset{count}.outdir ,'Data_',Dataset{count}.echo{count_echo},'_imag.nii.gz']
%         RealDataMC = [Dataset{count}.outdir ,'Data_',Dataset{count}.echo{count_echo},'_real_mcf.nii.gz'];
%         ImagDataMC = [Dataset{count}.outdir ,'Data_',Dataset{count}.echo{count_echo},'_imag_mcf.nii.gz'];
%         
%         MagDataMC = [Dataset{count}.outdir ,'Data_',Dataset{count}.echo{count_echo},'_mag_mcf.nii.gz'];
%         PhaseDataMC = [Dataset{count}.outdir ,'Data_',Dataset{count}.echo{count_echo},'_phase_mcf.nii.gz'];

debugmode = 0;
MagnData = input.MagName;
PhaseData = input.PhaseName;

StartMagString = findstr(MagnData,'-mag');
StartPhaseString = findstr(PhaseData,'-phase');

RealData = [MagnData(1:StartMagString-1),'-real.nii.gz'];
ImagData = [MagnData(1:StartMagString-1),'-imag.nii.gz'];
RealDataMC = [MagnData(1:StartMagString-1),'-real_mcf.nii.gz'];
ImagDataMC = [MagnData(1:StartMagString-1),'-imag_mcf.nii.gz'];

if and(isfield(input,'MagNameOutput'),isfield(input,'PhaseNameOutput'))
    MagDataMC = input.MagNameOutput;
    PhaseDataMC = input.PhaseNameOutput;
else
    MagDataMC = [MagnData(1:StartMagString-1),'-mag_mcf.nii.gz'];
    PhaseDataMC = [PhaseData(1:StartPhaseString-1),'-phase_mcf.nii.gz'];
end
complexdata = loadcomplex_nii( MagnData , PhaseData );
save_nii_quick(complexdata, real(complexdata.img) , RealData)
save_nii_quick(complexdata, imag(complexdata.img) , ImagData)

if isfolder(input.MotionMat) 
    unix([ 'applyxfm4D ',RealData, ' ',input.RefVol,' ', RealDataMC, ' ',input.MotionMat,' -interp sinc'] );
    unix([ 'applyxfm4D ',ImagData, ' ',input.RefVol,' ', ImagDataMC, ' ',input.MotionMat,' -interp sinc'] );
else
    unix([ 'flirt -in ',RealData, ' -ref ',input.RefVol,' -out ', RealDataMC, ' -applyxfm -init ',input.MotionMat,' -interp sinc'] );
    unix([ 'flirt -in ',ImagData, ' -ref ',input.RefVol,' -out ', ImagDataMC, ' -applyxfm -init ',input.MotionMat,' -interp sinc'] );
end

realdata = load_untouch_nii( RealDataMC );
imagdata = load_untouch_nii( ImagDataMC );
complexdata = realdata;
complexdata.img = (single(realdata.img) + i * single(imagdata.img));
save_nii_quick(complexdata, abs(complexdata.img) , MagDataMC);
save_nii_quick(complexdata, angle(complexdata.img) , PhaseDataMC);


% if debugmode
%     figure(1)
%     subplot(321)
%     Orthoview(realdata.img(:,:,:,1)); colorbar; title('real')
%     subplot(322)
%     Orthoview(imagdata.img(:,:,:,1)); colorbar; title('imag')
%     complexdata.img = angle(single(realdata.img) + i * single(imagdata.img)) ;
%     subplot(323)
%     Orthoview(complexdata.img(:,:,:,1)); colorbar; title('phase')
%     complexdata.img = abs(single(realdata.img) + i * single(imagdata.img)) ;
%     subplot(324)
%     Orthoview(complexdata.img(:,:,:,1)); colorbar; title('abs')
% 
%     realdata_old = load_untouch_nii( RealData,[1] );
%     imagdata_old = load_untouch_nii( ImagData,[1] );
%     subplot(325)
%     Orthoview(realdata_old.img(:,:,:,1)); colorbar; title('real pre mc')
%     subplot(326)
%     Orthoview(imagdata_old.img(:,:,:,1)); colorbar; title('imag pre mc')
%     magdata= load_untouch_nii( MagDataMC,[1] );
%     phasedata = load_untouch_nii( PhaseDataMC,[1] );
%     subplot(325)
%     Orthoview(phasedata.img(:,:,:,1)); colorbar; title('phase post saving')
%     subplot(326)
%     Orthoview(magdata.img(:,:,:,1)); colorbar; title('mag post saving')
% 
% end;
% 
unix(['rm ',RealData])
unix(['rm ',ImagData])
unix(['rm ',RealDataMC])
unix(['rm ',ImagDataMC])

res = 1;
