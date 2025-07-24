

bids_dir = '/project/3055010.04/RunningProjects/MyelinWaterImaging/bidsPhilipsVariants'

preprocessing = '0';
SepiaPrep = '0';
fittingMCR='3'; 
fittingMCRGPU='0'; 
writingMCR='3' ;
acqname = 'fl3d'
run_label = 'run-1';
subj_label = 'sub-003'
Macro_all(bids_dir, preprocessing, SepiaPrep, fittingMCR, fittingMCRGPU, writingMCR, acqname, run_label,subj_label)

