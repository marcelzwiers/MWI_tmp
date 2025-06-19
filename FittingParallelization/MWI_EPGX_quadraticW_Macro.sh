#! /bin/bash
#
# GRE_Macro.sh
#
# This is a shell script for Improved method fo in vivo MWI
# Perform following processing steps on GRE data:
# (1) rename to standard format;
# (2) brain extraction;
#
# Dependencies: (1)fsl
#
# Creator: Kwok-shing Chan @DCCN
# k.chan@donders.ru.nl
# Date created: 7 October 2018
# Date last edit:
############################################################

export subj=003
export run=1
export ses=mri02

export subj_label=sub-${subj}
export sess_label=ses-${ses}
export run_label=run-${run}

# call diretcory
. ../../subject_directory_master.sh 

export subj_script_dir=${code_dir}${subj_label}/${sess_label}/

# GRE Protocols and flip angles
protocol=( 'TR38NTE12' 'TR50NTE15' 'TR55NTE13' )
flipAngle=( 'FA5' 'FA10' 'FA20' 'FA50' 'FA70' )

# cprot=2

# ################################################################################
# ## Step 5: MWI

## EPGX, quadratic
for ((cprot=0;cprot<=2;cprot++)); do

export tmp_output_dir=${derivative_MWI_dir}MCR/no_preprocessing/using_5_flipangle/tmp_quadraticW_${protocol[${cprot}]}/

# create directories
mkdir -p ${tmp_output_dir}

for ((i=1;i<=160;i++)); do
echo ${protocol[${cprot}]}
echo ${i}
cat ${script_dir}MWI_01_MCR_5FA_quadraticW_01_slice.m > ${tmp_output_dir}MWI_01_MCR_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "4i subj='${subj}' ;" ${tmp_output_dir}MWI_01_MCR_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "6i ses='${ses}' ;" ${tmp_output_dir}MWI_01_MCR_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "26i protocol= { '${protocol[${cprot}]}' } ;" ${tmp_output_dir}MWI_01_MCR_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "28i slice=${i};" ${tmp_output_dir}MWI_01_MCR_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m

if [ $i -le 47 ]; then
matlab_sub --walltime 01:00:00 --mem 10gb ${tmp_output_dir}MWI_01_MCR_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m
else
matlab_sub --walltime 20:00:00 --mem 10gb ${tmp_output_dir}MWI_01_MCR_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m
sleep 1s
fi

if [ $i -eq 47 ]; then
sleep 3m
fi

if [ $i -eq 99 ]; then
sleep 3.5h
fi

if [ $i -eq 160 ]; then
sleep 3.5h
fi

done
done

# ################################################################################
# ## Step 5: MCR-DIMWI

## EPGX, quadratic
for ((cprot=0;cprot<=2;cprot++)); do

export tmp_output_dir=${derivative_MWI_dir}MCRDIMWI/no_preprocessing/using_5_flipangle/tmp_quadraticW_${protocol[${cprot}]}/

# create directories
mkdir -p ${tmp_output_dir}

for ((i=1;i<=160;i++)); do

echo ${i}
cat ${script_dir}MWI_01_MCRDIMWI_5FA_quadraticW_01_slice.m > ${tmp_output_dir}MWI_01_MCRDIMWI_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "4i subj='${subj}' ;" ${tmp_output_dir}MWI_01_MCRDIMWI_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "6i ses='${ses}' ;" ${tmp_output_dir}MWI_01_MCRDIMWI_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "26i protocol= { '${protocol[${cprot}]}' } ;" ${tmp_output_dir}MWI_01_MCRDIMWI_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "28i slice=${i};" ${tmp_output_dir}MWI_01_MCRDIMWI_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m

if [ $i -le 47 ]; then
matlab_sub --walltime 01:00:00 --mem 10gb ${tmp_output_dir}MWI_01_MCRDIMWI_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m
else
matlab_sub --walltime 25:00:00 --mem 10gb ${tmp_output_dir}MWI_01_MCRDIMWI_5FA_quadraticW_01_${protocol[${cprot}]}_01_slice_${i}.m
sleep 1s
fi

if [ $i -eq 47 ]; then
sleep 3m
fi

if [ $i -eq 99 ]; then
sleep 5h
fi

if [ $i -eq 160 ]; then
sleep 5h
fi

done
done

# ################################################################################
# ## Step 6: MWI, less echo for the second protocol
cprot=1

## EPGX, quadratic
export tmp_output_dir=${derivative_MWI_dir}MCR/no_preprocessing/using_5_flipangle/lessTE/tmp_quadraticW_${protocol[${cprot}]}/

# create directories
mkdir -p ${tmp_output_dir}

for ((i=1;i<=160;i++)); do

echo ${i}
cat ${script_dir}MWI_01_MCR_5FA_quadraticW_lessTE_01_slice.m > ${tmp_output_dir}MWI_01_MCR_5FA_qW_lessTE_01_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "4i subj='${subj}' ;" ${tmp_output_dir}MWI_01_MCR_5FA_qW_lessTE_01_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "6i ses='${ses}' ;" ${tmp_output_dir}MWI_01_MCR_5FA_qW_lessTE_01_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "26i protocol= { '${protocol[${cprot}]}' } ;" ${tmp_output_dir}MWI_01_MCR_5FA_qW_lessTE_01_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "28i slice=${i};" ${tmp_output_dir}MWI_01_MCR_5FA_qW_lessTE_01_${protocol[${cprot}]}_01_slice_${i}.m

if [ $i -le 47 ]; then
matlab_sub --walltime 01:00:00 --mem 10gb ${tmp_output_dir}MWI_01_MCR_5FA_qW_lessTE_01_${protocol[${cprot}]}_01_slice_${i}.m
else
matlab_sub --walltime 20:00:00 --mem 10gb ${tmp_output_dir}MWI_01_MCR_5FA_qW_lessTE_01_${protocol[${cprot}]}_01_slice_${i}.m
sleep 1s
fi

if [ $i -eq 47 ]; then
sleep 3m
fi

if [ $i -eq 99 ]; then
sleep 3.5h
fi

if [ $i -eq 160 ]; then
sleep 3.5h
fi

done

# ################################################################################
# ## Step 7: MWI, 3FA, first protocol
cprot=0
# # (0) [2,3,5]=10,20,70 degrees; 1st place
# # (1) [1,2,5]=5,10,70 degrees; 2nd place
# # (2) [2,3,4]=10,20,50 degrees; 3rd place
# # (3) [1,2,3]=5,10,20 degrees; really bad option :(
# # fa_options=( '[2,3,5]' '[1,2,5]' '[2,3,4]' '[1,2,3]' )
# cfa=3

## EPGX, quadratic
for ((cfa=0;cfa<=3;cfa++)); do

export tmp_output_dir=${derivative_MWI_dir}MCR/no_preprocessing/using_3_flipangle/tmp_quadraticW_${protocol[${cprot}]}_faused_option${cfa}/

# create directories
mkdir -p ${tmp_output_dir}

for ((i=1;i<=160;i++)); do

echo ${i}
cat ${script_dir}MWI_01_MCR_3FA_quadraticW_01_slice.m > ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "4i subj='${subj}' ;" ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "6i ses='${ses}' ;" ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "27i cfa= ${cfa} +1 ;" ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "28i slice=${i};" ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m

if [ $i -le 47 ]; then
matlab_sub --walltime 01:00:00 --mem 10gb ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m
else
matlab_sub --walltime 15:00:00 --mem 10gb ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m
sleep 1s
fi

if [ $i -eq 47 ]; then
sleep 3m
fi

if [ $i -eq 99 ]; then
sleep 2.5h
fi

if [ $i -eq 160 ]; then
sleep 2.5h
fi

done
done

################################################################################
## Step 7: MWI, 3FA, first protocol
cfa=0
# # (0) [2,3,5]=10,20,70 degrees; 1st place
# # (1) [1,2,5]=5,10,70 degrees; 2nd place
# # (2) [2,3,4]=10,20,50 degrees; 3rd place
# # (3) [1,2,3]=5,10,20 degrees; really bad option :(
# # fa_options=( '[2,3,5]' '[1,2,5]' '[2,3,4]' '[1,2,3]' )
# cfa=3

## EPGX, quadratic
for ((cprot=0;cprot<=2;cprot++)); do

export tmp_output_dir=${derivative_MWI_dir}MCR/no_preprocessing/using_3_flipangle/tmp_quadraticW_${protocol[${cprot}]}_faused_option${cfa}/

# create directories
mkdir -p ${tmp_output_dir}

for ((i=1;i<=160;i++)); do

echo ${i}
cat ${script_dir}MWI_01_MCR_3FA_quadraticW_01_slice.m > ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "4i subj='${subj}' ;" ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "6i ses='${ses}' ;" ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "27i protocol= { '${protocol[${cprot}]}' } ;" ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "28i cfa= ${cfa} +1 ;" ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m
sed -i "29i slice=${i};" ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m

if [ $i -le 47 ]; then
matlab_sub --walltime 01:00:00 --mem 10gb ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m
else
matlab_sub --walltime 15:00:00 --mem 10gb ${tmp_output_dir}MWI_01_MCR_3FA_qW_01_fa_option${cfa}_${protocol[${cprot}]}_01_slice_${i}.m
sleep 1s
fi

if [ $i -eq 47 ]; then
sleep 3m
fi

if [ $i -eq 99 ]; then
sleep 2.5h
fi

if [ $i -eq 160 ]; then
sleep 2.5h
fi

done
done

