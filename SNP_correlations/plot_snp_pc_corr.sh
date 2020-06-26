#$ -S /bin/sh
#$ -N plot_snp_pc_corr
#$ -cwd
#$ -j y

module load conda

R -q --no-save --args --infile-prefix /home/hkings/DATA/PC_correlations6,18/snp_corr --outfile /home/hkings/DATA/PC_correlations6,18/snp_corr.png < plot_snp_pc_corr.R

#Run with: qsub -q new.q -v R_LIBS=/home/amstilp/devel/analysis_pipeline_cff_wgs/R_library calculate_snp_pc_corr.sh plot_snp_pc_corr.sh
