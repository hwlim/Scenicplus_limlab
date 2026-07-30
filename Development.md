# Development Log

20260512: Initial test
  - /Volumes/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus
  - Failed in 07_run_scenicplus.log
    ```bash
    26 [Tue May 12 22:19:35 2026]
    27 localrule download_genome_annotations:
    28     output: /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/genome_annotation.tsv, /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/chromsizes.tsv
    29     jobid: 8
    30     reason: Missing output files: /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/genome_annotation.tsv, /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/chromsizes.tsv
    31     resources: tmpdir=/scratch/limc8h
    32 
    33 /data/limlab/Resource/conda_env/scenicplus_limlab/lib/python3.11/site-packages/pybiomart/dataset.py:269: DtypeWarning: Columns (0) have mixed types. Specify dtype option on import or set low_memory=False.
    34   result = pd.read_csv(StringIO(response.text), sep='\t')
    35 2026-05-12 22:20:09,890 Download gene annotation INFO     Using genome: GRCm39
    36 Could not find Id on https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&term=GRCm39
    37 Returning gene annotation without subestting for assembled chromosomesand converting to UCSC style. Please make sure that the chromosome namesin the returned object match with the chromosome names in the scplus_obj.Chromosome sizes will not be returned
    38 2026-05-12 22:20:10,163 SCENIC+      INFO     Chrosomome sizes was not found, please provide this information manually.
    39 2026-05-12 22:20:10,164 SCENIC+      INFO     Saving genome annotation to: /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/genome_annotation.tsv
    40 Waiting at most 5 seconds for missing files.
    41 MissingOutputException in rule download_genome_annotations in file /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_pipeline/Snakemake/workflow/Snakefile, line 221:
    42 Job 8  completed successfully, but some output files are missing. Missing files after 5 seconds. This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait:
    43 /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/chromsizes.tsv
    44 Removing output files of failed job download_genome_annotations since they might be corrupted:
    45 /data/parklabngs/Lim.0.Analysis/20240521-multiome-kidney/Seurat/merged_macs/ENP/ScenicPlus/results/scplus_out/genome_annotation.tsv
    46 Shutting down, this might take some time.
    47 Exiting because a job execution failed. Look above for error message
    48 Complete log: .snakemake/log/2026-05-12T221935.788913.snakemake.log
    49 WorkflowError:
    50 At least one job did not complete successfully.
    ```