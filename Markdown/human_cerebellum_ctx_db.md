# Creating custom cistarget database

In this tutorial we will create a custom cistarget database using consensus peaks.

This involves precomputed scores for all the motifs in our motif collection on a predefined set of regions

We provide precomputed databases for [human](https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/), [mouse](https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/) and [fly](https://resources.aertslab.org/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc_v10_clust/region_based/). These databases are computed on regulatory regions spanning the genome. Feel free to use these databases, however for the best results we recommend to generate a custom database given that it is highly likely that the precomputed databases don't cover all the regions in your consensus peak set.

## Download create_cistarget_database

We will start by downloading and installing the `create_cistarget_database` repository.


```python
cd /staging/leuven/stg_00002/lcb/sdewin/PhD/python_modules/scenicplus_development_tutorial/ctx_db
source /staging/leuven/stg_00002/mambaforge/vsc33053/etc/profile.d/conda.sh
conda activate scenicplus_development_tutorial
```


```python
git clone https://github.com/aertslab/create_cisTarget_databases
```

## Download cluster-buster

[Cluster-buster](https://github.com/weng-lab/cluster-buster) will be used to score the regions using our motif collection. We provide a precompiled binary of cluster buster.



```python
wget https://resources.aertslab.org/cistarget/programs/cbust
chmod a+x cbust
```


## Download motif collection

Next, we will download the motif collection.


```python
mkdir -p aertslab_motif_colleciton
wget -O aertslab_motif_colleciton/v10nr_clust_public.zip https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip
```


```python
cd aertslab_motif_colleciton; unzip -q v10nr_clust_public.zip
cd ..
```

These are the motif-to-TF annotations for:

- Chicken: motifs-v10-nr.chicken-m0.00001-o0.0.tbl
- fly: motifs-v10-nr.flybase-m0.00001-o0.0.tbl
- human: motifs-v10-nr.hgnc-m0.00001-o0.0.tbl
- mouse: motifs-v10-nr.mgi-m0.00001-o0.0.tbl


```python
ls aertslab_motif_colleciton/v10nr_clust_public/snapshots/
```

Here are some example motifs, they are stored in cb format.


```python
ls -l aertslab_motif_colleciton/v10nr_clust_public/singletons | head
```



```python
cat aertslab_motif_colleciton/v10nr_clust_public/singletons/bergman__Adf1.cb
```


## Prepare fasta from consensus regions

Next we will get sequences for all the consensus peaks. We will also add 1kb of background padding, this will be used as a background sequence for cluster-buster. It is completely optional to add this padding, we have noticed that it does not affect the analyses a lot.


```python
module load cluster/wice/bigmem
module load BEDTools/2.30.0-GCC-10.3.0

REGION_BED="/staging/leuven/stg_00002/lcb/sdewin/PhD/python_modules/pycisTopic_polars_tutorial/outs/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="/staging/leuven/res_00001/genomes/homo_sapiens/hg38_ucsc/fasta/hg38.fa"
CHROMSIZES="/staging/leuven/res_00001/genomes/homo_sapiens/hg38_ucsc/fasta/hg38.chrom.sizes"
DATABASE_PREFIX="10x_brain_1kb_bg_with_mask"
SCRIPT_DIR="/staging/leuven/stg_00002/lcb/sdewin/PhD/python_modules/scenicplus_development_tutorial/ctx_db/create_cisTarget_databases"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        hg38.10x_brain.with_1kb_bg_padding.fa \
        1000 \
        yes
```




```python
head -n 2 hg38.10x_brain.with_1kb_bg_padding.fa
```


## Create cistarget databases

Now we can create the ranking and score database. This step will take some time so we recommend to run it as a job (i.e. not in jupyter notebooks).


```python
ls aertslab_motif_colleciton/v10nr_clust_public/singletons > motifs.txt
```


```python
OUT_DIR=""${PWD}""
CBDIR="${OUT_DIR}/aertslab_motif_colleciton/v10nr_clust_public/singletons"
FASTA_FILE="${OUT_DIR}/hg38.10x_brain.with_1kb_bg_padding.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 20
```
