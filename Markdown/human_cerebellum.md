# Running SCENIC+

This tutorial will illustrate how to run the SCENIC+ pipeline via Snakemake.

Before running this pipeline you should:

- Preprocess the scATAC-seq side of the data using [pycisTopic](https://github.com/aertslab/pycisTopic), [click here](...) for a tutorial.
- Preprocess the scRNA-seq side of the data using [Scanpy](https://github.com/scverse/scanpy), [click here](...) for a tutorial.
- Optionally, but highly recommended, generate a cisTarget database using the consensus peaks specific to your dataset, [click here](...) for a tutorial.

In case you have human, mouse or fly data you can also use one of the precomputed cisTarget databases. These can be found on our [resources website](https://resources.aertslab.org/cistarget/databases/).


```python
import os
os.chdir("/staging/leuven/stg_00002/lcb/sdewin/PhD/python_modules/scenicplus_development_tutorial")
```

SCENIC+ can be run entirely using the command line. Please refer to the docummentation (by using the `--help` flag) for more information.

Below we will use the Snakemake pipeline, which is already included in SCENIC+ to perform the analysis.


```python
!scenicplus
```

    
       ____   ____ _____ _   _ ___ ____      
      / ___| / ___| ____| \ | |_ _/ ___|[31;1m _ [0m
      \___ \| |   |  _| |  \| || | |  [31;1m _|.|_[0m
       ___) | |___| |___| |\  || | |__[31;1m|_..._|[0m
      |____/ \____|_____|_| \_|___\____|[31;1m|_|[0m 
    
    
    scenicplus verions: 1.0a1
    usage: scenicplus [-h] {init_snakemake,prepare_data,grn_inference} ...
    
    Single-Cell Enhancer-driven gene regulatory Network Inference and Clustering
    
    positional arguments:
      {init_snakemake,prepare_data,grn_inference}
    
    options:
      -h, --help            show this help message and exit


### scRNA-seq preparation

scRNA-seq side of the experiment can be processed according to the regular [Scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) tutorial. Just make sure to store the raw gene expresison matrix in `adata.raw`.

**Call the following piece of code:**

```python

adata.raw = adata

```

**_BEFORE_** normalizing the data

```python

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

```


```python
ls data
```

    [0m[01;36m10x_multiome_brain_cisTopicObject_noDBL.pkl[0m@  [01;36madata.h5ad[0m@  [01;34mregion_sets[0m/


## Initialize Snakemake

To run Snakemake we first have to initialize the pipeline. This will create a folder named `Snakemake` containing a folder for the `config.yaml` file and a folder containing the actual workflow definition.


```python
!mkdir -p scplus_pipeline
!scenicplus init_snakemake --out_dir scplus_pipeline
```

    2024-03-06 14:28:10,766 SCENIC+      INFO     Creating snakemake folder in: scplus_pipeline
    [0m


```python
!tree scplus_pipeline/
```

    scplus_pipeline/
    └── Snakemake
        ├── config
        │   └── config.yaml
        └── workflow
            └── Snakefile
    
    3 directories, 2 files



```python
!mkdir -p outs
!mkdir -p tmp
```

Next, modify the `config.yaml` file.

The most important fields are `input_data` and `output_data`. For the other default values can be kept:

**input_data**:
This is the input to the pipeline, these files should already exist.

- `cisTopic_obj_fname`: the path to your cistopic object containing processed chromatin accessibility data.
- `GEX_anndata_fname`: the path to your scanpy h5ad file containing processed gene expression data.
- `region_set_folder`: the path to directory containing several directories with bed files. Differential motif enrichment (1 vs all) will be run within each sub folder. As an example the structure of the folder is shown below.
- `ctx_db_fname`: the path to the cisTarget **ranking** database.
- `dem_db_fname`: the path to the cisTarget **score** database.
- `path_to_motif_annotations`: the path to the motif-to-TF annotaiton. For human (hgnc), mouse (mgi), chicken and fly (flybase) these files can be downloaded from [our resources website](https://resources.aertslab.org/cistarget/motif2tf/), please download the relevant file starting with "motifs-v10nr_clust".

**output_data**:
This is the output of the pipeline, these files will be created.
If some of these files already exists (for example when the pipeline has only been partially run) some steps of the workflow might be skipped.

- `combined_GEX_ACC_mudata`: where the [MuData](https://github.com/scverse/mudata) object containing gene expression and imputed chromatin accessibility should be stored.
- `dem_result_fname`: where the h5 file containing DEM based enriched motifs should be stored.
- `ctx_result_fname`: where the h5 file containing cistarget based enriched motifs should be stored.
- `output_fname_dem_html`: where the html file containing DEM based enriched motifs should be stored.
- `output_fname_ctx_html`: where the html file containing cistarget based enriched motifs should be stored.
- `cistromes_direct`: where the [AnnData](https://github.com/scverse/anndata) h5ad file should be stored containing TF-to-region links based on direct motif-to-TF annotations.
- `cistromes_extended`: where the [AnnData](https://github.com/scverse/anndata) h5ad file should be stored containing TF-to-region links based on exteded (e.g. orthology based) motif-to-TF annotations.
- `tf_names`: where a text file containing TF names, based on the enriched motifs, should be stored.
- `genome_annotation`: where a data frame (tsv) should be stored containing genome annotation.
- `chromsizes`: where the chromsizes file should be stored.
- `search_space`: where the search space for each gene should be stored.
- `tf_to_gene_adjacencies`: where the TF-to-gene links, with importance scores, should be stored.
- `region_to_gene_adjacencies`: where the region-to-gene links, with importance scores, should be stored.
- `eRegulons_direct`: where the dataframe (tsv) containing eRegulons (TF-region-gene links) based on direct motif-to-TF annotations should be stored.
- `eRegulons_extended`: where the dataframe (tsv) containing eRegulons (TF-region-gene links) based on extended (e.g. orthology based) motif-to-TF annotations should be stored.
- `AUCell_direct`: where the [MuData](https://github.com/scverse/mudata) containing target gene and target region enrichement scores for each cells, based on direct motif-to-TF annotations should be stored.
- `AUCell_extended`: where the [MuData](https://github.com/scverse/mudata) containing target gene and target region enrichement scores for each cells, based on extended (e.g. orthology based) motif-to-TF annotations should be stored.
- `scplus_mdata`: where the final output [MuData](https://github.com/scverse/mudata) containing AUCell values and the (TF-region-gene links) based on both direct and extended motif-to-TF annotations should be stored.

**params_general**
General parameters.
- `temp_dir`: Directory to store temporary data.
- `n_cpu`: maximum number of CPU's to use.
- `seed`: seed to use to initialize the random state.

**params_data_preparation**
Parameters used for the data preparation step.
- `bc_transform_func`: lambda function to transform the scRNA-seq barcode so they match with the scATAC-seq ones
- `is_multiome`: boolean specifying wether the data is multiome or not.
- `key_to_group_by`: in case of non-multiome data, cell metadata variable to group cells by in order to generate metacells that can be matched across the scRNA-seq and scATAC-seq side of the data. This variable should be prefixed with eiter "GEX:" or "ACC:".
- `nr_cells_per_metacells`: in case of non-multione data, number of cells to sample to sample for each metacell.
- `direct_annotation`: Which annotations fields to use for generating direct motif-to-TF annotations.
- `extended_annotation`: Which annotations fields to use for generating extended motif-to-TF annotations
- `species`: Species name, for example "hsapiens"
- `biomart_host`: Biomart host to use for downloading genome annotations. Make sure that this host matches the genome reference you are using, please visit [this website](https://www.ensembl.org/info/website/archives/index.html) for more information.
- `search_space_upstream`: string in the form "\<minmal\> \<maximal\>" specifying the \<minimal\> and \<maximal\> search space to consider downstream of the TSS of each gene.
- `search_space_downstream`: string in the form "\<minmal\> \<maximal\>" specifying the \<minimal\> and \<maximal\> search space to consider upstrean of the TSS of each gene.
- `search_space_extend_tss`: string in the form "\<upstream\> \<downstream\>" specifying the amount of basepairs the TSS of each gene should be extended, \<upstream\>  and 

**params_motif_enrichment**
parameters for performing motif enrichment analysis.
- `species`: Species used for the analysis. This parameter is used to download the correct motif-to-TF annotations from the cisTarget webservers.
- `annotation_version`: Version of the motif-to-TF annotation to use. This parameter is used to download the correct motif-to-TF data from the cisTarget webservers.
- `motif_similarity_fdr`: Threshold on motif similarity scores for calling similar motifs.
- `orthologous_identity_threshold`: Threshold on the protein-protein orthology score for calling orthologous motifs
- `annotations_to_use`: Which annotations to use for annotation motifs to TFs.
- `fraction_overlap_w_dem_database`: Fraction of nucleotides, of regions in the bed file, that should overlap with regions in the scores database.
- `dem_max_bg_regions`: Maximum number of regions to use as background for DEM.
- `dem_balance_number_of_promoters`: Boolean specifying wether the number of promoters should be equalized between the foreground and background set of regions.
- `dem_promoter_space`: Number of basepairs up- and downstream of the TSS that are considered as being the promoter for that gene.
- `dem_adj_pval_thr`: Threshold on the Benjamini-Hochberg adjusted p-value from the Wilcoxon test performed on the motif score of foreground vs background regions for a motif to be considered as enriched.
- `dem_log2fc_thr`: Threshold on the log2 fold change of the motif score of foreground vs background regions for a motif to be considered as enriched.
- `dem_mean_fg_thr`: Minimul mean signal in the foreground to consider a motif enriched for DEM.
- `dem_motif_hit_thr`: Minimal CRM score to consider a region enriched for a motif for DEM.
- `fraction_overlap_w_ctx_database`: Fraction of nucleotides, of regions in the bed file, that should overlap with regions in the ranking database.
- `ctx_auc_threshold`: Threshold on the AUC value for calling significant motifs
- `ctx_nes_threshold`: Threshold on the NES value for calling significant motifs.
- `ctx_rank_threshold`: The total number of ranked regions to take into account when creating a recovery curves.

**params_inference**
Parameters for performing GRN inference.
- `tf_to_gene_importance_method`: Method to use to calculate TF-to-gene importance scores.
- `region_to_gene_importance_method`: Method to use to calculate region-to-gene importance scores.
- `region_to_gene_correlation_method`: Method to use to calculate region-to-gene correlation coefficients.
- `order_regions_to_genes_by`: Value to order region-to-gene scores by for selecting top regions per gene
- `order_TFs_to_genes_by`: value to order TF-to-gene scores by for selecting top TFs per gene
- `gsea_n_perm`: Number or permutations to perform for calculating GSEA enrichment scores.
- `quantile_thresholds_region_to_gene`: space seperated list containing quantile threshold to be used for binarizing region-to-gene links.
- `top_n_regionTogenes_per_gene`: space seperated list containing the number of top regions per gene for binarizing region-to-gene links.
- `top_n_regionTogenes_per_region`: space seperated list containging per region the number of top genes for binarizing region-to-gene links
- `min_regions_per_gene`: minimum number of regions per gene for the link to be included in eGRNs.
- `rho_threshold`: absolute threshold on the correlation coefficient to seperate positive and negative region-to-gene and TF-to-gene links
- `min_target_genes`: minimum number of target genes per TF for the link(s) to be includedin eGRNs.



```python
!tree /staging/leuven/stg_00002/lcb/sdewin/PhD/python_modules/pycisTopic_polars_tutorial/outs/region_sets
```

    /staging/leuven/stg_00002/lcb/sdewin/PhD/python_modules/pycisTopic_polars_tutorial/outs/region_sets
    ├── DARs_cell_type
    │   ├── AST.bed
    │   ├── BG.bed
    │   ├── COP.bed
    │   ├── ENDO.bed
    │   ├── GC.bed
    │   ├── GP.bed
    │   ├── INH_PVALB.bed
    │   ├── INH_SNCG.bed
    │   ├── INH_SST.bed
    │   ├── INH_VIP.bed
    │   ├── MG.bed
    │   ├── MGL.bed
    │   ├── MOL.bed
    │   ├── NFOL.bed
    │   ├── OPC.bed
    │   └── PURK.bed
    ├── Topics_otsu
    │   ├── Topic1.bed
    │   ├── Topic10.bed
    │   ├── Topic11.bed
    │   ├── Topic12.bed
    │   ├── Topic13.bed
    │   ├── Topic14.bed
    │   ├── Topic15.bed
    │   ├── Topic16.bed
    │   ├── Topic17.bed
    │   ├── Topic18.bed
    │   ├── Topic19.bed
    │   ├── Topic2.bed
    │   ├── Topic20.bed
    │   ├── Topic21.bed
    │   ├── Topic22.bed
    │   ├── Topic23.bed
    │   ├── Topic24.bed
    │   ├── Topic25.bed
    │   ├── Topic26.bed
    │   ├── Topic27.bed
    │   ├── Topic28.bed
    │   ├── Topic29.bed
    │   ├── Topic3.bed
    │   ├── Topic30.bed
    │   ├── Topic31.bed
    │   ├── Topic32.bed
    │   ├── Topic33.bed
    │   ├── Topic34.bed
    │   ├── Topic35.bed
    │   ├── Topic36.bed
    │   ├── Topic37.bed
    │   ├── Topic38.bed
    │   ├── Topic39.bed
    │   ├── Topic4.bed
    │   ├── Topic40.bed
    │   ├── Topic5.bed
    │   ├── Topic6.bed
    │   ├── Topic7.bed
    │   ├── Topic8.bed
    │   └── Topic9.bed
    └── Topics_top_3k
        ├── Topic1.bed
        ├── Topic10.bed
        ├── Topic11.bed
        ├── Topic12.bed
        ├── Topic13.bed
        ├── Topic14.bed
        ├── Topic15.bed
        ├── Topic16.bed
        ├── Topic17.bed
        ├── Topic18.bed
        ├── Topic19.bed
        ├── Topic2.bed
        ├── Topic20.bed
        ├── Topic21.bed
        ├── Topic22.bed
        ├── Topic23.bed
        ├── Topic24.bed
        ├── Topic25.bed
        ├── Topic26.bed
        ├── Topic27.bed
        ├── Topic28.bed
        ├── Topic29.bed
        ├── Topic3.bed
        ├── Topic30.bed
        ├── Topic31.bed
        ├── Topic32.bed
        ├── Topic33.bed
        ├── Topic34.bed
        ├── Topic35.bed
        ├── Topic36.bed
        ├── Topic37.bed
        ├── Topic38.bed
        ├── Topic39.bed
        ├── Topic4.bed
        ├── Topic40.bed
        ├── Topic5.bed
        ├── Topic6.bed
        ├── Topic7.bed
        ├── Topic8.bed
        └── Topic9.bed
    
    3 directories, 96 files



```python
!bat scplus_pipeline/Snakemake/config/config.yaml
```

## Run pipeline

Once the config file is filled in the pipeline can be run.


```python
cd scplus_pipeline/Snakemake/
```


```python
ls
```

    config	workflow



```python
source /staging/leuven/stg_00002/mambaforge/vsc33053/etc/profile.d/conda.sh
conda activate scenicplus_development_tutorial
snakemake --cores 20
```

## Main outputs

The main output of the pipeline is the `scplusmdata.h5mu` file. This is a MuData file containing the eRegulons and enrichment scores.


```python
import os
os.chdir("/staging/leuven/stg_00002/lcb/sdewin/PhD/python_modules/scenicplus_development_tutorial")
```


```python
import mudata
scplus_mdata = mudata.read("outs/scplusmdata.h5mu")
```

Direct and extended predicted TF-to-region-to-gene links. This dataframe contains also a ranking of each TF-region-gene triplet, based on its importance `triplet_rank`.


```python
scplus_mdata.uns["direct_e_regulon_metadata"]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Region</th>
      <th>Gene</th>
      <th>importance_R2G</th>
      <th>rho_R2G</th>
      <th>importance_x_rho</th>
      <th>importance_x_abs_rho</th>
      <th>TF</th>
      <th>is_extended</th>
      <th>eRegulon_name</th>
      <th>Gene_signature_name</th>
      <th>Region_signature_name</th>
      <th>importance_TF2G</th>
      <th>regulation</th>
      <th>rho_TF2G</th>
      <th>triplet_rank</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr20:35305658-35306158</td>
      <td>CEP250</td>
      <td>0.020002</td>
      <td>0.056479</td>
      <td>0.001130</td>
      <td>0.001130</td>
      <td>BCL11A</td>
      <td>False</td>
      <td>BCL11A_direct_+/+</td>
      <td>BCL11A_direct_+/+_(352g)</td>
      <td>BCL11A_direct_+/+_(641r)</td>
      <td>0.743348</td>
      <td>1</td>
      <td>0.149413</td>
      <td>1879</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr4:41408947-41409447</td>
      <td>UCHL1</td>
      <td>0.039919</td>
      <td>0.278051</td>
      <td>0.011099</td>
      <td>0.011099</td>
      <td>BCL11A</td>
      <td>False</td>
      <td>BCL11A_direct_+/+</td>
      <td>BCL11A_direct_+/+_(352g)</td>
      <td>BCL11A_direct_+/+_(641r)</td>
      <td>1.572051</td>
      <td>1</td>
      <td>0.309145</td>
      <td>1767</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr10:133327408-133327908</td>
      <td>CALY</td>
      <td>0.021935</td>
      <td>0.542038</td>
      <td>0.011889</td>
      <td>0.011889</td>
      <td>BCL11A</td>
      <td>False</td>
      <td>BCL11A_direct_+/+</td>
      <td>BCL11A_direct_+/+_(352g)</td>
      <td>BCL11A_direct_+/+_(641r)</td>
      <td>0.934755</td>
      <td>1</td>
      <td>0.369021</td>
      <td>7337</td>
    </tr>
    <tr>
      <th>3</th>
      <td>chr1:239814303-239814803</td>
      <td>CHRM3</td>
      <td>0.054480</td>
      <td>0.731600</td>
      <td>0.039857</td>
      <td>0.039857</td>
      <td>BCL11A</td>
      <td>False</td>
      <td>BCL11A_direct_+/+</td>
      <td>BCL11A_direct_+/+_(352g)</td>
      <td>BCL11A_direct_+/+_(641r)</td>
      <td>1.682106</td>
      <td>1</td>
      <td>0.497771</td>
      <td>4406</td>
    </tr>
    <tr>
      <th>4</th>
      <td>chr10:13987335-13987835</td>
      <td>FRMD4A</td>
      <td>0.002173</td>
      <td>0.372226</td>
      <td>0.000809</td>
      <td>0.000809</td>
      <td>BCL11A</td>
      <td>False</td>
      <td>BCL11A_direct_+/+</td>
      <td>BCL11A_direct_+/+_(352g)</td>
      <td>BCL11A_direct_+/+_(641r)</td>
      <td>0.764116</td>
      <td>1</td>
      <td>0.411129</td>
      <td>6392</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>7648</th>
      <td>chr12:78102733-78103233</td>
      <td>NAV3</td>
      <td>0.061226</td>
      <td>-0.598159</td>
      <td>-0.036623</td>
      <td>0.036623</td>
      <td>TCF12</td>
      <td>False</td>
      <td>TCF12_direct_-/-</td>
      <td>TCF12_direct_-/-_(123g)</td>
      <td>TCF12_direct_-/-_(136r)</td>
      <td>5.374822</td>
      <td>-1</td>
      <td>-0.474900</td>
      <td>791</td>
    </tr>
    <tr>
      <th>7649</th>
      <td>chr18:36376597-36377097</td>
      <td>FHOD3</td>
      <td>0.006648</td>
      <td>-0.060058</td>
      <td>-0.000399</td>
      <td>0.000399</td>
      <td>TCF12</td>
      <td>False</td>
      <td>TCF12_direct_-/-</td>
      <td>TCF12_direct_-/-_(123g)</td>
      <td>TCF12_direct_-/-_(136r)</td>
      <td>1.216466</td>
      <td>-1</td>
      <td>-0.411142</td>
      <td>4699</td>
    </tr>
    <tr>
      <th>7650</th>
      <td>chr14:77307021-77307521</td>
      <td>TMED8</td>
      <td>0.011415</td>
      <td>-0.258848</td>
      <td>-0.002955</td>
      <td>0.002955</td>
      <td>TCF12</td>
      <td>False</td>
      <td>TCF12_direct_-/-</td>
      <td>TCF12_direct_-/-_(123g)</td>
      <td>TCF12_direct_-/-_(136r)</td>
      <td>1.879806</td>
      <td>-1</td>
      <td>-0.217667</td>
      <td>6561</td>
    </tr>
    <tr>
      <th>7651</th>
      <td>chr14:103019978-103020478</td>
      <td>TRAF3</td>
      <td>0.035434</td>
      <td>-0.251028</td>
      <td>-0.008895</td>
      <td>0.008895</td>
      <td>TCF12</td>
      <td>False</td>
      <td>TCF12_direct_-/-</td>
      <td>TCF12_direct_-/-_(123g)</td>
      <td>TCF12_direct_-/-_(136r)</td>
      <td>1.481711</td>
      <td>-1</td>
      <td>-0.183403</td>
      <td>5162</td>
    </tr>
    <tr>
      <th>7652</th>
      <td>chr12:20587950-20588450</td>
      <td>PDE3A</td>
      <td>0.018121</td>
      <td>-0.314249</td>
      <td>-0.005695</td>
      <td>0.005695</td>
      <td>TCF12</td>
      <td>False</td>
      <td>TCF12_direct_-/-</td>
      <td>TCF12_direct_-/-_(123g)</td>
      <td>TCF12_direct_-/-_(136r)</td>
      <td>2.102950</td>
      <td>-1</td>
      <td>-0.393440</td>
      <td>2064</td>
    </tr>
  </tbody>
</table>
<p>7653 rows × 15 columns</p>
</div>




```python
scplus_mdata.uns["extended_e_regulon_metadata"]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Region</th>
      <th>Gene</th>
      <th>importance_R2G</th>
      <th>rho_R2G</th>
      <th>importance_x_rho</th>
      <th>importance_x_abs_rho</th>
      <th>TF</th>
      <th>is_extended</th>
      <th>eRegulon_name</th>
      <th>Gene_signature_name</th>
      <th>Region_signature_name</th>
      <th>importance_TF2G</th>
      <th>regulation</th>
      <th>rho_TF2G</th>
      <th>triplet_rank</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1:6475533-6476033</td>
      <td>TNFRSF25</td>
      <td>0.019441</td>
      <td>0.217651</td>
      <td>0.004231</td>
      <td>0.004231</td>
      <td>EGR3</td>
      <td>True</td>
      <td>EGR3_extended_+/+</td>
      <td>EGR3_extended_+/+_(24g)</td>
      <td>EGR3_extended_+/+_(30r)</td>
      <td>1.070784</td>
      <td>1</td>
      <td>0.231038</td>
      <td>4058</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr2:241803041-241803541</td>
      <td>ATG4B</td>
      <td>0.015792</td>
      <td>0.066621</td>
      <td>0.001052</td>
      <td>0.001052</td>
      <td>EGR3</td>
      <td>True</td>
      <td>EGR3_extended_+/+</td>
      <td>EGR3_extended_+/+_(24g)</td>
      <td>EGR3_extended_+/+_(30r)</td>
      <td>1.234643</td>
      <td>1</td>
      <td>0.074444</td>
      <td>5455</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr2:42101928-42102428</td>
      <td>EML4</td>
      <td>0.037126</td>
      <td>0.276649</td>
      <td>0.010271</td>
      <td>0.010271</td>
      <td>EGR3</td>
      <td>True</td>
      <td>EGR3_extended_+/+</td>
      <td>EGR3_extended_+/+_(24g)</td>
      <td>EGR3_extended_+/+_(30r)</td>
      <td>1.235069</td>
      <td>1</td>
      <td>0.202865</td>
      <td>2496</td>
    </tr>
    <tr>
      <th>3</th>
      <td>chr10:114303963-114304463</td>
      <td>TDRD1</td>
      <td>0.082349</td>
      <td>0.078376</td>
      <td>0.006454</td>
      <td>0.006454</td>
      <td>EGR3</td>
      <td>True</td>
      <td>EGR3_extended_+/+</td>
      <td>EGR3_extended_+/+_(24g)</td>
      <td>EGR3_extended_+/+_(30r)</td>
      <td>0.865910</td>
      <td>1</td>
      <td>0.071393</td>
      <td>1707</td>
    </tr>
    <tr>
      <th>4</th>
      <td>chr19:18161129-18161629</td>
      <td>IFI30</td>
      <td>0.027840</td>
      <td>0.106714</td>
      <td>0.002971</td>
      <td>0.002971</td>
      <td>EGR3</td>
      <td>True</td>
      <td>EGR3_extended_+/+</td>
      <td>EGR3_extended_+/+_(24g)</td>
      <td>EGR3_extended_+/+_(30r)</td>
      <td>0.962073</td>
      <td>1</td>
      <td>0.354533</td>
      <td>3209</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>6137</th>
      <td>chr22:25436989-25437489</td>
      <td>GRK3</td>
      <td>0.003550</td>
      <td>-0.114036</td>
      <td>-0.000405</td>
      <td>0.000405</td>
      <td>TCF12</td>
      <td>True</td>
      <td>TCF12_extended_-/-</td>
      <td>TCF12_extended_-/-_(174g)</td>
      <td>TCF12_extended_-/-_(202r)</td>
      <td>1.927643</td>
      <td>-1</td>
      <td>-0.302335</td>
      <td>4648</td>
    </tr>
    <tr>
      <th>6138</th>
      <td>chr12:104340926-104341426</td>
      <td>CHST11</td>
      <td>0.001517</td>
      <td>-0.152301</td>
      <td>-0.000231</td>
      <td>0.000231</td>
      <td>TCF12</td>
      <td>True</td>
      <td>TCF12_extended_-/-</td>
      <td>TCF12_extended_-/-_(174g)</td>
      <td>TCF12_extended_-/-_(202r)</td>
      <td>1.120303</td>
      <td>-1</td>
      <td>-0.326319</td>
      <td>4314</td>
    </tr>
    <tr>
      <th>6139</th>
      <td>chr20:41184836-41185336</td>
      <td>ZHX3</td>
      <td>0.001790</td>
      <td>-0.133616</td>
      <td>-0.000239</td>
      <td>0.000239</td>
      <td>TCF12</td>
      <td>True</td>
      <td>TCF12_extended_-/-</td>
      <td>TCF12_extended_-/-_(174g)</td>
      <td>TCF12_extended_-/-_(202r)</td>
      <td>1.006590</td>
      <td>-1</td>
      <td>-0.296368</td>
      <td>3888</td>
    </tr>
    <tr>
      <th>6140</th>
      <td>chr5:66827771-66828271</td>
      <td>MAST4</td>
      <td>0.002882</td>
      <td>-0.218997</td>
      <td>-0.000631</td>
      <td>0.000631</td>
      <td>TCF12</td>
      <td>True</td>
      <td>TCF12_extended_-/-</td>
      <td>TCF12_extended_-/-_(174g)</td>
      <td>TCF12_extended_-/-_(202r)</td>
      <td>2.580701</td>
      <td>-1</td>
      <td>-0.420271</td>
      <td>2502</td>
    </tr>
    <tr>
      <th>6141</th>
      <td>chr15:34337164-34337664</td>
      <td>GOLGA8A</td>
      <td>0.026024</td>
      <td>-0.375745</td>
      <td>-0.009779</td>
      <td>0.009779</td>
      <td>TCF12</td>
      <td>True</td>
      <td>TCF12_extended_-/-</td>
      <td>TCF12_extended_-/-_(174g)</td>
      <td>TCF12_extended_-/-_(202r)</td>
      <td>2.319901</td>
      <td>-1</td>
      <td>-0.286801</td>
      <td>2936</td>
    </tr>
  </tbody>
</table>
<p>6142 rows × 15 columns</p>
</div>



## Downstream analysis

## eRegulon dimensionality reduction

The eRegulon enrichment scores can be used to perform dimensionality reductions


```python
import scanpy as sc 
import anndata
eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,
)
```


```python
eRegulon_gene_AUC.obs = scplus_mdata.obs.loc[eRegulon_gene_AUC.obs_names]
```


```python
sc.pp.neighbors(eRegulon_gene_AUC, use_rep = "X")
```


```python
sc.tl.umap(eRegulon_gene_AUC)
```


```python
sc.pl.umap(eRegulon_gene_AUC, color = "scRNA_counts:Seurat_cell_type")
```

    
![png](output_29_1.png)
    


## eRegulon specificity score


```python
from scenicplus.RSS import (regulon_specificity_scores, plot_rss)
```


```python
rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scRNA_counts:Seurat_cell_type",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC"]
)
```


```python
plot_rss(
    data_matrix = rss,
    top_n = 3,
    num_columns = 5
)
```


    
![png](output_33_0.png)
    


## Plot eRegulon enrichment scores

eRegulon enrichment scores can be plotted on the UMAP.


```python
sc.pl.umap(eRegulon_gene_AUC, color = list(set([x for xs in [rss.loc[ct].sort_values()[0:2].index for ct in rss.index] for x in xs ])))
```


    
![png](output_35_1.png)
    


## Heatmap dotplot

We can draw a heatmap where the color represent target gene enrichment and the dotsize target region enrichment.


```python
from scenicplus.plotting.dotplot import heatmap_dotplot
```


```python
heatmap_dotplot(
    scplus_mudata = scplus_mdata,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scRNA_counts:Seurat_cell_type",
    eRegulon_metadata_key = "direct_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    orientation = "horizontal",
    figsize = (16, 5)
)
```


    
![png](output_38_0.png)
    





    <Figure Size: (1600 x 500)>



# Converting mudata output to old-style SCENIC+ object.

Not all functions in the original release of SCENIC+ are updated to use the new mudata output of SCENIC+. To be able to still use these old functions while they get updated we have a function to convert the mudata object the the old SCENIC+ object.


```python
import os
os.chdir("/staging/leuven/stg_00002/lcb/sdewin/PhD/python_modules/scenicplus_development_tutorial")
```


```python
import mudata
scplus_mdata = mudata.read("outs/scplusmdata.h5mu")
```

```python
from scenicplus.scenicplus_class import mudata_to_scenicplus
```

To regenerate the complete SCENIC+ object it is necessary to provide the paths to your motif enrichment results. However, this data is not necessary for most downstream functions so providing these paths is completely **optional**.


```python
scplus_obj = mudata_to_scenicplus(
    mdata = scplus_mdata,
    path_to_cistarget_h5 = "outs/ctx_results.hdf5",
    path_to_dem_h5 = "outs/dem_results.hdf5"
)
```


```python
scplus_obj
```

