# Preprocessing the scRNA-seq data

In this tutorial we will perform some very basic preprocessing steps.

## Download data

The data used for this tutorial is freely available, and can be downloaded from the [10x genomics website](https://www.10xgenomics.com/datasets/frozen-human-healthy-brain-tissue-3-k-1-standard-1-0-0).


```python
import os
os.chdir("/staging/leuven/stg_00002/lcb/sdewin/PhD/python_modules/scenicplus_development_tutorial/scRNA_seq_pp")
```


```python
!mkdir -p data
!wget -O data/filtered_feature_bc_matrix.tar.gz https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_filtered_feature_bc_matrix.tar.gz
```



```python
!cd data; tar -xzf filtered_feature_bc_matrix.tar.gz; cd ..
```


```python
!wget -O data/cell_data.tsv https://raw.githubusercontent.com/aertslab/pycisTopic/polars/data/cell_data_human_cerebellum.tsv
```


## Preprocessing

We will do some very basic preprocessing steps, for more information we refer the reader to the [Scanpy tutorials](https://scanpy.readthedocs.io/en/stable/tutorials.html).


```python
import scanpy as sc
```


```python
adata = sc.read_10x_mtx(
    "data/filtered_feature_bc_matrix/",
    var_names = "gene_symbols"
)
```


```python
adata.var_names_make_unique()
```


```python
adata
```


We have already annotated this dataset, this is beyond the scope of this tutorial. We will load this annotation here.


```python
import pandas as pd
cell_data = pd.read_table("data/cell_data.tsv", index_col = 0)
cell_data
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
      <th>VSN_cell_type</th>
      <th>VSN_leiden_res0.3</th>
      <th>VSN_leiden_res0.6</th>
      <th>VSN_leiden_res0.9</th>
      <th>VSN_leiden_res1.2</th>
      <th>VSN_sample_id</th>
      <th>Seurat_leiden_res0.6</th>
      <th>Seurat_leiden_res1.2</th>
      <th>Seurat_cell_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>AAACAGCCATTATGCG-1-10x_multiome_brain</th>
      <td>MOL_B</td>
      <td>MOL_B (0)</td>
      <td>MOL_B_1 (0)</td>
      <td>MOL_B_1  (1)</td>
      <td>MOL_B_3 (6)</td>
      <td>10x_multiome_brain</td>
      <td>NFOL (1)</td>
      <td>MOL (1)</td>
      <td>MOL</td>
    </tr>
    <tr>
      <th>AAACCAACATAGACCC-1-10x_multiome_brain</th>
      <td>MOL_B</td>
      <td>MOL_B (0)</td>
      <td>MOL_B_1 (0)</td>
      <td>MOL_B_3 (5)</td>
      <td>MOL_B_4 (4)</td>
      <td>10x_multiome_brain</td>
      <td>NFOL (1)</td>
      <td>NFOL (3)</td>
      <td>NFOL</td>
    </tr>
    <tr>
      <th>AAACCGAAGATGCCTG-1-10x_multiome_brain</th>
      <td>INH_VIP</td>
      <td>INH_VIP (6)</td>
      <td>INH_VIP (8)</td>
      <td>INH_VIP (8)</td>
      <td>INH_VIP (10)</td>
      <td>10x_multiome_brain</td>
      <td>INH_VIP (7)</td>
      <td>INH_VIP (6)</td>
      <td>INH_VIP</td>
    </tr>
    <tr>
      <th>AAACCGAAGTTAGCTA-1-10x_multiome_brain</th>
      <td>MOL_A</td>
      <td>MOL_A (1)</td>
      <td>MOL_A_2 (1)</td>
      <td>MOL_A_1  (0)</td>
      <td>MOL_A_2 (0)</td>
      <td>10x_multiome_brain</td>
      <td>NFOL (1)</td>
      <td>NFOL (3)</td>
      <td>NFOL</td>
    </tr>
    <tr>
      <th>AAACCGCGTTAGCCAA-1-10x_multiome_brain</th>
      <td>MGL</td>
      <td>MGL (7)</td>
      <td>MGL (10)</td>
      <td>MGL (10)</td>
      <td>MGL (12)</td>
      <td>10x_multiome_brain</td>
      <td>MGL (8)</td>
      <td>MGL (9)</td>
      <td>MGL</td>
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
    </tr>
    <tr>
      <th>TTTGTGAAGGGTGAGT-1-10x_multiome_brain</th>
      <td>INH_VIP</td>
      <td>INH_VIP (6)</td>
      <td>INH_VIP (8)</td>
      <td>INH_VIP (8)</td>
      <td>INH_VIP (10)</td>
      <td>10x_multiome_brain</td>
      <td>INH_SST (5)</td>
      <td>INH_SST (8)</td>
      <td>INH_SST</td>
    </tr>
    <tr>
      <th>TTTGTGAAGTCAGGCC-1-10x_multiome_brain</th>
      <td>AST_CER</td>
      <td>AST_CER (2)</td>
      <td>AST_CER (2)</td>
      <td>AST_CER (2)</td>
      <td>AST_CER_1 (7)</td>
      <td>10x_multiome_brain</td>
      <td>BG (2)</td>
      <td>BG (2)</td>
      <td>BG</td>
    </tr>
    <tr>
      <th>TTTGTGGCATGCTTAG-1-10x_multiome_brain</th>
      <td>MOL_B</td>
      <td>MOL_B (0)</td>
      <td>MOL_B_1 (0)</td>
      <td>MOL_B_1  (1)</td>
      <td>MOL_B_1 (1)</td>
      <td>10x_multiome_brain</td>
      <td>MOL (0)</td>
      <td>MOL (1)</td>
      <td>MOL</td>
    </tr>
    <tr>
      <th>TTTGTTGGTGATCAGC-1-10x_multiome_brain</th>
      <td>MOL_A</td>
      <td>MOL_A (1)</td>
      <td>MOL_A_2 (1)</td>
      <td>MOL_A_1  (0)</td>
      <td>MOL_A_1 (11)</td>
      <td>10x_multiome_brain</td>
      <td>NFOL (1)</td>
      <td>NFOL (3)</td>
      <td>NFOL</td>
    </tr>
    <tr>
      <th>TTTGTTGGTGATTTGG-1-10x_multiome_brain</th>
      <td>INH_SST</td>
      <td>INH_SST (5)</td>
      <td>INH_SST (7)</td>
      <td>INH_SST (7)</td>
      <td>INH_SST (9)</td>
      <td>10x_multiome_brain</td>
      <td>INH_SST (5)</td>
      <td>INH_SST (8)</td>
      <td>INH_SST</td>
    </tr>
  </tbody>
</table>
<p>2392 rows × 9 columns</p>
</div>



We modify the index of this cell type annotation dataframe so that the cell barcode names match with those in the AnnData object.


```python
cell_data.index = [cb.rsplit("-", 1)[0] for cb in cell_data.index]
```


```python
adata = adata[list(set(adata.obs_names) & set(cell_data.index))].copy()
```


```python
adata.obs = cell_data.loc[adata.obs_names]
```


```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
```

### Data normalization

It's important to save the **non normalized** and **non scaled** matrix in the raw slot!


```python
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
```


```python
adata.obs
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
      <th>VSN_cell_type</th>
      <th>VSN_leiden_res0.3</th>
      <th>VSN_leiden_res0.6</th>
      <th>VSN_leiden_res0.9</th>
      <th>VSN_leiden_res1.2</th>
      <th>VSN_sample_id</th>
      <th>Seurat_leiden_res0.6</th>
      <th>Seurat_leiden_res1.2</th>
      <th>Seurat_cell_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ACCTGTTGTGGATTCA-1</th>
      <td>MOL_A</td>
      <td>MOL_A (1)</td>
      <td>MOL_A_2 (1)</td>
      <td>MOL_A_1  (0)</td>
      <td>MOL_A_2 (0)</td>
      <td>10x_multiome_brain</td>
      <td>MOL (0)</td>
      <td>MOL (0)</td>
      <td>MOL</td>
    </tr>
    <tr>
      <th>CTAAAGCTCCCGCCTA-1</th>
      <td>INH_VIP</td>
      <td>INH_VIP (6)</td>
      <td>INH_VIP (8)</td>
      <td>INH_VIP (8)</td>
      <td>INH_VIP (10)</td>
      <td>10x_multiome_brain</td>
      <td>INH_VIP (7)</td>
      <td>INH_VIP (6)</td>
      <td>INH_VIP</td>
    </tr>
    <tr>
      <th>AAGGAAGCACATAACT-1</th>
      <td>AST_CER</td>
      <td>AST_CER (2)</td>
      <td>AST_CER (2)</td>
      <td>AST_CER (2)</td>
      <td>AST_CER_2 (5)</td>
      <td>10x_multiome_brain</td>
      <td>BG (2)</td>
      <td>BG (2)</td>
      <td>BG</td>
    </tr>
    <tr>
      <th>GAGCATGCAATATGGA-1</th>
      <td>MOL_A</td>
      <td>MOL_A (1)</td>
      <td>MOL_A_2 (1)</td>
      <td>MOL_A_1  (0)</td>
      <td>MOL_A_2 (0)</td>
      <td>10x_multiome_brain</td>
      <td>MOL (0)</td>
      <td>MOL (0)</td>
      <td>MOL</td>
    </tr>
    <tr>
      <th>GCTCACAAGACAAGTG-1</th>
      <td>MOL_A</td>
      <td>MOL_A (1)</td>
      <td>MOL_A_2 (1)</td>
      <td>MOL_A_1  (0)</td>
      <td>MOL_A_2 (0)</td>
      <td>10x_multiome_brain</td>
      <td>MOL (0)</td>
      <td>MOL (1)</td>
      <td>MOL</td>
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
    </tr>
    <tr>
      <th>ACATTGCAGCGGATTT-1</th>
      <td>MOL_B</td>
      <td>MOL_B (0)</td>
      <td>MOL_B_1 (0)</td>
      <td>MOL_B_1  (1)</td>
      <td>MOL_B_3 (6)</td>
      <td>10x_multiome_brain</td>
      <td>NFOL (1)</td>
      <td>NFOL (3)</td>
      <td>NFOL</td>
    </tr>
    <tr>
      <th>ACGCCTTTCGACCTGA-1</th>
      <td>AST</td>
      <td>AST+ENDO (9)</td>
      <td>AST+ENDO (6)</td>
      <td>AST+ENDO (13)</td>
      <td>AST (15)</td>
      <td>10x_multiome_brain</td>
      <td>AST+ENDO (6)</td>
      <td>AST (7)</td>
      <td>AST</td>
    </tr>
    <tr>
      <th>CTCAGGATCCACCTGT-1</th>
      <td>MOL_B</td>
      <td>MOL_B (0)</td>
      <td>MOL_B_1 (0)</td>
      <td>MOL_B_1  (1)</td>
      <td>MOL_B_1 (1)</td>
      <td>10x_multiome_brain</td>
      <td>NFOL (1)</td>
      <td>MOL (1)</td>
      <td>MOL</td>
    </tr>
    <tr>
      <th>CCAACATAGATAACCC-1</th>
      <td>MOL_A</td>
      <td>MOL_A (1)</td>
      <td>MOL_A_1 (9)</td>
      <td>MOL_A_1 (9)</td>
      <td>MOL_A_1 (11)</td>
      <td>10x_multiome_brain</td>
      <td>NFOL (1)</td>
      <td>COP (10)</td>
      <td>COP</td>
    </tr>
    <tr>
      <th>ATGTAAGCACTTAGGC-1</th>
      <td>MOL_B</td>
      <td>MOL_B (0)</td>
      <td>MOL_B_1 (0)</td>
      <td>MOL_B_1  (1)</td>
      <td>MOL_B_1 (1)</td>
      <td>10x_multiome_brain</td>
      <td>MOL (0)</td>
      <td>MOL (0)</td>
      <td>MOL</td>
    </tr>
  </tbody>
</table>
<p>2313 rows × 9 columns</p>
</div>




```python
sc.tl.pca(adata)
sc.pl.pca(adata, color = "Seurat_cell_type")
```


    
![png](output_21_1.png)
    



```python
sc.pp.neighbors(adata)
```


```python
sc.tl.umap(adata)
```


```python
sc.pl.umap(adata, color = "Seurat_cell_type")
```


    
![png](output_24_1.png)
    



```python
adata.write("adata.h5ad")
```
