# Examine overlaps with published data


```python
import numpy as np
import pandas as pd
import session_info
from pyhere import here
```

## Public si-eQTL analysis


```python
shen = ["GDAP2", "AIM2", "SLAMF6", "RLF", "ATG4C", "FUT7",
        "TMEM218", "C11orf74", "RAB35", "TMEM5", "HNRNPK",
        "CDCA3", "ERCC5", "GJB6", "SNTB2", "SPNS3", 
        "XAF1", "RBBP8", "RUFY4", "CA2", "RAPGEF1"]
print("Shen et al.:")
print(len(shen))

kukurba = ["NOD2", "WDR36", "BSCL2", "MAP7D3", "RHOXF1", "DNAH1"]
print("Kukurba et al.:")
print(len(kukurba))

yao = ["NOD2", "HLA-DRB5", "HLA-DRB5", "KIAA0586", "PPP2R5A", 
       "TSNAXIP1", "MUT", "GRIK2", "C15orf37", "LIMA1", "IL6ST", 
       "HCG8", "BLOC1S3", "NKX3-1", "CXorf23"]
print("Yao et al.:")
print(len(np.unique(yao)))
print("Total of Yao + Kukurba:")
len(set(yao) | set(kukurba))
```

    Shen et al.:
    21
    Kukurba et al.:
    6
    Yao et al.:
    14
    Total of Yao + Kukurba:





    19



## Load BrainSeq si-eQTL results

### Interacting variant-gene pairs


```python
bs0 = pd.read_csv("../../_m/BrainSeq_sexGenotypes_4features_3regions.txt.gz", sep='\t')
bs0["ensembl_gene_id"] = bs0.gene_id.str.replace("\\..*", "", regex=True)
biomart = pd.read_csv("../_h/biomart.csv", index_col=0)
bs = bs0.merge(biomart, on="ensembl_gene_id").drop_duplicates(subset="gene_id")
print(bs.shape)
bs.tail(2)
```

    (692, 15)





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
      <th>region</th>
      <th>gene_id</th>
      <th>variant_id</th>
      <th>gencode_id</th>
      <th>gene_name</th>
      <th>seqnames</th>
      <th>start</th>
      <th>end</th>
      <th>lfsr</th>
      <th>posterior_mean</th>
      <th>feature_type</th>
      <th>ensembl_gene_id</th>
      <th>external_gene_name</th>
      <th>entrezgene</th>
      <th>description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>8285</th>
      <td>Caudate</td>
      <td>ENSG00000270605.1</td>
      <td>chr1:28102893:G:C</td>
      <td>ENSG00000270605.1</td>
      <td>ENSG00000270605</td>
      <td>chr1</td>
      <td>28239509</td>
      <td>28241453</td>
      <td>0.049873</td>
      <td>-0.261352</td>
      <td>Gene</td>
      <td>ENSG00000270605</td>
      <td>AL353622.1</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>8286</th>
      <td>DLPFC</td>
      <td>ENSG00000187498.16</td>
      <td>chr13:109650494:C:T</td>
      <td>ENSG00000187498.16</td>
      <td>COL4A1</td>
      <td>chr13</td>
      <td>110148963</td>
      <td>110307157</td>
      <td>0.048863</td>
      <td>0.193807</td>
      <td>Gene</td>
      <td>ENSG00000187498</td>
      <td>COL4A1</td>
      <td>1282.0</td>
      <td>collagen type IV alpha 1 chain [Source:HGNC Sy...</td>
    </tr>
  </tbody>
</table>
</div>




```python
bs[(bs['external_gene_name'].isin(shen))].to_csv("siEQTL_Shen_comparison.csv", index=False)
```


```python
bs[(bs['external_gene_name'].isin(kukurba))]
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
      <th>region</th>
      <th>gene_id</th>
      <th>variant_id</th>
      <th>gencode_id</th>
      <th>gene_name</th>
      <th>seqnames</th>
      <th>start</th>
      <th>end</th>
      <th>lfsr</th>
      <th>posterior_mean</th>
      <th>feature_type</th>
      <th>ensembl_gene_id</th>
      <th>external_gene_name</th>
      <th>entrezgene</th>
      <th>description</th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>
</div>




```python
bs[(bs['external_gene_name'].isin(yao))]
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
      <th>region</th>
      <th>gene_id</th>
      <th>variant_id</th>
      <th>gencode_id</th>
      <th>gene_name</th>
      <th>seqnames</th>
      <th>start</th>
      <th>end</th>
      <th>lfsr</th>
      <th>posterior_mean</th>
      <th>feature_type</th>
      <th>ensembl_gene_id</th>
      <th>external_gene_name</th>
      <th>entrezgene</th>
      <th>description</th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>
</div>




```python
bs[(bs['external_gene_name'].isin(shen+kukurba+yao))]
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
      <th>region</th>
      <th>gene_id</th>
      <th>variant_id</th>
      <th>gencode_id</th>
      <th>gene_name</th>
      <th>seqnames</th>
      <th>start</th>
      <th>end</th>
      <th>lfsr</th>
      <th>posterior_mean</th>
      <th>feature_type</th>
      <th>ensembl_gene_id</th>
      <th>external_gene_name</th>
      <th>entrezgene</th>
      <th>description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>5787</th>
      <td>Caudate</td>
      <td>ENSG00000125703.15</td>
      <td>chr1:63060301:G:A</td>
      <td>ENSG00000125703.15</td>
      <td>ATG4C</td>
      <td>chr1</td>
      <td>62784132</td>
      <td>62865516</td>
      <td>0.019839</td>
      <td>0.094136</td>
      <td>Gene</td>
      <td>ENSG00000125703</td>
      <td>ATG4C</td>
      <td>84938.0</td>
      <td>autophagy related 4C cysteine peptidase [Sourc...</td>
    </tr>
    <tr>
      <th>8130</th>
      <td>Caudate</td>
      <td>ENSG00000104267.10</td>
      <td>chr8:84966439:A:T</td>
      <td>ENSG00000104267.10</td>
      <td>CA2</td>
      <td>chr8</td>
      <td>85463968</td>
      <td>85481493</td>
      <td>0.045578</td>
      <td>0.289557</td>
      <td>Gene</td>
      <td>ENSG00000104267</td>
      <td>CA2</td>
      <td>760.0</td>
      <td>carbonic anhydrase 2 [Source:HGNC Symbol;Acc:H...</td>
    </tr>
  </tbody>
</table>
</div>



## GTEx comparison


```python
gtex = pd.read_csv(here("input/public_results/gtex_results/_m",
                        "GTEx_Analysis_v8_sbeQTLs/GTEx_Analysis_v8_sbeQTLs.txt"), 
                   sep='\t')
gtex.iloc[0:2, 0:10]
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
      <th>ensembl_gene_id</th>
      <th>hugo_gene_id</th>
      <th>gene_type</th>
      <th>variant_id</th>
      <th>rs_id</th>
      <th>Tissue</th>
      <th>maf</th>
      <th>pval_nominal_sb</th>
      <th>slope_sb</th>
      <th>slope_se_sb</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000241860.6</td>
      <td>RP11-34P13.13</td>
      <td>processed_transcript</td>
      <td>chr1_14677_G_A_b38</td>
      <td>rs201327123</td>
      <td>Adipose_Subcutaneous</td>
      <td>0.051635</td>
      <td>0.847114</td>
      <td>0.055080</td>
      <td>0.285537</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000227232.5</td>
      <td>WASH7P</td>
      <td>unprocessed_pseudogene</td>
      <td>chr1_64764_C_T_b38</td>
      <td>rs769952832</td>
      <td>Adipose_Subcutaneous</td>
      <td>0.061102</td>
      <td>0.316881</td>
      <td>0.222928</td>
      <td>0.222511</td>
    </tr>
  </tbody>
</table>
</div>




```python
gtex.iloc[0:2, 10:14]
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
      <th>numtested</th>
      <th>pvals.corrected</th>
      <th>qval</th>
      <th>pval_nominal_f</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>0.847114</td>
      <td>1.000000</td>
      <td>0.022302</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>0.316881</td>
      <td>0.981254</td>
      <td>0.003978</td>
    </tr>
  </tbody>
</table>
</div>




```python
## qval threshold equal to number of published sb-eQTL
gtex[(gtex['qval'] < 0.25) & (gtex["Tissue"].str.contains("Brain"))]\
    .loc[:, ["ensembl_gene_id", "hugo_gene_id", "Tissue", "pvals.corrected", 'qval']].head(10)
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
      <th>ensembl_gene_id</th>
      <th>hugo_gene_id</th>
      <th>Tissue</th>
      <th>pvals.corrected</th>
      <th>qval</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>62155</th>
      <td>ENSG00000026025.15</td>
      <td>VIM</td>
      <td>Brain_Amygdala</td>
      <td>0.000004</td>
      <td>0.012836</td>
    </tr>
    <tr>
      <th>116842</th>
      <td>ENSG00000160818.16</td>
      <td>GPATCH4</td>
      <td>Brain_Nucleus_accumbens_basal_ganglia</td>
      <td>0.000088</td>
      <td>0.198445</td>
    </tr>
    <tr>
      <th>121904</th>
      <td>ENSG00000141562.17</td>
      <td>NARF</td>
      <td>Brain_Nucleus_accumbens_basal_ganglia</td>
      <td>0.000056</td>
      <td>0.198445</td>
    </tr>
    <tr>
      <th>122123</th>
      <td>ENSG00000267174.5</td>
      <td>CTC-510F12.4</td>
      <td>Brain_Nucleus_accumbens_basal_ganglia</td>
      <td>0.000083</td>
      <td>0.198445</td>
    </tr>
  </tbody>
</table>
</div>




```python
## qval threshold equal to number of published sb-eQTL
gtex[(gtex['qval'] < 0.25) & (gtex["Tissue"].str.contains("Whole"))]\
    .loc[:, ["ensembl_gene_id", "hugo_gene_id", "Tissue", "pvals.corrected", 'qval']].head(10)
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
      <th>ensembl_gene_id</th>
      <th>hugo_gene_id</th>
      <th>Tissue</th>
      <th>pvals.corrected</th>
      <th>qval</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>362961</th>
      <td>ENSG00000221571.3</td>
      <td>RNU6ATAC35P</td>
      <td>Whole_Blood</td>
      <td>0.000039</td>
      <td>0.139762</td>
    </tr>
    <tr>
      <th>365043</th>
      <td>ENSG00000196743.8</td>
      <td>GM2A</td>
      <td>Whole_Blood</td>
      <td>0.000011</td>
      <td>0.116825</td>
    </tr>
    <tr>
      <th>367164</th>
      <td>ENSG00000148459.15</td>
      <td>PDSS1</td>
      <td>Whole_Blood</td>
      <td>0.000027</td>
      <td>0.139762</td>
    </tr>
  </tbody>
</table>
</div>




```python
gtex_sig = gtex[(gtex['qval'] < 0.25)]
gtex_sig.shape
```




    (369, 22)




```python
gtex_sig.head(10)
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
      <th>ensembl_gene_id</th>
      <th>hugo_gene_id</th>
      <th>gene_type</th>
      <th>variant_id</th>
      <th>rs_id</th>
      <th>Tissue</th>
      <th>maf</th>
      <th>pval_nominal_sb</th>
      <th>slope_sb</th>
      <th>slope_se_sb</th>
      <th>...</th>
      <th>qval</th>
      <th>pval_nominal_f</th>
      <th>slope_f</th>
      <th>slope_se_f</th>
      <th>pval_nominal_m</th>
      <th>slope_m</th>
      <th>slope_se_m</th>
      <th>pval_nominal</th>
      <th>slope</th>
      <th>slope_se</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1096</th>
      <td>ENSG00000076356.6</td>
      <td>PLXNA2</td>
      <td>protein_coding</td>
      <td>chr1_208030492_G_A_b38</td>
      <td>rs3811383</td>
      <td>Adipose_Subcutaneous</td>
      <td>0.123924</td>
      <td>5.391600e-05</td>
      <td>0.338278</td>
      <td>0.083064</td>
      <td>...</td>
      <td>0.121068</td>
      <td>1.718880e-08</td>
      <td>0.456729</td>
      <td>0.075705</td>
      <td>9.155700e-01</td>
      <td>0.009739</td>
      <td>0.091682</td>
      <td>2.747400e-05</td>
      <td>0.171830</td>
      <td>0.040604</td>
    </tr>
    <tr>
      <th>5262</th>
      <td>ENSG00000170632.13</td>
      <td>ARMC10</td>
      <td>protein_coding</td>
      <td>chr7_103076937_C_T_b38</td>
      <td>rs6958836</td>
      <td>Adipose_Subcutaneous</td>
      <td>0.169535</td>
      <td>5.011130e-05</td>
      <td>0.357403</td>
      <td>0.087384</td>
      <td>...</td>
      <td>0.192900</td>
      <td>4.933240e-01</td>
      <td>-0.054539</td>
      <td>0.079379</td>
      <td>3.219220e-07</td>
      <td>-0.429800</td>
      <td>0.079545</td>
      <td>8.797530e-08</td>
      <td>-0.216374</td>
      <td>0.039857</td>
    </tr>
    <tr>
      <th>5644</th>
      <td>ENSG00000120907.17</td>
      <td>ADRA1A</td>
      <td>protein_coding</td>
      <td>chr8_26839198_G_A_b38</td>
      <td>rs117380715</td>
      <td>Adipose_Subcutaneous</td>
      <td>0.216867</td>
      <td>1.045890e-05</td>
      <td>-0.323552</td>
      <td>0.072676</td>
      <td>...</td>
      <td>0.084548</td>
      <td>4.637410e-18</td>
      <td>-0.779707</td>
      <td>0.076596</td>
      <td>3.976660e-10</td>
      <td>-0.469672</td>
      <td>0.069091</td>
      <td>5.637370e-52</td>
      <td>-0.568916</td>
      <td>0.033334</td>
    </tr>
    <tr>
      <th>6414</th>
      <td>ENSG00000136830.11</td>
      <td>FAM129B</td>
      <td>protein_coding</td>
      <td>chr9_127584339_G_A_b38</td>
      <td>rs10739693</td>
      <td>Adipose_Subcutaneous</td>
      <td>0.304647</td>
      <td>7.387010e-07</td>
      <td>-0.283660</td>
      <td>0.056579</td>
      <td>...</td>
      <td>0.004976</td>
      <td>1.978000e-06</td>
      <td>-0.333315</td>
      <td>0.066772</td>
      <td>1.653380e-01</td>
      <td>-0.082625</td>
      <td>0.059205</td>
      <td>1.393160e-08</td>
      <td>-0.168762</td>
      <td>0.029260</td>
    </tr>
    <tr>
      <th>7220</th>
      <td>ENSG00000166787.3</td>
      <td>SAA3P</td>
      <td>transcribed_unprocessed_pseudogene</td>
      <td>chr11_18269355_T_C_b38</td>
      <td>rs34068567</td>
      <td>Adipose_Subcutaneous</td>
      <td>0.278830</td>
      <td>2.207290e-05</td>
      <td>0.323030</td>
      <td>0.075427</td>
      <td>...</td>
      <td>0.074347</td>
      <td>6.409400e-08</td>
      <td>0.453034</td>
      <td>0.078725</td>
      <td>3.138000e-01</td>
      <td>0.063002</td>
      <td>0.062292</td>
      <td>2.433600e-08</td>
      <td>0.211910</td>
      <td>0.037395</td>
    </tr>
    <tr>
      <th>8540</th>
      <td>ENSG00000183463.5</td>
      <td>URAD</td>
      <td>protein_coding</td>
      <td>chr13_27990205_T_A_b38</td>
      <td>rs7335293</td>
      <td>Adipose_Subcutaneous</td>
      <td>0.500000</td>
      <td>9.078700e-09</td>
      <td>-0.444892</td>
      <td>0.076123</td>
      <td>...</td>
      <td>0.000122</td>
      <td>9.982650e-21</td>
      <td>-0.887723</td>
      <td>0.078738</td>
      <td>1.892290e-09</td>
      <td>-0.457733</td>
      <td>0.070571</td>
      <td>3.077310e-53</td>
      <td>-0.640604</td>
      <td>0.036976</td>
    </tr>
    <tr>
      <th>9191</th>
      <td>ENSG00000282651.2</td>
      <td>IGHV5-10-1</td>
      <td>IG_V_gene</td>
      <td>chr14_106114510_A_G_b38</td>
      <td>rs4573838</td>
      <td>Adipose_Subcutaneous</td>
      <td>0.419105</td>
      <td>2.025150e-05</td>
      <td>-0.406760</td>
      <td>0.094541</td>
      <td>...</td>
      <td>0.074347</td>
      <td>5.248710e-12</td>
      <td>-0.682629</td>
      <td>0.089412</td>
      <td>2.805060e-03</td>
      <td>-0.289091</td>
      <td>0.094806</td>
      <td>3.458420e-21</td>
      <td>-0.445408</td>
      <td>0.045073</td>
    </tr>
    <tr>
      <th>14611</th>
      <td>ENSG00000143933.16</td>
      <td>CALM2</td>
      <td>protein_coding</td>
      <td>chr2_46225349_C_T_b38</td>
      <td>rs12477148</td>
      <td>Adipose_Visceral_Omentum</td>
      <td>0.072495</td>
      <td>4.497930e-05</td>
      <td>-0.480557</td>
      <td>0.116471</td>
      <td>...</td>
      <td>0.161955</td>
      <td>4.747150e-04</td>
      <td>-0.491287</td>
      <td>0.134732</td>
      <td>9.165740e-01</td>
      <td>0.013043</td>
      <td>0.124116</td>
      <td>2.197750e-05</td>
      <td>-0.246023</td>
      <td>0.057281</td>
    </tr>
    <tr>
      <th>15082</th>
      <td>ENSG00000144410.4</td>
      <td>CPO</td>
      <td>protein_coding</td>
      <td>chr2_206822186_C_T_b38</td>
      <td>rs12470278</td>
      <td>Adipose_Visceral_Omentum</td>
      <td>0.097015</td>
      <td>3.204120e-05</td>
      <td>0.682291</td>
      <td>0.162191</td>
      <td>...</td>
      <td>0.115370</td>
      <td>1.165430e-01</td>
      <td>0.280837</td>
      <td>0.176978</td>
      <td>4.431060e-06</td>
      <td>-0.558002</td>
      <td>0.113158</td>
      <td>7.896400e-06</td>
      <td>-0.320288</td>
      <td>0.070745</td>
    </tr>
    <tr>
      <th>17452</th>
      <td>ENSG00000211698.2</td>
      <td>TRGV4</td>
      <td>TR_V_gene</td>
      <td>chr7_38361995_A_C_b38</td>
      <td>rs10233345</td>
      <td>Adipose_Visceral_Omentum</td>
      <td>0.335821</td>
      <td>6.438100e-05</td>
      <td>0.427491</td>
      <td>0.105837</td>
      <td>...</td>
      <td>0.139089</td>
      <td>7.011230e-06</td>
      <td>-0.481758</td>
      <td>0.100091</td>
      <td>8.857670e-15</td>
      <td>-1.068840</td>
      <td>0.112111</td>
      <td>1.630390e-49</td>
      <td>-0.838766</td>
      <td>0.049090</td>
    </tr>
  </tbody>
</table>
<p>10 rows × 22 columns</p>
</div>



### mashr


```python
gtex_overlap = bs[(bs['gene_id'].isin(gtex_sig.ensembl_gene_id))].drop_duplicates()
print(gtex_overlap.shape)
gtex_overlap
```

    (2, 15)





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
      <th>region</th>
      <th>gene_id</th>
      <th>variant_id</th>
      <th>gencode_id</th>
      <th>gene_name</th>
      <th>seqnames</th>
      <th>start</th>
      <th>end</th>
      <th>lfsr</th>
      <th>posterior_mean</th>
      <th>feature_type</th>
      <th>ensembl_gene_id</th>
      <th>external_gene_name</th>
      <th>entrezgene</th>
      <th>description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>4943</th>
      <td>Caudate</td>
      <td>ENSG00000272977.1</td>
      <td>chr22:25059120:A:C</td>
      <td>ENSG00000272977.1</td>
      <td>ENSG00000272977</td>
      <td>chr22</td>
      <td>25476218</td>
      <td>25479971</td>
      <td>0.011928</td>
      <td>0.323847</td>
      <td>Gene</td>
      <td>ENSG00000272977</td>
      <td>AL008721.2</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>8285</th>
      <td>Caudate</td>
      <td>ENSG00000270605.1</td>
      <td>chr1:28102893:G:C</td>
      <td>ENSG00000270605.1</td>
      <td>ENSG00000270605</td>
      <td>chr1</td>
      <td>28239509</td>
      <td>28241453</td>
      <td>0.049873</td>
      <td>-0.261352</td>
      <td>Gene</td>
      <td>ENSG00000270605</td>
      <td>AL353622.1</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>




```python
gtex_overlap.shape[0]/bs.shape[0] * 100
```




    0.2890173410404624




```python
gtex_sig[(gtex_sig['ensembl_gene_id'].isin(bs.gene_id))]
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
      <th>ensembl_gene_id</th>
      <th>hugo_gene_id</th>
      <th>gene_type</th>
      <th>variant_id</th>
      <th>rs_id</th>
      <th>Tissue</th>
      <th>maf</th>
      <th>pval_nominal_sb</th>
      <th>slope_sb</th>
      <th>slope_se_sb</th>
      <th>...</th>
      <th>qval</th>
      <th>pval_nominal_f</th>
      <th>slope_f</th>
      <th>slope_se_f</th>
      <th>pval_nominal_m</th>
      <th>slope_m</th>
      <th>slope_se_m</th>
      <th>pval_nominal</th>
      <th>slope</th>
      <th>slope_se</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>297207</th>
      <td>ENSG00000270605.1</td>
      <td>RP5-1092A3.4</td>
      <td>antisense</td>
      <td>chr1_28223937_C_T_b38</td>
      <td>rs481640</td>
      <td>Skin_Not_Sun_Exposed_Suprapubic</td>
      <td>0.323985</td>
      <td>0.000237</td>
      <td>-0.241524</td>
      <td>0.065161</td>
      <td>...</td>
      <td>0.208195</td>
      <td>4.046590e-21</td>
      <td>-0.760755</td>
      <td>0.063142</td>
      <td>6.395880e-13</td>
      <td>-0.599143</td>
      <td>0.07247</td>
      <td>9.917080e-60</td>
      <td>-0.591873</td>
      <td>0.031047</td>
    </tr>
    <tr>
      <th>338770</th>
      <td>ENSG00000272977.1</td>
      <td>CTA-390C10.10</td>
      <td>sense_intronic</td>
      <td>chr22_25459662_G_A_b38</td>
      <td>rs6004655</td>
      <td>Spleen</td>
      <td>0.167401</td>
      <td>0.000038</td>
      <td>0.413128</td>
      <td>0.097856</td>
      <td>...</td>
      <td>0.212883</td>
      <td>7.931210e-14</td>
      <td>-1.129030</td>
      <td>0.107327</td>
      <td>5.155260e-16</td>
      <td>-1.276590</td>
      <td>0.10469</td>
      <td>4.896640e-51</td>
      <td>-1.193670</td>
      <td>0.056659</td>
    </tr>
  </tbody>
</table>
<p>2 rows × 22 columns</p>
</div>




```python
gtex_sig[(gtex_sig['ensembl_gene_id'].isin(bs.gene_id))]\
    .to_csv("siEQTL_gtex_comparison.csv", index=False)
```

# Session information


```python
session_info.show()
```




<details>
<summary>Click to view session information</summary>
<pre>
-----
numpy               1.20.2
pandas              1.4.3
pyhere              1.0.0
session_info        1.0.0
-----
</pre>
<details>
<summary>Click to view modules imported as dependencies</summary>
<pre>
anyio               NA
arrow               1.2.3
asttokens           NA
attr                21.4.0
babel               2.10.3
backcall            0.2.0
certifi             2021.10.08
chardet             4.0.0
charset_normalizer  2.0.12
cloudpickle         2.1.0
cython_runtime      NA
dateutil            2.8.2
debugpy             1.6.3
decorator           5.1.1
executing           0.10.0
fastjsonschema      NA
fqdn                NA
google              NA
idna                2.10
importlib_metadata  NA
ipykernel           6.15.1
ipython_genutils    0.2.0
isoduration         NA
jedi                0.18.1
jinja2              3.1.1
json5               NA
jsonpointer         2.3
jsonschema          4.17.3
jupyter_server      1.23.5
jupyterlab_server   2.19.0
markupsafe          2.1.1
mpl_toolkits        NA
nbformat            5.4.0
packaging           21.3
parso               0.8.3
pexpect             4.8.0
pickleshare         0.7.5
pkg_resources       NA
platformdirs        3.0.0
prometheus_client   NA
prompt_toolkit      3.0.30
psutil              5.9.1
ptyprocess          0.7.0
pure_eval           0.2.2
pvectorc            NA
pyarrow             7.0.0
pydev_ipython       NA
pydevconsole        NA
pydevd              2.8.0
pydevd_file_utils   NA
pydevd_plugins      NA
pydevd_tracing      NA
pygments            2.13.0
pyrsistent          NA
pytz                2022.1
requests            2.28.2
rfc3339_validator   0.1.4
rfc3986_validator   0.1.1
send2trash          NA
six                 1.16.0
sniffio             1.2.0
stack_data          0.4.0
terminado           0.15.0
tornado             6.2
traitlets           5.9.0
typing_extensions   NA
uri_template        NA
urllib3             1.26.8
wcwidth             0.2.5
webcolors           1.11.1
websocket           1.4.0
zipp                NA
zmq                 25.0.0
zstandard           0.18.0
</pre>
</details> <!-- seems like this ends pre, so might as well be explicit -->
<pre>
-----
IPython             8.4.0
jupyter_client      8.0.2
jupyter_core        5.2.0
jupyterlab          3.6.1
notebook            6.4.12
-----
Python 3.9.10 (main, Feb 22 2022, 16:34:24) [GCC 4.8.5 20150623 (Red Hat 4.8.5-44)]
Linux-3.10.0-1160.el7.x86_64-x86_64-with-glibc2.17
-----
Session information updated at 2023-06-08 09:10
</pre>
</details>


