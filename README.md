# San Pedro Ocean Time-series 18S-V4 Amplicon and Network Analyses 
**SPOT_Protist_Net**  
**By: Samantha J. Gleich**  
**Last updated: November 10, 2022**  


![](static/slide1.png)


# Make ASVs in Qiime2
First we will use our raw, 18S-V4 reads to create amplicon sequence variants (ASVs). This pipeline was adapted from a pipeline documented by Sarah Hu (https://github.com/shu251/V4_tagsequencing_18Sdiversity_qiime2)
## Make V4 classifier using the newest version of the PR2 database (version 14)
Import PR2 sequences as qiime artifact:
```
qiime tools import --type 'FeatureData[Sequence]' --input-path pr2_version_4.14.0_SSU_mothur.fasta --output-path pr2_v4_14.qza
```
Import PR2 taxonomy information as qiime artifact: 
```
qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path pr2_version_4.14.0_SSU_mothur.txt --output-path pr2_v4_14_tax.qza
```
Extract V4 reads from PR2 sequences:
```
qiime feature-classifier extract-reads --i-sequences pr2_v4_14.qza --p-f-primer CCAGCASCYGCGGTAATTCC --p-r-primer ACTTTCGTTCTTGATYRA --o-reads v4_extracts_v14.qza
```
Train the V4 classifier: 
```
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads v4_extracts_v14.qza --i-reference-taxonomy pr2_v4_14_tax.qza --o-classifier pr2_v4_v14_classifier.qza
```
## Trim reads using Trimmomatic (version 0.38)
Used Sarah Hu's run_trim.py script (https://github.com/shu251/V4_tagsequencing_18Sdiversity_qiime2)
```
python run_trim.py manifest.txt
```
## Import sequence files in as qiime artifact
```
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path trimmed_manifest.txt --output-path demux.qza --input-format PairedEndFastqManifestPhred33
```
## Remove 18S-V4 primers 
```
qiime cutadapt trim-paired --i-demultiplexed-sequences demux.qza --p-cores 8 --p-front-f CCAGCASCYGCGGTAATTCC --p-front-r ACTTTCGTTCTTGATYRA --o-trimmed-sequences demux_trimmed.qza
```
## Join paired-end sequences and remove merged sequences that are < 300 bp
```
qiime vsearch join-pairs --i-demultiplexed-seqs demux_trimmed.qza --o-joined-sequences demux-joined.qza --p-minmergelen 300
```
## Quality filter
```
qiime quality-filter q-score-joined --i-demux demux-joined.qza --o-filtered-sequences demux-joined-filtered.qza --o-filter-stats demux-joined-filter-stats.qza --p-min-quality 20
```
## Dereplicate sequences (i.e. remove dups)
```
qiime vsearch dereplicate-sequences --i-sequences demux-joined-filtered.qza --o-dereplicated-table derep_table.qza --o-dereplicated-sequences derep_seqs.qza
```
## Make ASVs
```
mkdir ASVs

qiime dada2 denoise-paired
	--i-demultiplexed-seqs demux_trimmed.qza \
	--o-table ASVs/table \
	--o-representative-sequences ASVs/rep-seqs \
	--p-n-threads 8 \
	--p-trunc-len-f 200 \
	--p-trunc-len-r 200 \
	--p-max-ee-f 2 \
	--p-max-ee-r 2 \
	--p-n-reads-learn 1000000 \
	--p-chimera-method pooled \
	--o-denoising-stats ASVs/stats-dada2.qza
  ```
  ## Assign taxonomy 
  ```
 qiime feature-classifier classify-sklearn --i-classifier pr2_v4_v14_classifier.qza --i-reads ASVs2/rep-seqs.qza --o-classification ASVs2/tax_sklearn.qza
  ```
  ## Convert qza files to TSV
  ```
  qiime tools export --input-path table.qza --output-path export_dir
  ```
  ```
  cd export_dir
  biom convert -i feature-table.biom -o feature-table.tsv --to-tsv
  ```
  ```
  qiime tools export --input-path tax_sklearn.qza --output-path export_dir
  ```
 
  
