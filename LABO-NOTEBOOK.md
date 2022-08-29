# THESIS-LICHENS: LABO NOTEBOOK

# Reads Trimming
## Fastp - FastQ
```javascript=
fastp --in1 6126-S14_R1.fastq --in2 6126-S14_R2.fastq --out1 6126-S14_R1-fastp.fastq --out2 6126-S14_R2-fastp.fastq --thread 40
#!nbreux default parameters de fastP: automatic adapters trimming, slidewindow, ...
fastqc 6126-S14_R1-fastp.fastq -t 40
fastqc 6126-S14_R2-fastp.fastq -t 40
```
# Taxonomic assignment of the reads

## Kraken2

### Script Taxo reads
```javascript=
#dans LICHENS-2021/
mkdir KRAKEN_taxo_reads/
cd ..

kraken2 --use-names --paired --db /media/vol1/databases/kraken2-nt/ reads/6126-S59_R1_001.fastq reads/6126-S59_R2_001.fastq --threads 40 --report LICHENS-2021/KRAKEN_taxo_reads/6126-S59.report > LICHENS-2021/KRAKEN_taxo_reads/6126-S59.kraken
kraken2 --use-names --paired --db /media/vol1/databases/kraken2-nt/ reads/6126-S61_R1_001.fastq reads/6126-S61_R2_001.fastq --threads 40 --report LICHENS-2021/KRAKEN_taxo_reads/6126-S61.report > LICHENS-2021/KRAKEN_taxo_reads/6126-S61.kraken
kraken2 --use-names --paired --db /media/vol1/databases/kraken2-nt/ reads/6126-S68_R1_001.fastq reads/6126-S68_R2_001.fastq --threads 40 --report LICHENS-2021/KRAKEN_taxo_reads/6126-S68.report > LICHENS-2021/KRAKEN_taxo_reads/6126-S68.kraken
kraken2 --use-names --paired --db /media/vol1/databases/kraken2-nt/ reads/6126-S64_R1_001.fastq reads/6126-S64_R2_001.fastq --threads 40 --report LICHENS-2021/KRAKEN_taxo_reads/6126-S64.report > LICHENS-2021/KRAKEN_taxo_reads/6126-S64.kraken
kraken2 --use-names --paired --db /media/vol1/databases/kraken2-nt/ reads/6126-S14_R1_001.fastq reads/6126-S14_R2_001.fastq --threads 40 --report LICHENS-2021/KRAKEN_taxo_reads/6126-S14.report > LICHENS-2021/KRAKEN_taxo_reads/6126-S14.kraken
kraken2 --use-names --paired --db /media/vol1/databases/kraken2-nt/ reads/6126-S23_R1_001.fastq reads/6126-S23_R2_001.fastq --threads 40 --report LICHENS-2021/KRAKEN_taxo_reads/6126-S23.report > LICHENS-2021/KRAKEN_taxo_reads/6126-S23.kraken
```


## Krona
```
cut -f2,3 ${f}.kraken > ${f}.krona.input
ktImportTaxonomy ${f}.krona.input -o ${f}.krona.out.html
```
Possibilité de visualiser les kraken.report sur Pavien web également.

# Assembly
## MetaSPADES
Metagenome assembly :
```javascript=
pyenv local 2.7.6
spades.py --pe1-1 `ls 6126-S14_R1-fastp.fastq` \
          --pe1-2 `ls 6126-S14_R2-fastp.fastq` \
          --meta -t 40 -m 160 -o 6126-S14-meta
```
Coverage bam file :
```javascript=
bamm make -d 6126-S14-meta/scaffolds.fasta -i 6126-S14_R1-fastp.fastq 6126-S14_R2-fastp.fastq -t 40
samtools sort scaffolds.6126-S14_R1-fastpstq.bam > scaffolds.6126-S14_R1-fastpstq.sorted.bam
samtools sort scaffolds.6126-S14_R2-fastpstq.bam > scaffolds.6126-S14_R2-fastpstq.sorted.bam
mkdir mapping
mv *.bam mapping/
cd mapping/
samtools index scaffolds.6126-S14_R1-fastpstq.sorted.bam
samtools index scaffolds.6126-S14_R2-fastpstq.sorted.bam
cd ../
```
#### MetaQUAST
Calcule stats d'assembly et recherche de novo reférencces dans la database SILVA (RNA S16, pas de cyanobactéries ni de fungi détectés)
```
metaquast.py 6126-S15-scaffolds.fasta --threads 40 --memory-efficient
```
#### QUAST ASSEMBLY
```
dans LICHENS_2021/: mkdir assembly-quast_results/
for f in `cat data.ids`;do cd assembly-quast_results/;mkdir ${f}_results/;cd ../;done
for f in `cat data.ids`;do quast.py ${f}/${f}-meta/${f}-scaffolds.fasta --threads 40 --memory-efficient -o assembly-quast_results/${f}_results;done
```
Outputs divers: summary.txt (duplication ratio, genome fraction, largest alignment, largest contig, num contigs, nums indels, num Ns per 100kbp, total aligned, total length,..), icarus.html, report.html, krona charts

# Binning

## Concoct
```javascript=
pyenv local 3.6.3
/media/vol1/apps/CONCOCT-1.1/scripts/cut_up_fasta.py 6126-S14-meta/scaffolds.fasta -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
/media/vol1/apps/CONCOCT-1.1/scripts/concoct_coverage_table.py contigs_10K.bed mapping/*.sorted.bam > coverage_table.tsv
concoct -t 40 --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/
/media/vol1/apps/CONCOCT-1.1/scripts/merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
mkdir concoct_output/fasta_bins
/media/vol1/apps/CONCOCT-1.1/scripts/extract_fasta_bins.py 6126-S14-meta/scaffolds.fasta concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins
```
## Metabat2
```javascript=
runMetaBat.sh -t 40 6126-S14-meta/scaffolds.fasta mapping/scaffolds.6126-S14_R1-fastpstq.bam mapping/scaffolds.6126-S14_R2-fastpstq.bam
```
Bins directory creation, prefix with binning software:
```
mkdir CONCOCT-bins
mkdir METABAT-bins
cp scaffolds.fasta.metabat-bins*/*.fa METABAT-bins/
cp concoct_output/fasta_bins/*.fa CONCOCT-bins/
cd METABAT-bins/
for f in *.fa; do mv $f metabat-$f; done
cd ../
cd CONCOCT-bins/
for f in *.fa; do mv $f concoct-bin.$f; done
find *.fa > ../concoct.list
cd ../
```
## Taxonomy
### GTDB
```javascript=
gtdbtk classify_wf --genome_dir METABAT-bins/ --out_dir metabat-GTDB -x fa --cpus 40
for f in `cat concoct.list`; do mkdir $f-temp/; cp CONCOCT-bins/$f $f-temp/; echo $f; gtdbtk classify_wf --genome_dir $f-temp --out_dir concoct-GTDB-$f -x fa --cpus 40; done
cat concoct-GTDB-*/gtdbtk.gtdbtk.ar122.summary.tsv > concoct-GTDB-ar122.summary.tsv
cat concoct-GTDB-*/gtdbtk.bac120.summary.tsv > concoct-GTDB-bac120.summary.tsv
cat concoct-GTDB-ar122.summary.tsv concoct-GTDB-bac120.summary.tsv > CONCOCT-GTDB.tsv
cat metabat-GTDB/gtdbtk.ar122.summary.tsv metabat-GTDB/gtdbtk.bac120.summary.tsv > METABAT-GTDB.tsv
```
## Quality : completeness & contamination
### CheckM
```javascript=
#Checkm
checkm lineage_wf -t 40 -x fa CONCOCT-bins runa > checkm-CONCOCT.out
checkm lineage_wf -t 40 -x fa METABAT-bins runb > checkm-METABAT.out
grep concoct-bin checkm-CONCOCT.out > temp
tr -s " " < temp > CONCOCT.checkm
grep metabat-bin checkm-METABAT.out > temp
tr -s " " < temp > METABAT.checkm

#fastqc
for f in CONCOCT-bins/*; do quast.py $f -t 40; done
for f in METABAT-bins/*; do quast.py $f -t 40; done
cat quast_results/*/report.txt > all-quast.results
tr -s " "  < all-quast.results > temp
mv temp all-quast.result
```
PS: bin_qa_plot retiré
```
checkm coding_plot -x fa checkm_cyanos_sel/ checkm_cyanos_sel/ bins_cyano_selected/ 
#the following arguments are required: dist_value
```
### EukCC
```javascript=
##to install via conda https://eukcc.readthedocs.io/en/latest/install.html
##wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
##bash Miniconda3-latest-Linux-x86_64.sh
##conda config --set auto_activate_base false
##conda install -c bioconda  -c conda-forge eukcc
##wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc_db_v1.1.tar.gz
##tar -xzvf eukcc_db_v1.1.tar.gz
##Use the license for gene mark .gm_key in your home dir, see install process of eukcc 
export PATH=/media/vol2/home/lcornet/miniconda3/bin:$PATH
perlbrew use perl-5.22.0-thread-multi
mkdir all-bins
cp CONCOCT-bins/* all-bins/
cp METABAT-bins/* all-bins/
cd all-bins
find *.fa > ../list
cd ../
for f in `cat list`; do eukcc --db /media/vol2/scratch/lcornet/LICHENS-2021/eukccdb --ncores 40 --ncorespplacer 10 --outdir eukcc_$f /media/vol2/scratch/lcornet/LICHENS-2021/6126-S14/all-bins/$f; done
cat eukcc_*/eukcc.tsv > all-bins_eukcc.tsv
```
### Generate Tables
```
#Generate Table
/media/vol2/scratch/lcornet/LICHENS-2021/generate-table.py list
mv Result-table.txt  6126-S14-Result-table.tsv
```
Cleaning :
```
#Cleaning
rm -rf runa runb *-temp/
rm -rf *.fastq
#rm -rf all-bins
#rm -rf list 
rm -rf concoct-GTDB_ar122.summary.tsv
rm -rf concoct-GTDB_bac120.summary.tsv
rm -rf concoct-GTDB-*/
rm -rf metabat-GTDB/
rm -rf eukcc_*
rm -rf scaffolds.fasta.metabat-bins*/
rm -rf quast_results/
```

# Bin refinement
## DAS_Tool

### Installation
```javascript=
$ conda create --name dast
$ conda activate dast
# $ conda config --add channels defaults
# $ conda config --add channels bioconda
# $ conda config --add channels conda-forge
$ conda install -c bioconda das_tool
   # environment location: /media/vol2/home/mdebras/anaconda3/envs/dast
   # DAS Tool version 1.1.3
$ DAS_Tool -h
$ Fasta_to_Scaffolds2Bin.sh -h
```

```javascript=
mkdir DAST_summary/
for f in `cat data.ids`;do cp ${f}/DASTool_Run1_DASTool_summary.txt DAST_summary/${f}_DAST_summary.txt;done
```
```
for f in `cat data.ids`;do cd ${f}/;mkdir DASTool_Run;mv DASTool_*.* DASTool_Run/.;cd ../;done
```

### Commandes
```javascript=
# Preparation of input files:

$ conda activate dast
$ cd LICHENS_2021/
$ for f in `cat data.ids`;do cd ${f}/METABAT-bins/; Fasta_to_Scaffolds2Bin.sh -e fa > metabat_scaffolds2bin.tsv;cd ../../;done
$ for f in `cat data.ids`;do cd ${f}/;mv METABAT-bins/metabat_scaffolds2bin.tsv ./metabat_scaffolds2bin.tsv;cd ../;done
$ for f in `cat data.ids`;do perl -pe "s/,/\tconcoct./g;" ${f}/concoct_output/clustering_merged.csv > concoct_scaffolds2bin.tsv;done
$ for f in `cat data.ids`;do cd ${f}/;mv concoct_output/concoct_scaffolds2bin.tsv ./concoct_scaffolds2bin.tsv;cd ../;done
```
```javascript=+
$ for f in `cat data.ids`;do cd ${f}/;
DAS_Tool -i concoct_scaffolds2bin.tsv,metabat_scaffolds2bin.tsv -l concoct,metabat -c ${f}-meta/${f}-scaffolds.fasta --search_engine diamond -o DASTool_Run1 --threads 20 --score_threshold 0.5;cd ../;done
```
```
# Labtop: 
# à partir de Collections/ qui contient le script py et les trois fichiers (tab sep scaff2bins de Concoct, Metabat, DasT) en txt : exemple:
$ python anvi-script-gen-alluvial.py --algorithm  DASTOOL --bin metabat-bin.7 > test.txt
```
A visualiser sur https://app.rawgraphs.io/ (aluvigraph)

## Bins selection
J'ai renommé les bins avec préfixe du sample
```
for f in `cat data.ids`;do cd ${f}/all-bins/;rename '' ${f}- *;cd ../../;done
```
Sub-directory pour les bins cyanobactéries:
```
mkdir cyano-bins/
for f in `cat cyano.ids`;do cp ${f} cyano-bins/.;done

cd cyano-bins/;find . -iname '*.fa' | wc -l
# ==> 54 fichiers (S36 retiré car pas de cyano)
```
Sub-direcctory pour les bins fungi:
```
mkdir fungi-bins/
for f in `cat fungi.ids`; do cp ${f} fungi-bins/.;done
# ==> 49 files (retirés: S25,30,43,45,51 + S91 très limite avec completness 70%)
```

:::warning
*Remarque* : fungi metabat souvent "splités" en ≠ bins <70% completness
\+ voir notes pour remarque de chaque sample
:::
:::info   
Q: Scaffolder les bins sélectionnées?? (SSPACE) non, augmente les erreurs d'assemblage
   :::
```javascript=
# dans reads/trimmed-reads/ :
echo PE bowtie 6126-S14_R1_001-fastp.fastq 6126-S14_R2_001-fastp.fastq 550 0.9 FR > 6126-S14_libraries.txt
# à partir directory maud-node:
mv reads/trimmed-reads/6126-S14_libraries.txt LICHENS_2021/6126-S14/all-bins/.

SSPACE.pl -l 6126-S14_libraries.txt -s 6126-S14-concoct-bin.16.fa -x 1 -m 32 -k 3 -n 15 -p 1 -b sspace
#ERROR: Invalid file in library PE: 6126-S14_R1_001-fastp.fastq -- fatal
```

## BUSCO

### Installation

```javascript=
$ conda create --name busco
$ conda activate busco
# $ conda config --add channels defaults
# $ conda config --add channels bioconda
# $ conda config --add channels conda-forge
$ conda install -c conda-forge -c bioconda busco=5.2.2
   # environment location: /media/vol2/home/mdebras/anaconda3/envs/busco
# $ busco -h
```
###  1) Cyanobacteries Nostoc
```javascript=
#dans LICHENS-2021:
mkdir BUSCO_cyanos/
cd BUSCO_cyanos/

conda activate busco
busco -m genome -i ../cyano-bins/${f}-*-bin.*.fa -o ${f} --auto-lineage --cpu 20

mkdir BUSCO_summaries/
cd ..
for f in `cat data.ids`;do cp BUSCO_cyanos/${f}/short_summary.*.nostocales_odb10.${f}.txt BUSCO_cyanos/BUSCO_summaries/.;done
python ../busco/scripts/generate_plot.py –wd ./BUSCO_cyanos/BUSCO_summaries/
```
###  2) Fungi
```javascript=
#dans LICHENS-2021/
mkdir BUSCO_fungi/
cd BUSCO_fungi/

conda activate busco
busco -m genome --augustus -i ../fungi-bins/${f}-*-bin.*.fa -o ${f} --auto-lineage --cpu 20

mkdir BUSCO_summaries/
cd ..
for f in `cat data.ids`;do cp ${f}/short_summary.*.    .${f}.txt BUSCO_fungi/BUSCO_summaries/.;done
python ../busco/scripts/generate_plot.py –wd ./BUSCO_fungi/BUSCO_summaries/

```

# BGCs Analysis
## AntiSMASH 6.0.1
### Installation
```
conda create -n antismash
```

### Commands
1) Bins cyanobacteria
```javascript=
dans LICHENS_2021/:
mkdir ANTISMASH_cyano/
for f in `cat data.ids`;do cd ANTISMASH_cyano/;mkdir ${f}/;cd ../;done

conda activate antismash
antismash --genefinding-tool prodigal --cb-general --cb-subclusters --cb-knownclusters --asf --pfam2go --smcog-trees --output-dir GENOME-MINING_cyano/antiSMASH_${f}/ cyano-bins/${f}-*.fa


```
```
scp mdebras@139.165.112.247:LICHENS_2021/ANTISMASH_cyano/6126-S14/index.html .
```

## Palantir
```
#dans LICHENS_2021/: 
mkdir PALANTIR_cyano/
```
1. Export SQL tables structuring the BGC data from AntiSMASH reports and annotated with Palantir:
:::info   
   *Q: introduire stats QUAST au tables?
   Q: répéter avec ≠ threshold ou pas? (ici: e-val par defaut)
   Q: script pour chaque samples dans une même db, pour chaque espèces, autre...?*
   Q: Trace du sample dans la db?? il faut surement renommer regions.js par ${f}_regions.js: 
   :::
    for f in `cat data.ids`;do mv ANTISMASH_cyano/${f}/regions.js ANTISMASH_cyano/${f}/${f}_regions.js;done
```javascript=
# 1 script pour chaque samples
cd PALANTIR_cyano/
nano palantir_sql_2

perlbrew use Bio-MUST-TRIAL
for f in `cat data2.ids`;do export_bgc_sql_tables.pl --infiles ../ANTISMASH_cyano/${f}/${f}_regions.js --db-name cyano_bgc_db --evalue-threshold 1e-4 --cpu 40;done

# faire -pe snode 40 sinon 1 cpu
```
2. Extract sequences from Palantir annotations in FASTA format: 
```javascript=
perlbrew use Bio-MUST-TRIAL
for f in `cat smash3.ids`;do extract_bgc_sequences.pl --report ANTISMASH_cyano/${f}/${f}_regions.js --prefix ${f} --outfile PALANTIR_cyano/${f}_extract_bcg_sequences.fasta;done
# pas d'option cpus donc qsub -pe snode 40
```
*Refaire lorsque certains types de clusters ciblés*

3. Draws NRPS/PKS gene clusters in PNG format
```
mkdir png/
```
```javascript=
perlbrew use Bio-MUST-TRIAL
for f in `cat data.ids`;do draw_bgc_maps.pl --report-file ../ANTISMASH_cyano/${f}/${f}_regions.js --outdir png/ --prefix ${f}_ --verbose;done
# pas d'option cpu
```
#658 fichiers obtenus ! (*NRPS,NRPS-like,T1PKS.png, NRPS,betalactone.png, T1PKS,hglE-KS.png, NRPS.png, ...* )

:::info   
 *Questions Palantir globales:
 Q: possible d'obtenir un format .gbk à partir des résultats Palantir?* 
 :::

Explorer les tables de la database :
1. L'importer sur mon labtop: 
scp mdebras@139.165.112.247:LICHENS_2021/PALANTIR_cyano/cyano_bgc_db .
scp -rf mdebras@139.165.112.247:LICHENS_2021/PALANTIR_cyano/cyano_bgc_db_tables/ .

## Big-SCAPE
:::info   
Notes: d'abord mettre les résultats regionXXX.gbk dans un dossier à créer (mkdir regions_gbk/) pour chaque sample.
Possibilité d'utiliser les résultats palantir plutot que antismash??
Evaluer mode mix et glocal
pour les cyano, quel output et chemin input pour bien ranger ? 
Attention: Lancer la commande python à partie du dir BIGSCAPE/ où se trouve le script bigscape.py
Info: script functions.py sert à quoi et comment ???
 :::
 ```
for f in `cat data.ids`;do cd ANTISMASH_cyano/${f}/;mkdir ${f}_regions_gbk/;cd ../../;done
for f in `cat data.ids`;do cd ANTISMASH_cyano/${f}/;mv c00*.gbk ${f}_regions_gbk/.;cd ../../;done
for f in `cat data.ids`;do cd ANTISMASH_cyano/${f}/;mkdir Big-SCAPE/;cd ../../;done
```
```javascript=
conda activate bigscape
for f in `cat data.ids`;do python bigscape.py --inputdir ../LICHENS_2021/ANTISMASH_cyano/${f}/${f}_regions_gbk/ --outputdir ../LICHENS_2021/ANTISMASH_cyano/${f}/Big-SCAPE/ --mode glocal --mix --cutoffs 0.3 0.5 0.75;done
```

outputs: cache/	html_content/  index.html*  logs/  network_files/  SVG/


# MG-RAST

Pour MG-RAST, fournir un fichier metadata avec la soumission de multiples sequences. Un metadata peut être construit avec Metazen. Metazen recquiert des modules perl:
enter the Perl shell to install the modules you want: perl -MCPAN -e shell

install JSON
install Spreadsheet::WriteExcel
install LWP::UserAgent
install CGI
Clone this GitHub repository to a location within your cgi-bin directory. The GitHub repository is located at: https://github.com/MG-RAST/metazen
Edit the config file conf/metazen_config_template.pm and save it as conf/metazen_config.pm Edit the perl header file tool_hdr_template and save it as tool_hdr
From the metazen directory run 'make'


