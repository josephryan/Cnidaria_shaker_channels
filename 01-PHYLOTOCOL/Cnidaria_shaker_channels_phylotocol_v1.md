# PLANNED ANALYSES FOR PHYLOGENETICS OF CNIDARIAN SHAKER K+ CHANNELS
 Principle Investigators: Timothy Jegla, Joseph Ryan, and Adolfo Lara
 Draft or Version Number: v.1.0  
 Date: 28 Sept 2020  
 Note: updates to this document will be tracked through github
 
## 1 INTRODUCTION: BACKGROUND INFORMATION AND SCIENTIFIC RATIONALE  

### 1.1 _Background Information_
Voltage-gated potassium channels, including shaker K+ channels, are critical for tuning the intrinsic excitability of neurons. 

### 1.2 _Rationale_ 
The model anthozoan, Nematostella vectensis and the model medusozoan, Hydra magnipapillata, both have a large and diverse set of shaker Voltage-gated potassium channels (Jegla et al. 2012, Li et al. 2015). However, there are clear differences in the diversity of shaker channels between these distantly related two model cnidarians. 

### 1.3 _Objectives_   
We will identify and categorize the shaker channels from a wide range of cnidarians with representatives from each of the major cnidarian lineages (ie. Scleractinia, Actinaria, Ceriantharia, Octocorallia, Cubozoa, Scyphozoa, Staurozoa, Hydrozoa). Species with published genomes and/or deeply sequenced transcriptomes were prioritized.

## 2 STUDY DESIGN AND ENDPOINTS  

#### 2.1 The following comma separated list includes: (1) taxa included in the study, major lineage, data types (gene models are predictions from the genome)

```
Nematostella vectensis, Actinaria, gene models
Exaptasia pallida, Actinaria, gene models
Acropora digitifera, Scleractinia, gene models
Stylophora pistillata, Scleractinia, gene models
Ceriantheopsis americanus (Gulf of Mexico), Ceriantharia, transcriptome
Renilla muelleri, Octocorallia, gene models
Heliopora coerulea, Octocorallia, transcriptome
Eleutherobia rubra, Octocorallia, transcriptome
Hydra vulgaris, Hydrozoa, gene models
Clytia hemisphaerica, Hydrozoa, gene models
Morbakka virulenta, Cubozoa, gene models
Chiropsalmus quadrumanus, Cubozoa, transcriptome
Rhopilema esculentum, Scyphozoa, gene models
Sanderia malayensis, Scyphozoa, gene models
Atolla vanhoeffeni, Scyphozoa, transcriptome
Calvadosia cruxmeletensis, Staurozoa, transcriptome
Haliclystus sanjuanensis, Staurozoa, transcriptome
```

#### 2.2 Identify shaker channels

2.2.1 Reciprocal Best BLAST 
We will use the script recip.pl to conduct a reciprocal best BLAST search 
We will use BLASTP (for translated CDS of gene models) or TBLASTN (for transcriptomes) in a reciprocal framework to identify potential shaker channels in the data from the species above. As queries to BLAST, we will use all Nematostella shaker channels from Li et al. (2015). Any hits to these queries with E-values < 0.01 will be BLASTed back against the Nematostella protein set and if they hit one of the shaker channel sequences they will be included in our phylogenetic analyzes (below).

#### 2.2 Phylogenetically classify cnidarian KV channels

2.3.1 Alignment
We will align the region encompassing the BTB/POZ domain (BTB_2; PFAM:PF02214) and ion transport protein domain (Ion_trans; PFAM:PF00520) of each shaker channel sequence (including the intervening region) with default parameters in MAFFT. We will create a second alignment that includes only those sequences with site occupancy (i.e. % of non-gapped characters) of at least 75% in the regions of each of the two target domains. We will use this pruned alignment as our main data source for downstream phylogenetic classification, but we will also produce an IQTREE analysis from the “complete” dataset and report this tree along with the corresponding alignment as supplementary material in the resulting publication.

#### 2.2 Run preliminary ML tree

```
iqtree-omp -s [infile.mafft-gb] -nt AUTO -bb 1000 -m TEST -pre [output prefix] > iq.out 2> iq.err
```

#### 2.3 Prune non-KV genes

Use make_subalignment to remove non-KV genes. This script takes a prefix a subset of (our ingroup) taxa and will return an alignment of only those genes within the clade descended from the most recent common ancestor of all genes with the specified prefix. We will use Nematostella as the ingroup prefix. 
```
./make_subalignment --tree=<newick_treefile> --aln=<phylip_alignment> --root=<root_taxa> --pre=<prefix>
```

### 2.4 IQTREE
```
iqtree -m TEST -s [alignment_file] -pre [prefix] -bb 1000
```

#### 2.4  RAXML with 25 starting parsimony trees and 25 random starting trees; we will use model determined by IQTREE (in 2.4 above)

```
raxmlHPC-SSE3.PTHREADS -T 25 -p [random_number] -# 25 -m PROTGAMMA[best-fit_model] -s [alignment_file] -n [name]_mp
```
```
raxmlHPC-SSE3.PTHREADS -T 25 -d -p [random_number] -# 25 -m PROTGAMMA[best-fit_model] -s [alignment_file] -n [name]_rt
```

#### COMPARE Iqtree and 50 rax trees using rax to report the likelihood values; generate a likelihood score using RAxML for Iq-tree and grep for the likelihood values from RAxML_info files for RAxML runs.

```
raxmlHPC-SSE3 -f e -m PROTGAMMA[best-fit_model] -t [single-gene_tree] -s [alignment_file] -n [output_name]
grep 'Starting final GAMMA-based' *info*
```

### Run bootstraps and apply them to the best tree above
```
raxmlHPC -m PROTGAMMA[best-fit_model] -s [alignment_file] -p 12345 -x 12345 -# 100 -n [name]
```

```
raxmlHPC -m PROTGAMMA[best-fit_model] -p 12345 -f b -t RAxML_bestTree.[name] -z RAxML_bootstrap.[name] -n T15
```

#### 2.6 Mr. Bayes

We will run 5 MrBayes runs with the following command:

```mpirun -np 25 mb cnid_kv.nex```

and the following execution block:

```prset aamodelpr = fixed(BEST-FIT-MODEL); lset rates = gamma;
mcmcp ngen=10000000 samplefreq=10000 mcmcdiagn=yes stoprule=yes stopval=0.01
      nruns=2 nchains=5 savebrlens=yes;
mcmc;
sumt filename=FILE.nex nRuns=2 Relburnin=YES BurninFrac=.25 Contype=Allcompat;
```

### CHOOSE TREE THAT WILL BE MAIN FIGURE

Since we cannot use Bayesian principles to evaluate our ML trees, we will use ML principles to evaluate the Bayes trees. If Bayes tree has better likelihood score than the ML tree we will report the Bayes tree as our main figure with BS values from above. If our ML tree has a higher likelhood than our Bayes tree we will report the ML tree with Bayesian log likelihood scores on the branches.
All trees will be presented supplement.  And differences between Bayes and ML will be discussed.

#### 4 Work completed to-date
Sept 28 2020 -- Nothing has been completed

## 5 REFERENCES

Jegla T, Marlow HQ, Chen B, Simmons DK, Jacobo SM, Martindale MQ. Expanded functional diversity of shaker K+ channels in cnidarians is driven by gene expansion. PloS one. 2012 Dec 10;7(12):e51366.

Li X, Liu H, Luo JC, Rhodes SA, Trigg LM, Van Rossum DB, Anishkin A, Diatta FH, Sassic JK, Simmons DK, Kamel B. Major diversification of voltage-gated K+ channels occurred in ancestral parahoxozoans. Proceedings of the National Academy of Sciences. 2015 Mar 3;112(9):E1010-9.

## APPENDIX

Version : Date : Significant Revisions 


