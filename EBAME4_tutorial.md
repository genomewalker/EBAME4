In this tutorial we will pick a random metagenome from [**MGnify**](https://www.ebi.ac.uk/metagenomics) and we will explore the functional fraction that has not been annotated.

-   [Refining the unknown](#refining-the-unknown)
-   [Exploring the non-annotated of the unannotated fraction](#exploring-the-non-annotated-of-the-unannotated-fraction)
-   [Exploring for remote homologies](#exploring-for-remote-homologies)
    -   [First unannotated sequence](#first-unannotated-sequence)
    -   [Second unannotated sequence](#second-unannotated-sequence)
    -   [Third unannotated sequence](#third-unannotated-sequence)
-   [Adding contextual data to the unknowns](#adding-contextual-data-to-the-unknowns)

Let's go to [**MGnify**](https://www.ebi.ac.uk/metagenomics) and have a look at the study **[MGYS00002304](https://www.ebi.ac.uk/metagenomics/studies/MGYS00002304)** and the sample **[ERS614327](https://www.ebi.ac.uk/metagenomics/samples/ERS614327)**

> This projects explores the functional diversity and activity of rocky subseafloor microbial communities in hydrothermal vent systems. Samples were collected in 2013 from a number of low temperature diffuse fluid vents at Axial seamount, located in the northeast Pacific Ocean. Shotgun metagenomics and metatranscriptomics were performed on four diffuse vent samples. Previous work at this site determined the taxonomic structure and distribution of microbial communities in venting fluids, but the contribution and mechanisms of the different redox driven metabolisms and the impact these reactions have on vent chemical signatures have not been fully characterized. This study helps to determine the genetic potential and expression patterns of the largely uncharacterized subseafloor microbial community and shows how these patterns change across the complicated biogeochemical gradients of hydrothermal vent systems.
>
> **NOTE**: Here is the original publication \[[**link**](https://onlinelibrary.wiley.com/doi/abs/10.1111/1462-2920.14011)]

This sample was assembled and analysed with the **[MGnify pipeline](https://www.ebi.ac.uk/metagenomics/pipelines/4.1)** (4.1). Let's have  look at the **[functional fraction](https://www.ebi.ac.uk/metagenomics/analyses/MGYA00153302#functional)** of the MGnify pipeline:

In total the pipeline identified 207,375 CDS and was able to assign a function to 122,456. The MGnify pipeline did a good job assigning a function to a 59% of the predicted CDS. What ca we do with the remaining 40%, this is a quite chunk of sequences that in most of the cases will be discarded. Let's see what we can do with those uncharacterised sequences.

> Use the Appliance [**Marine Microbiome**](https://biosphere.france-bioinformatique.fr/catalogue/appliance/125/), and load the cloud flavour **ifb.m4.2xlarge** in ifb-core-cloud

We will create a new folder for the tutorial where we will run all the analyses:

```bash
$ mkdir EBAME4_unknowns
$ cd  EBAME4_unknowns
```

Time to get the [**ORFs**](https://www.ebi.ac.uk/metagenomics/analyses/MGYA00153302#download) that couldn't be annotated by MGnify

```bash
$ wget https://www.ebi.ac.uk/metagenomics/api/v1/analyses/MGYA00153302/file/OCRB01_FASTA_CDS_unannotated.faa.gz
```

We will retrieve the [**consensus**](https://en.wikipedia.org/wiki/Consensus_sequence) sequences of the approach we showed during the lectures to tackle the unknown:

```bash
$ wget http://ebame4.metagenomics.eu/data/cluster_consensus.fasta.gz
```

## Refining the unknown

MMseqs2 needs to convert the amino acid FASTA file to their DB format \[[**+info**](https://github.com/soedinglab/MMseqs2/wiki#mmseqs2-database-format)]

```bash
$ mmseqs createdb OCRB01_FASTA_CDS_unannotated.faa.gz OCRB01_FASTA_CDS_unannotated_db
$ mmseqs createdb cluster_consensus.fasta.gz cluster_cons_db
```

Once the database is created, we will perform an iterative search, similar to a PSI-BLAST search using MMseqs2 with two iterations and an e-value cutoff of 1e-05:

```bash
$ mmseqs search OCRB01_FASTA_CDS_unannotated_db cluster_cons_db results_db tmp -e 1e-05 --num-iterations 2
```

> This would take around 8 minutes  
> **NOTE:** In a real world example we would use a more stringent threshold to assign a hit. We would use a combination of coverage and bitscore/e-value. The bitscore per column is a good alternative to evaluate the homology observed

Let's see how many hits we get:

```bash
$ wc -l results_db
1727417
```

MMseqs2 has a way to identify the best hit using an iteration within different sensitivity levels \[[**+info**](https://github.com/soedinglab/MMseqs2/wiki#how-to-find-the-best-hit-the-fastest-way)]. This takes a little longer and for the sake of time we will take a shortcut and get the first hit (ordered by e-value and bit score):

```bash
$ mmseqs filterdb results_db firsthits_db --extract-lines 1
```

> **DISCLAIMER: NOT A GOOD PRACTICE**

Next step is to export the results in a BLAST-like tabular format with **[convertalis](https://github.com/soedinglab/MMseqs2/wiki#alignment-format)**:

```bash
$ mmseqs convertalis OCRB01_FASTA_CDS_unannotated_db cluster_cons_db firsthits_db firsthits.tsv --format-mode 2
$ head firsthits.tsv
ENA-OCRB01125456-OCRB01125456.1-metagenome-genome-assembly--contig:-NODE-125456-length-241-cov-1.10753_1        5115198;kwp;kwp_c_138134;kwp_sc_1       0.566   83      33      1       1       80      89      171     3.37E-26        106     80      182
ENA-OCRB01151517-OCRB01151517.1-metagenome-genome-assembly--contig:-NODE-151517-length-222-cov-1.04192_1        15796456;gu;gu_c_278302;gu_sc_1 0.551   49      22      0       1       49      185     233     6.41E-09        54      49      233
ENA-OCRB01170121-OCRB01170121.1-metagenome-genome-assembly--contig:-NODE-170121-length-210-cov-0.941935_1       9964903;gu;gu_c_412598;gu_sc_1  0.887   62      7       0       1       62      34      95      7.45E-30        115     70      95
ENA-OCRB01034200-OCRB01034200.1-metagenome-genome-assembly--contig:-NODE-34200-length-383-cov-1.45427_2 31044474;k;HTH_AsnC-type;NA     0.433   30      17      0       1       30      26      55      9.25E-03        36      33      429
ENA-OCRB01153721-OCRB01153721.1-metagenome-genome-assembly--contig:-NODE-153721-length-220-cov-1.39394_1        25631458;k;DALR_1__tRNA-synt_1d;NA      0.684   73      23      0       1       73      416     488     5.30E-28        110     73      638
ENA-OCRB01106820-OCRB01106820.1-metagenome-genome-assembly--contig:-NODE-106820-length-256-cov-1.0995_1 1395605;k;DAGK_cat;NA   0.587   80      32      1       6       84      127     206     2.10E-25        104     85      250
ENA-OCRB01102060-OCRB01102060.1-metagenome-genome-assembly--contig:-NODE-102060-length-260-cov-1.46341_1        25568887;kwp;kwp_c_228835;kwp_sc_23896  0.350   80      52      0       3       82      86      165     8.84E-07        50      86      196
ENA-OCRB01098458-OCRB01098458.1-metagenome-genome-assembly--contig:-NODE-98458-length-264-cov-0.698565_1        8924309;gu;gu_c_52327;gu_sc_22538       0.920   88      7       0       1       88      95      182     1.01E-48        171     88      356
ENA-OCRB01101445-OCRB01101445.1-metagenome-genome-assembly--contig:-NODE-101445-length-261-cov-1.02427_2        30408861;k;Glyphos_transf;NA    0.953   43      2       0       1       43      2       44      8.21E-21        88      43      373
ENA-OCRB01043521-OCRB01043521.1-metagenome-genome-assembly--contig:-NODE-43521-length-347-cov-1.42808_1 24627654;gu;gu_c_16633;gu_sc_1  0.655   90      31      0       1       90      1       90      4.44E-31        120     90      128
```

Let's count how many unannotated ORFs have a hit to our database:

```bash
$ wc -l firsthits.tsv
45660 firsthits.tsv
```

As you can see 53% (45,660) of the ORFs that couldn't be annotated (84,919) have a hit in our database. Let's see the distribution between the different categories:

```bash
$ cut -f2 -d ';' firsthits.tsv | sort | uniq -c
  760 eu
16609 gu
15337 k
12954 kwp
```

Let's extract the clusters that had a match:

```bash
$ awk '!seen[$2]++{gsub(";","\t",$2); print $2}' firsthits.tsv  > matching_clusters.tsv
$ wc -l matching_clusters.tsv
36992 matching_clusters.tsv
```

We will use this file to explore the environmental context of these clusters (specifically GUs and EUs) in R, but first let's see what is in the ORFs we were not able to assign in one of our clusters.

## Exploring the non-annotated of the unannotated fraction

We still have a 47% of the unannotated ORFs that didn't have a hit in our DB, let's have a look what's going on...

First we are going to get the sequences that didn't have a hit:

```bash
$ awk '$3==1 {print $1}' results_db.index >  noHitSeqList
$ wc -l  noHitSeqList
39259 noHitSeqList
```

In total we have 39,259 ORFs without a hit, let's create a database that only contains these ORFs with the **createsubdb** command and a little bit of bash magic:

```bash
$ mmseqs createsubdb noHitSeqList OCRB01_FASTA_CDS_unannotated_db test_noHit
$ mmseqs createsubdb noHitSeqList OCRB01_FASTA_CDS_unannotated_db_h test_noHit_h

$ join -1 1 -2 1 <(sort -k1,1 noHitSeqList) <(sort -k1,1 OCRB01_FASTA_CDS_unannotated_db.lookup) | sort -k1,1V > test_noHit.lookup
```

Now that we have our _subdb_ we can cascade cluster it at 30% and explore the newly created clusters

```bash
$ mmseqs cluster test_noHit test_noHit_clu tmp -c 0.8 --cov-mode 0 --min-seq-id 0.3 -s 5
```

> This is would the step where you would update our clusters with the new sequences using **[clusterupdate](https://github.com/soedinglab/MMseqs2/wiki#updating-a-database-clustering-using-mmseqs-clusterupdate)** and running our curation/annotation workflow

Let's see how many clusters we have and the size:

```bash
$ mmseqs result2stats test_noHit test_noHit test_noHit_clu test_noHit_clu.sizes --stat linecount
wc -l test_noHit_clu.sizes
38411 test_noHit_clu.sizes
```

Wow, after clustering down to 30% of similarity we still have **38,411** clusters... Let's check the cluster size distribution:

```bash
$ tr -d '\0' < test_noHit_clu.sizes | sort -n | uniq -c

37738 1
  563 2
   82 3
   13 4
    6 5
    4 6
    1 7
    1 8
    2 9
    1 10
```

As you can see most of the clusters we were not able to annotate using our database are **singletons**. As we are only analysing one sample (and we don't have any abundance information) is difficult to assess if they are a real or just artefacts (we can look at the protein disorder, amino acid distribution).

The next step is to dig into the clusters we were not able to find any homology with the profile search.

First we will export the cluster results and add the number of sequences in each cluster:

```bash
$ mmseqs createtsv test_noHit test_noHit test_noHit_clu test_noHit_clu.tsv
$ cut -f1 test_noHit_clu.tsv | sort | uniq -c > test_noHit_clu_sizes.tsv
```

And then we will create consensus sequences and representative sequences for each cluster.

```bash
$ mmseqs result2profile test_noHit test_noHit test_noHit_clu test_noHit_profiles
$ mmseqs result2flat test_noHit test_noHit test_noHit_profiles_consensus test_noHit_cons.fasta

$ mmseqs result2repseq test_noHit test_noHit_clu test_noHit_clu_rep
$ mmseqs result2flat test_noHit test_noHit test_noHit_clu_rep test_noHit_rep.fasta
```

> In this case we will use the representative sequences because the clusters are very tiny. In a real world example with multiple samples and larger clusters one would prefer to use the consensus sequences

## Exploring for remote homologies

In order to dig deeper on those ORFs that we cannot annotate we will use [**JackHMMER**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-431) and [**HHBLITS**](https://www.nature.com/articles/nmeth.1818) to look for remote homologies. Although both methods are slightly different

> Rob Finn has a very nice blog post explaining the utility of Jackhmmer on **[Cryptogenomicon](https://cryptogenomicon.org/2012/04/16/interactive-iterative-searches-using-jackhmmer/)**

First we will identify some ORF clusters that have a certain size in contigs with a length larger than 300bp:

```bash
$ awk '$1>=4{split($0,a,"-"); if (a[12] >= 300){print $0}}' test_noHit_clu_sizes.tsv | sort -k1n

*9 ENA-OCRB01017322-OCRB01017322.1-metagenome-genome-assembly--contig:-NODE-17322-length-513-cov-3.52402_1
6 ENA-OCRB01014802-OCRB01014802.1-metagenome-genome-assembly--contig:-NODE-14802-length-552-cov-0.995976_2
6 ENA-OCRB01046301-OCRB01046301.1-metagenome-genome-assembly--contig:-NODE-46301-length-338-cov-12.371_1
4 ENA-OCRB01000480-OCRB01000480.1-metagenome-genome-assembly--contig:-NODE-480-length-2263-cov-12.3383_1
4 ENA-OCRB01006694-OCRB01006694.1-metagenome-genome-assembly--contig:-NODE-6694-length-795-cov-1.73919_2
*4 ENA-OCRB01014407-OCRB01014407.1-metagenome-genome-assembly--contig:-NODE-14407-length-558-cov-74.0139_1
4 ENA-OCRB01022547-OCRB01022547.1-metagenome-genome-assembly--contig:-NODE-22547-length-456-cov-2.61845_1
*4 ENA-OCRB01030292-OCRB01030292.1-metagenome-genome-assembly--contig:-NODE-30292-length-402-cov-29.9222_1
4 ENA-OCRB01053066-OCRB01053066.1-metagenome-genome-assembly--contig:-NODE-53066-length-320-cov-6.59623_1
```

We will select few of the sequences to explore the limitations of our database and approach. The main limiting factor is the coverage that we have in terms of different environments, if you remember, this metagenome comes from **hydrothermal vent systems**.

Let's take a look to three of the unannotated ORFs. Here are the parameters to use in HHblits and Jackhmmer:

-   **HHblits** will run with 3 iterations and using the [**Uniclust**](https://uniclust.mmseqs.com/) database
-   **Jackhmmer** will run with 3 iterations and using the [**UniprotKB**](https://www.uniprot.org/help/uniprotkb) database

A reminder from the HHblits manual for understanding the results and select the potential remote hits:

> Check probability and E-value: HHsearch and HHblits can detect homologous relationships far beyond the twilight zone, i.e., below 20% sequence identity. Sequence identity is therefore not an appropriate measure of relatedness anymore. The estimated probability of the template to be (at least partly) homologous to your query sequence is the most important criterion to decide whether a template HMM is actually homologous or just a high-scoring chance hit. When it is larger than 95%, say, the homology is nearly certain. Roughly speaking, one should give a hit serious consideration (i.e., check the other points in this list) whenever (1) the hit has > 50% probability, or (2) it has > 30% probability and is among the top three hits. The E-value is an alternative measure of statistical significance. It tells you how many chance hits with a score better than this would be expected if the database contained only hits unrelated to the query. At E-values below one, matches start to get marginally significant. Contrary to the probability, when calculating the E-value HHsearch and HHblits do not take into account the secondary structure similarity. Therefore, the probability is a more sensitive measure than the E-value.

## First unannotated sequence

```bash
$ grep -A1 'ENA-OCRB01017322-OCRB01017322.1-metagenome-genome-assembly--contig:-NODE-17322-length-513-cov-3.52402_1' test_noHit_rep.fasta

>ENA-OCRB01017322-OCRB01017322.1-metagenome-genome-assembly--contig:-NODE-17322-length-513-cov-3.52402_1
TMSKLMIKDLNEAQSMDHAAMSAVRGGLTVGAMTFAADQSQTIGGPGSANVGNTTAVNAATFAPSTSLTEVSPISYTEMDMATLTNVANTGVSFA
```

> **HHBLITS** results: [**here**](https://toolkit.tuebingen.mpg.de/#/jobs/2068744)  
> **Jackhmmer**  results: [**here**](https://www.ebi.ac.uk/Tools/hmmer//results/7B1998A0-C54B-11E8-8D7E-ACFCDBC3747A.3/score)  
> **Blastp** results: [**here**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=V3K1G1UH015)

This example shows a cluster that is missing in our database and in case we would update our DB it would be classified as **genomic unknown**

## Second unannotated sequence

```bash
$ grep -A1 'ENA-OCRB01014407-OCRB01014407.1-metagenome-genome-assembly--contig:-NODE-14407-length-558-cov-74.0139_1' test_noHit_rep.fasta

>ENA-OCRB01014407-OCRB01014407.1-metagenome-genome-assembly--contig:-NODE-14407-length-558-cov-74.0139_1
IATVVAATSASAFMHDNNNGWGGNNMGPFSGGNNWGPMTGGNNMGPFMGGSNAGPFSGGQNMGPFGGGSNAGPFSGAQNWGPFQGGNNWLNNTDFGTKFNTSNKTDSNASGAADGSASGSGDANAKGVADAYAKGVADAQAQAQADAYAKGYADAQASGTASSDDAAK
```

> **HHBLITS** results: [**here**](https://toolkit.tuebingen.mpg.de/#/jobs/8756404)  
> **Jackhmmer**  results: [**here**](https://www.ebi.ac.uk/Tools/hmmer/results/9604472E-C53B-11E8-81AA-E1FBDBC3747A.2/score)  
> **Blastp** results: [**here**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&VIEW_RESULTS=FromRes&RID=V3CGA9B3014&UNIQ_OBJ_NAME=A_SearchResults_1g6qrK_2ccD_dvuLMePo3sj_GTXQl_12pd7u&QUERY_INDEX=0)

Another example of the limitations of our database and the importance of having specific databases. Have a look at Mick Watson's [**behind the paper article**](https://naturemicrobiologycommunity.nature.com/users/83344-mick-watson/posts/30668-microbiome-2-0-or-what-to-do-when-you-have-no-hits-to-public-databases)  (**[Stewart et al. 2018](https://www.nature.com/articles/s41467-018-03317-6?WT.mc_id=COM_NComms_1802_Watson)**). Our workflow would classify this cluster as **genomic unknonw**

## Third unannotated sequence

```bash
$ grep -A1 'ENA-OCRB01030292-OCRB01030292.1-metagenome-genome-assembly--contig:-NODE-30292-length-402-cov-29.9222_1' test_noHit_rep.fasta

>ENA-OCRB01030292-OCRB01030292.1-metagenome-genome-assembly--contig:-NODE-30292-length-402-cov-29.9222_1
PLSVDDATALLKDVTNNANLIAHLTYGNQIVTPLTVDHILQLMNGVATHNRSWVITNLVNKNLMPSDLTP
NQVRALLGVQNNEGSVAAIKSLTDQGLMQNDLSIADAFNILGDLSDISSVGNQRQQDRDNAIRV
```

> **HHBLITS** results: [**here**](https://toolkit.tuebingen.mpg.de/#/jobs/7453964)  
> **Jackhmmer**  results: [**here**](https://www.ebi.ac.uk/Tools/hmmer/results/D81B8DE0-C54D-11E8-9F41-1E0CE976C163/score)  
> **Blastp** results: [**here**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=V3M8MWZB015)

This cluster is an example of an **environmental unknown**, there are no clear homologies in any of the databases even we use very sensible methods that look for very remote homologies.

## Adding contextual data to the unknowns

For this part we will use the contextual data and taxonomy to get some context of the clusters we have a hit.

For this part we will use the code found here: \[[**R code**](https://raw.githubusercontent.com/genomewalker/EBAME4/master/r_tutorial/ebame4_r_tutorial.R)]\[[**Notebook**](https://htmlpreview.github.io/?https://raw.githubusercontent.com/genomewalker/EBAME4/master/r_tutorial/ebame4_r_tutorial.html)]
