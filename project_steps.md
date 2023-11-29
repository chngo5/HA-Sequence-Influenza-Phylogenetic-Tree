This project involves creating a phylogenetic tree from 17 strains of influenza found in humans.
To begin organization, first I created a Final Project directory in my terminal. 

```
#create Final Project directory
mkdir FinalProject
```
![Final Project directory steps.PNG](https://github.com/chngo5/HA-Sequence-Influenza-Phylogenetic-Tree/blob/c04a582ab02ea4747313aa90a2e7dcb9dbffcfee/Final%20Project%20directory%20steps.PNG)

Then, within that directory, I created a separate directory labeled "Strains."

```
#create Strains directory within FinalProject
mkdir Strains
```
![make Strains folder in Final Project directory.PNG]
Then, I imported data for the 17 strains used. All data was retrieved from the NCBI Influenza Virus nucleotide sequences search located at this link: https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database

The specific filters used for the data were: Human for Host, HA for Segment, 2022 in both Collection Date fields. Other fields were varied to collect needed sequences.

Then, I made text files containing the FASTA sequences for each strain. I labeled the strains with headers matching identifying information in the sequences found on NCBI.

The data was input into text files using 'nano' for 17 HA nucleotide sequences.

```
#nucleotide sequence names for 17 strain files
A_India_CG-AIIMSR-292_2022_H1N1
A_Tokyo_13424_2022_H3N2
A_China_CSKFQ-22-5_2022_H3N8
A_China_ZMD-22-2_2022_H3N8
A_Yangzhou_125_2022_H5N6
A_Alaska_USAFSAM-14036_2022_H1N1
A_Arizona_67_2022_H1N1
A_California_168_2022_H1N1
A_New_York_52_2022_H1N1
A_Texas_80_2022_H1N1
A_Friuli-Venezia_Giulia_USAFSAM-13527_2022_H3N2
A_Rheinland-Pfalz_USAFSAM-13452_2022_H1N1
A_Suffolk_13426_2022_H1N1
A_Rheinland-Pfalz_USAFSAM-13689_2022_H1N1
A_Victoria_4897_2022_H1N1
A_Wisconsin_588_2019_H1N1_vacc
A_Darwin_6_2021_H3N2_vacc
```

Then, I concatenated all sequences in the files above into a single file named "all_strains."

```
#cat all files into one file named all_strains
cat A_India_CG-AIIMSR-292_2022_H1N1 A_Tokyo_13424_2022_H3N2 A_China_CSKFQ-22-5_2022_H3N8 A_China_ZMD-22-2_2022_H3N8 A_Yangzhou_125_2022_H5N6 A_Alaska_USAFSAM-14036_2022_H1N1 A_Arizona_67_2022_H1N1 A_California_168_2022_H1N1 A_New_York_52_2022_H1N1 A_Texas_80_2022_H1N1 A_Friuli-Venezia_Giulia_USAFSAM-13527_2022_H3N2 A_Rheinland-Pfalz_USAFSAM-13452_2022_H1N1 A_Suffolk_13426_2022_H1N1 A_Rheinland-Pfalz_USAFSAM-13689_2022_H1N1 A_Victoria_4897_2022_H1N1 A_Wisconsin_588_2019_H1N1_vacc A_Darwin_6_2021_H3N2_vacc > all_strains
```

Then, I activated the conda environment where the mafft and iqtree software installs are located.
```
#activate conda environment with mafft/iqtree
conda activate lab9_path/lab9_conda
```

Then, I specified the path to mafft.
```
#mafft path
/project/stuckert/chngo5/Lab9/lab9_path/lab9_conda/bin/mafft
```
Then, I specified the path to iqtree.
```
#iqtree path
/project/stuckert/chngo5/Lab9/lab9_path/lab9_conda/bin/iqtree
```
Then, I ran mafft.
```
#mafft
/project/stuckert/chngo5/Lab9/lab9_path/lab9_conda/bin/mafft --auto all_strains > output.fa
```
Then, I ran iqtree.
```
#iqtree
/project/stuckert/chngo5/Lab9/lab9_path/lab9_conda/bin/iqtree -s output.fa -m HKY -bb 1000 -pre result
```
The output of the iqtree.contree file needed for the phylogenetic tree I saved as follows:
```
cat result.contree
```
(A_Alaska_USAFSAM-14036_2022_H1N1:0.0045713784,((((A_Arizona_67_2022_H1N1:0.0005712333,(A_New_York_52_2022_H1N1:0.0017160797,A_Texas_80_2022_H1N1:0.0034348686)67:0.0000020787)65:0.0000020787,((A_California_168_2022_H1N1:0.0000064057,A_Victoria_4897_2022_H1N1:0.0005881445)70:0.0000020801,A_Rheinland-Pfalz_USAFSAM-13689_2022_H1N1:0.0000020787)84:0.0005631675)90:0.0011814886,A_Suffolk_13426_2022_H1N1:0.0029522434)96:0.0044771878,((A_India_CG-AIIMSR-292_2022_H1N1:0.0000020787,A_Wisconsin_588_2019_H1N1_vacc:0.0000020787)100:0.0074721767,A_Rheinland-Pfalz_USAFSAM-13452_2022_H1N1:0.0074143619)78:0.0035659752)54:0.0016326210,(((A_China_CSKFQ-22-5_2022_H3N8:0.0052482923,A_China_ZMD-22-2_2022_H3N8:0.0085978397)100:0.1140944664,(A_Darwin_6_2021_H3N2_vacc:0.0031525641,(A_Friuli-Venezia_Giulia_USAFSAM-13527_2022_H3N2:0.0035211137,A_Tokyo_13424_2022_H3N2:0.0047360511)95:0.0046053682)100:0.1307725254)100:0.5594455464,A_Yangzhou_125_2022_H5N6:0.2789883607)100:0.2176715540);

Then, I moved to R Studio to make the phylogenetic tree. 
```
#R studio with 17 strains
library(ape)
library(ggtree)
shh_tree <- "(A_Alaska_USAFSAM-14036_2022_H1N1:0.0045713785,((((A_Arizona_67_2022_H1N1:0.0005712333,(A_New_York_52_2022_H1N1:0.0017160799,A_Texas_80_2022_H1N1:0.0034348684)77:0.0000020787)66:0.0000020787,((A_California_168_2022_H1N1:0.0000064057,A_Victoria_4897_2022_H1N1:0.0005881445)73:0.0000020801,A_Rheinland-Pfalz_USAFSAM-13689_2022_H1N1:0.0000020787)68:0.0005631675)89:0.0011814886,A_Suffolk_13426_2022_H1N1:0.0029522434)97:0.0044771878,((A_India_CG-AIIMSR-292_2022_H1N1:0.0000020787,A_Wisconsin_588_2019_H1N1_vacc:0.0000020787)100:0.0074721767,A_Rheinland-Pfalz_USAFSAM-13452_2022_H1N1:0.0074143619)75:0.0035659752)55:0.0016326209,(((A_China_CSKFQ-22-5_2022_H3N8:0.0052482920,A_China_ZMD-22-2_2022_H3N8:0.0085978400)100:0.1140944768,(A_Darwin_6_2021_H3N2_vacc:0.0031525639,(A_Friuli-Venezia_Giulia_USAFSAM-13527_2022_H3N2:0.0035211137,A_Tokyo_13424_2022_H3N2:0.0047360511)95:0.0046053685)100:0.1307725204)100:0.5594456325,A_Yangzhou_125_2022_H5N6:0.2789883719)100:0.2176715749);"
shh_tree <- read.tree(text=shh_tree)
shh_tree <- ape::root.phylo(shh_tree, outgroup = "A_Yangzhou_125_2022_H5N6")
ggtree(shh_tree) + geom_tiplab()+ xlim(0, 1.75)
```
The resulting tree looks as follows:
