**<center>AcrFinder</center>**

<center>(c) <a href='http://bcb.unl.edu'>Yin Lab</a>@<a href='https://www.unl.edu'>UNL</a>2019</center>

## Contents:

<a href='#installation'>I. Installation / Dependencies</a>

<a href='#about'>II. About</a>

<a href='#using_acrfinder'>III. Using AcrFinder</a>

<a href='#docker_support'>IV. Docker Support</a>

<a href='#examples'>V. Examples</a>

<a href='#faq'>VIII. FAQ</a>

****

<div id='installation' />

## I. Installation / Dependencies

Clone/download the repository. Some dependencies are included and can be found in the <span style='color:tomato'>dependencies/</span> directory. Program expects these versions and using other versions can result in unexpected behavior.

`CRISPRCasFinder` - Already in <span style='color:tomato'>dependencies/</span> directory. To use `CRISPRCasFinder` on your machine make sure you run its install script. The manual can be found <a href='https://crisprcas.i2bc.paris-saclay.fr/Home/Download'>here</a>. Running the install script will setup paths for all the dependencies of `CRISPRCasFinder`.

It is a common problem to forget to install `CRISPRCasFinder`, so ensure that `CRISPRCasFinder` runs properly before executing <span style='color:RebeccaPurple'>acr_aca_cri_runner.py</span> to avoid errors.

`blastn` - <span style='color:RebeccaPurple'>acr_aca_cri_runner.py</span> will call/use `blastn` to search a genome. Install `blastn` from <a href='https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download'>NCBI</a>.

`psiblast+` - Used with CDD to find mobilome proteins. Install at <a href='https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download'>NCBI</a>.

`python3` - For all scripts with .py extension. Use any version at or above 3.4.

`PyGornism` - Already in <span style='color:tomato'>dependencies/</span> directory. Used to parse organism files and generate organism files in certain formats.

****


<div id='about' />

## II. About

### AcrFinder is a tool used to identify Anti-CRISPR proteins (Acr) using both sequence homology and guilt-by-association approaches.

This README file contains information about only the python scripts found in the current directory. These are the scripts that are used to identify genomic loci that contain Acr and/or Aca homologs.

To find out how to use other dependencies look at online sources:

`CRISPRCasFinder` - https://crisprcas.i2bc.paris-saclay.fr/CrisprCasFinder/Index

*`CRISPRCasFinder` is used to identify CRISPR Cas systems. This will then be used to <a href='#classification'>classify</a> the genomic loci that contain Acr and/or Aca homologs. If no CRISPR Cas systems are found within a genome, then only homology based search will be implemented for Acr homologs.

****

<div id='using_acrfinder' />

### **III. <span style='color:RebeccaPurple'>Using AcrFinder</span>**

#### Input

AcrFinder needs **.fna**, **.gff** and **.faa** as input. Only **.fna** file as input is also acceptable; in that case, the **.gff** and **.faa** file will be generated by running <a href='https://github.com/hyattpd/Prodigal'>Prodigal</a>.

#### List of Options

| Option | Alternative | Purpose                     |
| -----  | ----------- | --------------------------- |
| -h     | --help      | Shows all available options |
| -n     | --inFNA      | <span style="color:red">Required</span> fna file |
| -f     | --inGFF     | <span style="color:red">Required</span> Path to gff file to use/parse |
| -a     | --inFAA     | <span style="color:red">Required</span> Path to faa file to use/parse |
| -m     | --aaThresh  | Max size of a protein in order to be considered Aca/Acr (aa) {default = 150} [integer] |
| -d     | --distThresh      | Max intergenic distance between proteins (bp) {default = 250} [integer] |
| -r     | --minProteins      | Min number of proteins needed per locus {default = 2} [integer] |
| -y     | --arrayEvidence      | Minimum evidence level needed of a CRISPR spacer to use {default = 3} [integer] |
| -o     | --outDir    | Path to output directory to store results in. If not provided, the program will attempt to create a new one with given path |
| -t     | --aca    | Known Aca file (.faa) to diamond candidate aca in candidate Acr-Aca loci |
| -u     | --acr    | Known Acr file (.faa) to diamond the homolog of Acr |
| -z     | --genomeType      | How to treat the genome. There are three options: **V**irus, **B**acteria and **A**rchaea. Viruses will not run `CRISPRCasFinder`, Archaea will run `CRISPRCasFinder` with a special Archaea flag (-ArchaCas), Bacteria will use `CRISPRCasFinder` without the Archaea flag {default = V} [string] |
| -e     | --proteinUpDown | Number of surrounding (up- and down-stream) proteins to use when gathering a neighborhood {default = 5} [integer] |
| -c     | --minCDDProteins | Minimum number of proteins in neighborhood that must have a CDD mobilome hit so the Acr/Aca locus can be attributed to a CDD hit {default = 2} [integer] |
| -g     | --gi        | Uses IslandViewer (GI) database. {default = false} [boolean] |
| -p     | --pai       |  Uses PHASTER (prophage) database. {default = false} [boolean] |
| -s     | --strict    | All proteins in locus must lie within a region found in DB(s) being used {default = false} [boolean] |
| -l     | --lax       | Only one protein must lie within a region found in DB(s) being used {default = true} [boolean] |

#### Output 

#### Classification

There are three levels of classification in output:

| Classification       | Meaning     |
| -------------------- | ----------- |
| Low Confidence       | If this Acr-Aca locus has a CRISPR-Cas locus but no self-targeting spacers in the genome, it is labeled as “low confidence” and inferred to target the CRISPR-Cas locus.                   |
| Medium Confidence    | If this Acr-Aca locus has a self-targeting spacer target in the genome but not nearby, it is labeled as “medium confidence” and inferred to target the CRISPR-Cas locus with the self-targeting spacer. "Nearby" means within 5,000 BP.                                                                                                       |
| High Confidence      | If this Acr-Aca locus has a nearby self-targeting spacer target, it is labeled as “high confidence” and inferred to target the CRISPR-Cas locus with the self-targeting spacer.  |

#### Ouput files

| Name                 | Meaning     |
| -------------------- | ----------- |
|*<output_dir>*/CRISPRCas_OUTPUT   | The output folder of CRISPRCasFinder |
|*<output_dir>*/subjects           | The folder that contains the input files |
|*<output_dir>*/intermediates      | The folder that contains intermediate result files |
|*<output_dir>*/blast_out.txt      | Results from blast+ |
|*<output_dir>*/*<organism_id>*_final_acr_aca.txt | The final set of Acr/Aca regions that passed the initial filters as well as CDD and pathogencity filtering. |
|*<output_dir>*/*<organism_id>*_final_acr_homolog.fasta | The final set of proteins that are aligned to Acr database within given similarity threshold in orgainism. |
|*<output_dir>*/intermediates/masked_db/ | The directory containing the db to be used by blast+ |
|*<output_dir>*/intermediates/spacers_with_desired_evidence.fna | The file containing CRISPR spacers found in the organism that have the desired evidence level. The query for blast+ |
|*<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca.txt | Potential Acr/Aca regions. |
|*<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca.faa | Potential Acr/Aca regions in an faa format. |
|*<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca_neighborhood.faa | An extension of the previous file that also inludes the neighboring proteins of the potential Acr/Aca. Used with psiblast+ and CDD's |
|*<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca_cdd_results.txt | Results from psiblast+ on potential Acr/Aca regions using CDD's |
|*<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca_diamond_result.txt | Results of diamond. These results can be aligned to **Aca database** within given similarity threshold using protein in *<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca.faa. |
|*<output_dir>*/intermediates/*<organism_id>*_candidate_acr_homolog_result.txt | Results of diamond. These results can be aligned to **Acr database** within given similarity threshold using protein in *<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca.faa. |
|*<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca_diamond_database.dmnd | Database of diamond made from *<organism_id>*_candidate_acr_aca.faa file. |
|*<output_dir>*/intermediates/*<organism_id>*_acr_homolog_result.txt | Results of diamond. These results can be aligned to **Acr database** within given similarity threshold using proteins in *<output_dir>*/subjects/*<organism_id>*_protein.faa. |
|*<output_dir>*/intermediates/*<organism_id>*_acr_homolog_result.fasta | *Protein Sequence* file (*.faa*) of protein in *<output_dir>*/intermediates/*<organism_id>*_acr_homolog_result.txt |
|*<output_dir>*/intermediates/*<organism_id>*_acr_diamond_database.dmnd | Database of diamond made from *<output_dir>*/subjects/*<organism_id>*_protein.faa file |

****

<!-- 
Uses the following filters to identify Acr/Aca proteins:

1) Acr/Aca proteins must exist on the same strand in order to be grouped together (group = locus).
2) Acr/Aca proteins must be a length <= **m** (aa) {default = 150}.
3) Proteins must be adjacent (next to each other).
4) Adjacent proteins must be seperated by a certain threshold, distance <= **d** (bp) {default = 250}.
5) There must be at least **r** total proteins for a group to be considred an Acr/Aca locus {default = 2}.
6) There must be at least one Aca* (protein with an HTH domain) in the locus.

*Note: Aca proteins are identified using HTH HMM's found in <span style='color:tomato'>dependencies/hth_hmm</span> -->

<!-- Using the above filters will result in candidate Acr/Aca regions and these candidates can be further filtered. This is done by using PSSM's and looking up pathogenicity databases. The databases used are IslandViewer and PHASTER. The files used can be found in <span style='color:tomato'>dependencies/all_gis_islandviewer_iv4.txt</span>  and <span style='color:tomato'>dependencies/z_DNA_fragment_DB.header</span> respectfully.

The PSSM's used were chosen because they imply Anti-CRISPR functionality. They can be found in <span style='color:tomato'>dependencies/cdds</span>.

By default the pathogenicity databases aren't used to filter the candidate Acr/Aca regions. However, CDD's are used by default.

If a pathogenicy database is specified for use, the default method used to check is "lax". This method will include candidate Acr/Aca regions in the final result set if at least one protein in the Acr/Aca locus is found within a pathogenic region. "strict" will include candidate Acr/Aca regions in the final result set if all proteins are found within a pathogenic region.

Both the selected results from the pathogenicity databases and the CDD's are used. -->

<!-- **For example**

| Locus # | Protein | Within Pathogenic region?   | CDD hit?   |
| -----  | ----------- | ------------------------ | ---------- |
| 1     | WP00001      | NO | NO |
| 1     | WP00002      | NO | NO |
| 1     | WP00003      | NO | NO |
| 2     | WP00010      | YES | NO |
| 2     | WP00011      | NO | NO |
| 3     | WP00021      | YES | NO |
| 3     | WP00022      | YES | NO |
| 3     | WP00022      | YES | NO |
| 4     | WP00031      | NO | YES |
| 4     | WP00032      | NO | YES |
| 4     | WP00033      | NO | YES |

If using no pathogenicity databases: Locus 4 alone is chosen

If using pathogenic databases (lax): locus 2,3 and 4 are chosen

If using pathogenic databases (strict): locus 3 and 4 are chosen -->

<!-- ### Usage -->
<!-- 
#### Options <div id='acr_aca_options'/>

##### Acr/Aca Loci Identification -->

<!-- | Option | Alternative | Purpose                     |
| -----  | ----------- | --------------------------- |
| -h     | --help      | Shows all available options |
| -m     | --aaThresh  | Max size of a protein in order to be considered Aca/Acr (aa) {default = 150} [integer] |
| -d     | --distThresh      | Max intergenic distance between proteins (bp) {default = 250} [integer] |
| -r     | --minProteins      | Min number of proteins needed per locus {default = 2} [integer] |
| -f     | --inGFF     | <span style="color:red">Required</span> Path to gff file to use/parse |
| -a     | --inFAA     | <span style="color:red">Required</span> Path to faa file to use/parse |
| -o     | --outDir    | Path to output directory to store results in. If directory doesn't exist program will attempt to create a new one with given path |
| -t     | --hthDB    | Path to file containing pressed HTH HMM obtained by running `hmmpress`. File will be used with `hmmscan`. Defualt file is provided in <span style='color:tomato'>dependencies/hth_hmm</span> |

*Note: A gff and faa file need to be specified. If no output directory is specified the current directory will be used. -->

<!-- ##### Limiting Regions Using CDD

| Option | Alternative     | Purpose                     |
| -----  | --------------- | --------------------------- |
| -e     | --proteinUpDown | Number of surrounding proteins to use when gathering a neighborhood {default = 5} [integer] |
| -c     | --minCDDProteins | Minimum number of proteins in neighborhood that must have a CDD hit so the Acr/Aca locus can be attributed to a CDD hit {default = 2} [integer] | -->

<!-- *Note: CDD is used by default. To avoid using CDD use -c 0 -->

<!-- ##### Limiting Regions Using IslandViewer and/or PHASTER

| Option | Alternative | Purpose                     |
| -----  | ----------- | --------------------------- |
| -g     | --gi        | Uses IslandViewer (GI) database. {default = false} [boolean] |
| -p     | --pai       |  Uses PHASTER (PAI) database. {default = false} [boolean] |
| -s     | --strict    | All proteins in locus must lie within a region found in DB(s) being used {default = false} [boolean] |
| -l     | --lax       | Only one protein must lie within a region found in DB(s) being used {default = true} [boolean] | -->
<!-- 
*Note: Both databases can be used -->

<!-- #### Scripts Used

**<span style='color:RebeccaPurple'>command_options.py</span>** - defines options users can use to customize the scripts functionality

**<span style='color:RebeccaPurple'>find_candidate_acr_aca.py</span>** - finds candidate Acr/Aca regions using the preliminary filters

**<span style='color:RebeccaPurple'>parse_acr_aca_with_cdd.py</span>** - parses candidate Acr/Aca regions using cdd's

**<span style='color:RebeccaPurple'>parse_acr_aca_with_db.py</span>** - parses candidate Acr/Aca regions using pathogenicy database(s) -->

<div id='docker_support' />

### **IV. <span style='color:RebeccaPurple'>Docker Support</span>**

To help users to configure the environment to use the software easily, we provide the _.Dockerfile_ can be used using the command:

```bash
cd path/to/repo
docker build -t [tag name] .
```

or you can pull the image from docker hub using the command:

```bash
docker pull [OPTIONS] haidyi/AcrFinder:latest
```

<!-- Runs `CRISPRCasFinder`.

After execution is complete the file <span style='color:tomato'>result.json</span> is parsed. The CRISPR spacers with evidence level greater than or equal to the evidence level passed as an arguemnt (see below) are selected. An fna file is created of just these spacers.

***Note: If no spacers (of any confidence level) are found then execution terminates for this script and any calling scripts**

A new fna file is created from the provided organism fna file. This fna file will have all spacers 'blanked out' with the letter N <+- 500> surounding BP. These 'blanking out' will ensure that any blast queries will not target the spacers themselves. A useless result. This newly created file is called <span style='color:tomato'>masked.fna</span>.

A blast database is created using <span style='color:tomato'>masked.fna</span>.

blastn is used with the new organism fna file as the query and uses the database created from the spacers. Only results with identity >= 95% are kept. -->

<!-- #### Options <div id='crispr_cas_runner_options'/>

#### Scripts Used By <span style='color:RebeccaPurple'>crispr_cas_runner.py</span>

**<span style='color:RebeccaPurple'>command_options.py</span>** - defines options users can use to customize the scripts functionality

**<span style='color:RebeccaPurple'>fastafy_select_spacers.py</span>** - parses <span style='color:tomato'>result.json</span> produced by `CRISPRCasFinder`. If no spacers are found the program terminates with code 0. All spacers with desired evidence level are then put into a new file, <span style='color:tomato'>spacers_with_desired_evidence.fna</span>

**<span style='color:RebeccaPurple'>mask_fna_with_spacers.py</span>** - uses the generated spacer file to mask or "blank out" CRISPR spacers to create a new FNA file that doesn't have references of the spacers. This is used to avoid self targeting the spacers with blastn. A region of +- 500 BP sorounding each spacer is also masked. The new FNA file is named <span style='color:tomato'>masked.fna</span> -->

****

<!-- <br>

<div id='other_scripts' />

#### **<span style='color:RebeccaPurple'>acr_aca.py</span>**

Used to parse Acr/Aca results file.
Contains methods to obtain all loci one at a time as well as getting the start and ending position of a certain locus.

**** -->


<div id='examples' />

### **V. Examples**

```bash
python3 acr_aca_cri_runner.py -n sample_organisms/GCF_000210795.2/GCF_000210795.2_genomic.fna -f sample_organisms/GCF_000210795.2/GCF_000210795.2_genomic.gff -a sample_organisms/GCF_000210795.2/GCF_000210795.2_protein.faa -o output_dir -z B
```
or you can only use **.fna** file as input.

```bash
python3 acr_aca_cri_runner.py -n sample_organisms/GCF_000210795.2/GCF_000210795.2_genomic.fna -o output_dir -z B
```

#### Use Docker Image
<!-- 
Increased number of sorounding (neighbors) proteins to use with CDD to 20

Increased the number of proteins that must have a CDD hit in the neighbors to 4 -->
```bash
docker run [OPTIONS] [NAME:TAG] python3 acr_aca_cri_runner.py -n sample_organisms/GCF_000210795.2/GCF_000210795.2_genomic.fna -f sample_organisms/GCF_000210795.2/GCF_000210795.2_genomic.gff -a sample_organisms/GCF_000210795.2/GCF_000210795.2_protein.faa -o output_dir -z B
```

****

<!-- <br>

Uses both pathogenicty databases

```bash
./acr_aca_cri_runner.py -n sample_organisms/ncbi_sequences/GCF_900090055.1.fna -a sample_organisms/ncbi_sequences/GCF_900090055.1.faa -f sample_organisms/ncbi_sequences/GCF_900090055.1.gff -o output_dir -g t -p t
```

**Doesn't have CRISPR Cas system**

```bash
./acr_aca_cri_runner.py -n sample_organisms/ncbi_sequences/GCF_000006175.1.fna -a sample_organisms/ncbi_sequences/GCF_000006175.1.faa -f sample_organisms/ncbi_sequences/GCF_000006175.1.gff  -o output_dir
```
****

<br> -->

<div id='faq' />

### **VI. FAQ**

**Q) I ran <span style='color:RebeccaPurple'>acr_aca_cri_runner.py</span> and I got errors that pertain to CRISPR/Cas. Whats the issue?**

A) Make sure `CRIPSRCasFinder` is installed properly. `CRIPSRCasFinder` has many dependencies of its own and will only work if they are all installed correctly. A good indicator of a correctly installed `CRIPSRCasFinder` is the following terminal output:

```bash
################################################################
# --> Welcome to dependencies/CRISPRCasFinder/CRISPRCasFinder.pl (version 4.2.17)
################################################################


vmatch2 is...............OK
mkvtree2 is...............OK
vsubseqselect2 is...............OK
fuzznuc (from emboss) is...............OK
needle (from emboss) is...............OK
```

<!-- **Q) `CRIPSRCasFinder` is running correctly and I get the following message, "No CRISPRCas systems found. Terminating...". Why is that?**

A) There were no CRISPR Cas systems found within the same NCID/sequence. The program will not try to find Acr/Aca proteins in an organism that has no CRISPR Cas systems in the same NCID/sequence.

You can run <span style='color:RebeccaPurple'>acr_aca_finder.py</span> to just execute the Acr/Aca identification with no restrictions due to CRISPR Cas systems. -->

