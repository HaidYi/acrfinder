**<center>Acr Aca Finder</center>**

<center>(c) Yin Lab, 2019</center>

## Contents:

<a href='#installation'>I. Installation / Dependencies</a>

<a href='#about'>II. About</a>

<a href='#acr_aca_cri_runner'>III. acr_aca_cri_runner.py</a>

<a href='#acr_aca_finder'>IV. acr_aca_finder.py</a>

<a href='#crispr_cas_runner'>V. crispr_cas_runner.py</a>

<a href='#other_scripts'>VI. Other Scripts</a>

<a href='#examples'>VII. Examples</a>

<a href='#faq'>VIII. FAQ</a>

****

<br>

<div id='installation' />

## I. Installation / Dependencies

Clone/download the repository. Some dependencies are included and can be found in the <span style='color:tomato'>dependencies/</span> directory. Program expects these versions and using other versions can result in unexpected behavior.

`CRISPRCasFinder` - Already in <span style='color:tomato'>dependencies/</span> directory. To use `CRISPRCasFinder` on your machine make sure you run its install script. The manual can be found <a href='https://crisprcas.i2bc.paris-saclay.fr/Home/Download'>here</a>. Running the install script will setup paths for all the dependencies of `CRISPRCasFinder`.

It is a common problem to forget to install `CRISPRCasFinder`, so ensure that `CRISPRCasFinder` runs properly before executing <span style='color:RebeccaPurple'>acr_aca_cri_runner.py</span> to avoid errors.

`blastn` - <span style='color:RebeccaPurple'>acr_aca_cri_runner.py</span> will call/use `blastn` to search a genome. Install `blastn` from <a href='https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download'>NCBI</a>.

`psiblast+` - Used with CDD's to find Acr/Aca proteins. Install at <a href='https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download'>NCBI</a>.

`python3` - For all scripts with .py extension. Use any version at or above 3.4.

`PyGornism` - Already in <span style='color:tomato'>dependencies/</span> directory. Used to parse organism files and generate organism files in certain formats.

****

<br>


<div id='about' />

## II. About

### Acr Aca Finder is a tool used to identify Anti-CRISPR proteins (Acr) using Anti-CRISPR associative proteins (Aca)

This README file contains information about only the python scripts found in the immediate directory. These are the scripts that make up the program that identifies Acr/Aca regions.

To find out how to use other dependencies look at online sources:

`CRISPRCasFinder` - https://crisprcas.i2bc.paris-saclay.fr/CrisprCasFinder/Index

`PyGornism` - https://github.com/rtomyj/PyGornism

*`CRISPRCasFinder` is used to identify CRISPR Cas systems. This will then be used to <a href='#classification'>classify</a> Acr/Aca regions. If no CRISPR Cas systems are found within an NCID then Acr/Aca identification will not proceed (when using <span style='color:RebeccaPurple'>acr_aca_cri_runner.py</span>).

****

<br>

<div id='acr_aca_cri_runner' />

### **III. <span style='color:RebeccaPurple'>acr_aca_cri_runner.py</span>**

#### Uses <span style='color:RebeccaPurple'>acr_aca_finder.py</span> and <span style='color:RebeccaPurple'>crispr_cas_runner.py</span> to identify CRISPR Cas systems and then classify Acr/Aca proteins

### Steps

Parses the various user options.

It executes <span style='color:RebeccaPurple'>crispr_cas_runner.py</span>. It will call `CRISPRCasFinder` to identify CRISPR Cas systems. **If no CRISPRCas Systems are found in at least one NCID then the program will terminate and not continue with Acr/Aca identification**. Spacers are then used with `blastn`. The `blastn` results file contains self targeting CRISPRs avoiding targeting the CRISPRs themselves.

<span style='color:RebeccaPurple'>acr_aca_finder.py</span> is then executed. This attempts to find Acr/Aca proteins using <a href='#acr_aca_finder_filters'>filters</a>. If proteins are found then a file will be generated and the path to the file will be obtained. There will be many intermediate Acr/Aca protein files.

Using the blastn generated file the final Acr/Aca file the proteins are classified. There are three levels of classification:

<div id='classification' />

| Classification       | Meaning     |
| -------------------- | ----------- |
| Low Confidence       | There were no CRISPR Cas systems found that self targeted the organism.                   |
| Medium Confidence    | There were self targeting CRISPR Cas systems. However, the region of self targeting was not near the Acr/Aca locus. "Near" means within 5,000 BP                                                                                                       |
| High Confidence      | There were self targeting spacers. There was also a region of self targeting nearby.  |

These results are concatenated on top of the previous Acr/Aca results.

### Usage

#### Options

##### <a href='#acr_aca_options'>Acr Aca Options</a>

##### <a href='#crispr_cas_runner_options'>CRISPR Cas Options</a>

| Option | Alternative | Purpose                          |
| -----  | ----------- | -------------------------------- |
| -z     | --genomeType      | How to treat the genome. There are three options: **V**irus, **B**acteria and **A**rchaea. Viruses will not run `CRISPRCasFinder`, Archaea will run `CRISPRCasFinder` with a special Archaea flag (-ArchaCas), Bacteria will use `CRISPRCasFinder` without the Archaea flag {default = V} [string] |

#### Scripts Used

**<span style='color:RebeccaPurple'>command_options.py</span>** - defines options users can use to customize the scripts functionality

**<a href='#acr_aca_finder'>acr_aca_finder.py</a>** - finds Acr/Aca loci using predifined filters. It then further discriminates these loci using CDD's and pathogenicity databases.

**<a href='#crispr_cas_runner'>crispr_cas_runner.py</a>** -uses `CRISPRCasFinder` on the fna file provided to find CRISPR arrays. blastn is then used on the same organism to find occurences of CRISPR arrays with desired evidence level.

#### Output Files Generated

None

#### Files Updated

*<output_dir>*/*<organism_id>*_final_acr_aca.txt - updated with Acr/Aca classification. If Acr/Aca locus is high confidence it also adds the potential Cas system(s) the locus targets.

****

<br>
<div id='acr_aca_finder' />

### **IV. <span style='color:RebeccaPurple'>acr_aca_finder.py</span>**

#### Used to identify Acr and Aca proteins

##### *CAN BE RUN INDEPENDENT OF acr_aca_cri_runner

### Steps

<div id='acr_aca_finder_filters' />

Uses the following filters to identify Acr/Aca proteins:

1) Acr/Aca proteins must exist on the same strand in order to be grouped together (group = locus).
2) Acr/Aca proteins must be a length <= **m** (aa) {default = 150}.
3) Proteins must be adjacent (next to each other).
4) Adjacent proteins must be seperated by a certain threshold, distance <= **d** (bp) {default = 250}.
5) There must be at least **r** total proteins for a group to be considred an Acr/Aca locus {default = 2}.
6) There must be at least one Aca* (protein with an HTH domain) in the locus.

*Note: Aca proteins are identified using HTH HMM's found in <span style='color:tomato'>dependencies/hth_hmm</span>

Using the above filters will result in candidate Acr/Aca regions and these candidates can be further filtered. This is done by using PSSM's and looking up pathogenicity databases. The databases used are IslandViewer and PHASTER. The files used can be found in <span style='color:tomato'>dependencies/all_gis_islandviewer_iv4.txt</span>  and <span style='color:tomato'>dependencies/z_DNA_fragment_DB.header</span> respectfully.

The PSSM's used were chosen because they imply Anti-CRISPR functionality. They can be found in <span style='color:tomato'>dependencies/cdds</span>.

By default the pathogenicity databases aren't used to filter the candidate Acr/Aca regions. However, CDD's are used by default.

If a pathogenicy database is specified for use, the default method used to check is "lax". This method will include candidate Acr/Aca regions in the final result set if at least one protein in the Acr/Aca locus is found within a pathogenic region. "strict" will include candidate Acr/Aca regions in the final result set if all proteins are found within a pathogenic region.

Both the selected results from the pathogenicity databases and the CDD's are used.

**For example**

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

If using pathogenic databases (strict): locus 3 and 4 are chosen

### Usage

#### Options <div id='acr_aca_options'/>

##### Acr/Aca Loci Identification

| Option | Alternative | Purpose                     |
| -----  | ----------- | --------------------------- |
| -h     | --help      | Shows all available options |
| -m     | --aaThresh  | Max size of a protein in order to be considered Aca/Acr (aa) {default = 150} [integer] |
| -d     | --distThresh      | Max intergenic distance between proteins (bp) {default = 250} [integer] |
| -r     | --minProteins      | Min number of proteins needed per locus {default = 2} [integer] |
| -f     | --inGFF     | <span style="color:red">Required</span> Path to gff file to use/parse |
| -a     | --inFAA     | <span style="color:red">Required</span> Path to faa file to use/parse |
| -o     | --outDir    | Path to output directory to store results in. If directory doesn't exist program will attempt to create a new one with given path |
| -t     | --hthDB    | Path to file containing pressed HTH HMM obtained by running `hmmpress`. File will be used with `hmmscan`. Defualt file is provided in <span style='color:tomato'>dependencies/hth_hmm</span> |

*Note: A gff and faa file need to be specified. If no output directory is specified the current directory will be used.

##### Limiting Regions Using CDD

| Option | Alternative     | Purpose                     |
| -----  | --------------- | --------------------------- |
| -e     | --proteinUpDown | Number of surrounding proteins to use when gathering a neighborhood {default = 5} [integer] |
| -c     | --minCDDProteins | Minimum number of proteins in neighborhood that must have a CDD hit so the Acr/Aca locus can be attributed to a CDD hit {default = 2} [integer] |

*Note: CDD is used by default. To avoid using CDD use -c 0

##### Limiting Regions Using IslandViewer and/or PHASTER

| Option | Alternative | Purpose                     |
| -----  | ----------- | --------------------------- |
| -g     | --gi        | Uses IslandViewer (GI) database. {default = false} [boolean] |
| -p     | --pai       |  Uses PHASTER (PAI) database. {default = false} [boolean] |
| -s     | --strict    | All proteins in locus must lie within a region found in DB(s) being used {default = false} [boolean] |
| -l     | --lax       | Only one protein must lie within a region found in DB(s) being used {default = true} [boolean] |

*Note: Both databases can be used

#### Scripts Used

**<span style='color:RebeccaPurple'>command_options.py</span>** - defines options users can use to customize the scripts functionality

**<span style='color:RebeccaPurple'>find_candidate_acr_aca.py</span>** - finds candidate Acr/Aca regions using the preliminary filters

**<span style='color:RebeccaPurple'>parse_acr_aca_with_cdd.py</span>** - parses candidate Acr/Aca regions using cdd's

**<span style='color:RebeccaPurple'>parse_acr_aca_with_db.py</span>** - parses candidate Acr/Aca regions using pathogenicy database(s)


#### Output Files Generated

*<output_dir>*/*<organism_id>*_final_acr_aca.txt - the final set of Acr/Aca regions that passed the initial filters as well as CDD and pathogencity filtering.

*<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca.txt - potential Acr/Aca regions.

*<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca.faa - potential Acr/Aca regions in an faa format.

*<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca_neighborhood.faa - an extension of the previous file that also inludes the neighboring proteins of the potential Acr/Aca. Used with psiblast+ and CDD's

*<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca_cdd_results.txt - results from psiblast+ on potential Acr/Aca regions using CDD's

*<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca_hth_results.txt - results of hmmscan. These results contain potential HTH domains.

*<output_dir>*/intermediates/*<organism_id>*_candidate_acr_aca_parsed_hth_results.txt - a modified version of the previous version. hmmscan parser is used to clean up and obtain the highest rated results.

****

<br>
<div id='crispr_cas_runner' />

### **V. <span style='color:RebeccaPurple'>crispr_cas_runner.py</span>**

#### Executes `CRISPRCasFinder` and parses results to classify Acr/Aca loci

### Steps

Runs `CRISPRCasFinder`.

After execution is complete the file <span style='color:tomato'>result.json</span> is parsed. The CRISPR spacers with evidence level greater than or equal to the evidence level passed as an arguemnt (see below) are selected. An fna file is created of just these spacers.

***Note: If no spacers (of any confidence level) are found then execution terminates for this script and any calling scripts**

A new fna file is created from the provided organism fna file. This fna file will have all spacers 'blanked out' with the letter N <+- 500> surounding BP. These 'blanking out' will ensure that any blast queries will not target the spacers themselves. A useless result. This newly created file is called <span style='color:tomato'>masked.fna</span>.

A blast database is created using <span style='color:tomato'>masked.fna</span>.

blastn is used with the new organism fna file as the query and uses the database created from the spacers. Only results with identity >= 95% are kept.

#### Options <div id='crispr_cas_runner_options'/>

##### CRISPR Spacer Options

| Option | Alternative | Purpose                          |
| -----  | ----------- | -------------------------------- |
| -y     | --arrayEvidence      | Minimum evidence level needed of a CRISPR spacer to use {default = 3} [integer] |
| -n     | --inFNA      | <span style="color:red">Required</span> FNA file of organism being used |

#### Scripts Used By <span style='color:RebeccaPurple'>crispr_cas_runner.py</span>

**<span style='color:RebeccaPurple'>command_options.py</span>** - defines options users can use to customize the scripts functionality

**<span style='color:RebeccaPurple'>fastafy_select_spacers.py</span>** - parses <span style='color:tomato'>result.json</span> produced by `CRISPRCasFinder`. If no spacers are found the program terminates with code 0. All spacers with desired evidence level are then put into a new file, <span style='color:tomato'>spacers_with_desired_evidence.fna</span>

**<span style='color:RebeccaPurple'>mask_fna_with_spacers.py</span>** - uses the generated spacer file to mask or "blank out" CRISPR spacers to create a new FNA file that doesn't have references of the spacers. This is used to avoid self targeting the spacers with blastn. A region of +- 500 BP sorounding each spacer is also masked. The new FNA file is named <span style='color:tomato'>masked.fna</span>

#### Output Files Generated

*<output_dir>*/intermediates/masked_db/ - directory containing the db to be used by blast+

*<output_dir>*/intermediates/spacers_with_desired_evidence.fna - file containing CRISPR spacers found in the organism that have the desired evidence level. The query for blast+

*<output_dir>*/blast_out.txt - results from blast+

****

<br>

<div id='other_scripts' />

### **VI. Other scripts**

#### **<span style='color:RebeccaPurple'>acr_aca.py</span>**

Used to parse Acr/Aca results file.
Contains methods to obtain all loci one at a time as well as getting the start and ending position of a certain locus.

****

<br>

<div id='examples' />

### **VII. Examples**

**Has CRISPR Cas system**

Increased length of allowable proteins to be 300

Increased intergenic distance to 500

At the least one protein per loci

```bash
./acr_aca_cri_runner.py -f sample_organisms/ncbi_sequences/GCF_000006645.1.gff -a sample_organisms/ncbi_sequences/GCF_000006645.1.faa -n sample_organisms/ncbi_sequences/GCF_000006645.1.fna -o output_dir -m 300 -d 500 -r 1
```

<br>

Increased number of sorounding (neighbors) proteins to use with CDD to 20

Increased the number of proteins that must have a CDD hit in the neighbors to 4

```bash
./acr_aca_cri_runner.py -f sample_organisms/ncbi_sequences/GCF_000006645.1.gff -a sample_organisms/ncbi_sequences/GCF_000006645.1.faa -n sample_organisms/ncbi_sequences/GCF_000006645.1.fna -o output_dir -e 20 -c 4
```

<br>

Uses both pathogenicty databases

```bash
./acr_aca_cri_runner.py -n sample_organisms/ncbi_sequences/GCF_900090055.1.fna -a sample_organisms/ncbi_sequences/GCF_900090055.1.faa -f sample_organisms/ncbi_sequences/GCF_900090055.1.gff -o output_dir -g t -p t
```

**Doesn't have CRISPR Cas system**

```bash
./acr_aca_cri_runner.py -n sample_organisms/ncbi_sequences/GCF_000006175.1.fna -a sample_organisms/ncbi_sequences/GCF_000006175.1.faa -f sample_organisms/ncbi_sequences/GCF_000006175.1.gff  -o output_dir
```
****

<br>

<div id='faq' />

### **VIII. FAQ**

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

**Q) `CRIPSRCasFinder` is running correctly and I get the following message, "No CRISPRCas systems found. Terminating...". Why is that?**

A) There were no CRISPR Cas systems found within the same NCID/sequence. The program will not try to find Acr/Aca proteins in an organism that has no CRISPR Cas systems in the same NCID/sequence.

You can run <span style='color:RebeccaPurple'>acr_aca_finder.py</span> to just execute the Acr/Aca identification with no restrictions due to CRISPR Cas systems.
