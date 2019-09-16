
 # ARF: Algorithmic rRNA Filtering
========================

## Rationale
Determining the taxonomic composition of a microbial community by 16s rRNA sequencing requires a high-quality repository of 16s rRNA alleles with accurate taxonomic identificaiton.

Hand-curated repositories, such as [SILVA](https://www.arb-silva.de/), [RDP](https://rdp.cme.msu.edu), and [greengenes](http://greengenes.secondgenome.com/) are high-quality sources. ARF is a purely algorithmic approach to identifying 16S rRNA alleles with taxonomic annotation that are internally consistent and is meant as an adjunct to these hand-curated sources.


## Approach
In broad strokes, ARF downloads 16s all bacterial and archaeal rRNA sequences<sup id="a1">[1](#f1)</sup> directly from [NCBI NT database](https://www.ncbi.nlm.nih.gov/nucleotide/), verifies that sequences annotated as being 16s rRNA actually are 16s rRNA genes, full length (at least 1200bp), and with minimal ambiguous bases. This becomes the `1200bp` set of rRNA reads, which are valid rRNA alleles, some with and some without verified taxonomic annotations.

The `1200bp` set is further subsetted into a `named` set of 16s rRNA genes that have some _species_ level taxonomic assignment (versus those rRNA alleles assigned to _uncultured bacterium_ or such).

The `named` set is further subsetted into a `filtered` subset using the [deenurp](https://github.com/fhcrc/deenurp) `filter-outliers` mode. Alleles are grouped by their taxonomic assignment. Each group is clustered at the sequence level, and outliers are filtered out.(Multiple centroids are tolerated.) The basic concept is majority rules.

In parallel type strains are found in the `named` set and placed into the `types` subset.

`ARF` is implemented as a [nextflow](https://nextflow.io) workflow. Nextflow can (and should) iteratively update an extant library of sequences. 

## Usage:
```
nextflow run jgolob/arf <ARGUMENTS>

Required Arguments:
--repo                          path to directory holding the current repo (default = './arf')
--out                           path where refreshed repo should be placed (default = './refreshed')
--email                         Valid email to use with NCBI
--ncbi_concurrent_connections   Number of concurrent connections (default = 3)
--retry_max                     Max retries with NCBI requests (default = 1)
--retry_delay                   Delay (ms) between retries (default=60000)
--min_len                       Minimum annotated length of 16s rRNA to even be downloaded (default 500)
--species_cap                   Maximum number of 16s rRNA genes per annotated species (default 5000)


Options:
--api_key                       NCBI api key (will increase download rate)
--debug                         Cuts down the number of records for testing
```

## Output
```
arf/
├── dedup
│   └── 1200bp
│       ├── named
│       │   ├── blast.nhr
│       │   ├── blast.nin
│       │   ├── blast.nsq
│       │   ├── filtered
│       │   │   ├── blast.nhr
│       │   │   ├── blast.nin
│       │   │   ├── blast.nsq
│       │   │   ├── lineages.csv
│       │   │   ├── lineages.txt
│       │   │   ├── outliers.csv
│       │   │   ├── seq_info.csv
│       │   │   ├── seqs.fasta
│       │   │   └── taxonomy.csv
│       │   ├── lineages.csv
│       │   ├── lineages.txt
│       │   ├── seq_info.csv
│       │   ├── seqs.fasta
│       │   └── taxonomy.csv
│       ├── seq_info.csv
│       ├── seqs.fasta
│       └── types
│           ├── blast.nhr
│           ├── blast.nin
│           ├── blast.nsq
│           ├── lineages.csv
│           ├── lineages.txt
│           ├── seq_info.csv
│           ├── seqs.fasta
│           └── taxonomy.csv
├── pubmed_info.csv
├── records.txt
├── references.csv
├── refseq_info.csv
├── seq_info.csv
├── seqs.fasta
├── taxdmp.zip
├── taxonomy.csv
└── taxonomy.db
```

Outputs are intended to be compatible with [mothur](https://www.mothur.org) and [MaLiAmPi](https://github.com/jgolob/maliampi).

## Disclosures
ARF makes tails wag and is a collaboration between Noah Hoffman and Chris Rosenthal of Laboratory Medicine at the University of Washington and Fred Hutch and Jonathan Golob of the University of Michigan Department of Medicine and Division of Infectious Diseases. Use it at your own risk.

<b id="f1">1</b>  [↩](#a1) The exact NCBI NT search strings are: 

1) Bacteria: `16s[All Fields] AND rRNA[Feature Key] AND Bacteria[Organism] AND 500 : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism])`

2) Archaea: `16s[All Fields] AND rRNA[Feature Key] AND Archaea[Organism] AND 500 : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism])`

3) Type strains: `16s[All Fields] AND rRNA[Feature Key] AND (Bacteria[Organism] OR  Archaea[Organism]) AND 500 : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism]) AND sequence_from_type[Filter]`


