#!/usr/bin/env nextflow

/*
  ya16sdb: Download and curate the 16S rRNA sequences from NCBI
    Originally developed by Noah Hoffman and Chris Rosenthal
    Adapted to nextflow from sconstruct by Jonathan Golob

    This workflow interatively updates a cache of 16S rRNA genes
    extracted from the NCBI NT database, verifies the annotations,
    sequence quality (i.e. gene length and ambiguous bases), and
    then validates the annotations of the 16S rRNA genes by
    'majority rules'.

    This is an approach to generate algorithmically a large,
    frequently updated and conservatively annotated 16S rRNA 
    library for use in 16S rRNA amplicon classification.

    Steps:
    1) Retrieve an updated list of record ids (i.e. versions) 
        with annotated 16S rRNA gene from NCBI NT, 
        as well as the seq start and stop when relevant.
        -> Do this for all bacteria, type strains, and archaea
    2) (If available), compare this to the starting list of record ids
        (versions) *already* in the repo.
    3) Mask away any known bad/malformed/invalid record ids
    4) Retrieve the new (i.e not in the repo already) records
    5) Parse the records, extract the 16S rRNA sequences
    6) Validate the annotations via alignment and/or RNA structure
        comparison to a library of true 16s rRNA genes
    7) Validate the length and percent of ambiguous bases in the rRNA gene to:
        /dedup/1200bp/
    8) Segregate out 16s rRNA genes from type strains to:
        /dedup/1200bp/types/
    9) Regex based parsing of taxonomic annotations to find out named, to:
        /dedup/1200bp/named/
    10) Cluster and remove outliers to generate:
        /dedup/1200bp/named/filtered/
    11) Add in 'trusted' seqs (manually curated versions in a text file.)
        /dedup/1200bp/named/filtered/trusted/
    12) Append the new entries to the extant library
    13) Use git to version and commit the changes.
*/


// User params initialization
params.help = false
params.testing = true


params.repo = './ya16sdb'
params.out = './refreshed'

params.email = false
params.ncbi_concurrent_connections = 3
params.retry_max = 1
params.retry_delay = 60000
params.min_len = 500
params.api_key = false
params.debug = false
// maximum number of 16S rRNA to return for a given species
params.species_cap = 5000

params.min_seqs_for_filtering = 5


if (params.debug == true){
    params.rRNA16S_bact_search = """16s[All Fields] AND rRNA[Feature Key] AND Bacteria[Organism] AND ${params.min_len} : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism]) AND sequence_from_type[Filter] AND ("2019/06/01"[Publication Date] : "2019/06/02"[Publication Date])"""
    params.rRNA16S_arch_search = """16s[All Fields] AND rRNA[Feature Key] AND Archaea[Organism] AND ${params.min_len} : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism]) AND ("2019/06/01"[Publication Date] : "2019/06/02"[Publication Date])"""
    params.rRNA16S_type_search = """16s[All Fields] AND rRNA[Feature Key] AND (Bacteria[Organism] OR  Archaea[Organism]) AND ${params.min_len} : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism]) AND sequence_from_type[Filter] AND ("2019/06/01"[Publication Date] : "2019/06/02"[Publication Date])"""

}
else {
    params.rRNA16S_bact_search = """16s[All Fields] AND rRNA[Feature Key] AND Bacteria[Organism] AND ${params.min_len} : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism])"""
    params.rRNA16S_arch_search = """16s[All Fields] AND rRNA[Feature Key] AND Archaea[Organism] AND ${params.min_len} : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism])"""
    params.rRNA16S_type_search = """16s[All Fields] AND rRNA[Feature Key] AND (Bacteria[Organism] OR  Archaea[Organism]) AND ${params.min_len} : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism]) AND sequence_from_type[Filter]"""
}

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run jgolob/ya16sdb <ARGUMENTS>
    
    Required Arguments:
    --repo                          path to directory holding the current repo (default = './ya16sdb')
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

    """.stripIndent()
}

// parameter validation

/*
    Repo structure:
    ./seqs.fasta  <- *all* 16S rRNA to date, valid, invalid, etc
    ./seq_info.csv <- Extract seq info, including version / accession
    ./refseq_info.csv <- If seqs are from refseq, put info here
    ./pubmed_info.csv <- Pubmed sources for the sequences
    ./references.csv  <- references for the sequences
    ./records.txt   <- All record IDs that meet our criteria

    ./dedup/1200bp/ DIR: holds full length 16S rRNA with limited ambiguous bases
    ./dedup/1200bp/seqs.fasta
    ./dedup/1200bp/seq_info.csv
    
    ./dedup/1200bp/types/ DIR: For 1200bp / type strain sourced 16S rRNA
    ./dedup/1200bp/types/seqs.fasta
    ./dedup/1200bp/types/seq_info.csv
    
    ./dedup/1200bp/named/ DIR: for 'named' 16sRNA (with a taxonomic annotation)
    ./dedup/1200bp/named/seqs.fasta
    ./dedup/1200bp/named/seq_info.csv
    ./dedup/1200bp/named/tax_id_map.csv 
    ./dedup/1200bp/named/taxonomy.csv <- holds taxonomic info for each gene

    ./dedup/1200bp/named/filtered/ DIR: for verified (filtered) taxonomic annotations 
    ./dedup/1200bp/named/filtered/seqs.fasta
    ./dedup/1200bp/named/filtered/seq_info.csv
    ./dedup/1200bp/named/filtered/tax_id_map.csv 
    ./dedup/1200bp/named/filtered/taxonomy.csv <- holds taxonomic info for each gene
*/
def paramInvalid() {
    if (
        (!file("${params.repo}").exists()) || 
        (!file("${params.repo}").isDirectory())
    ) { 
        log.error "Repo path ${params.repo} doesn't exist or isn't a directory.";
        return true;
    }
    if (params.email == false || params.email == null){
        log.error "You must supply a valid email.";
        return true;
    }
    
    // implicit else
    return false
}

if (params.help || paramInvalid()){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Step 0: Load in prior files from the repo
prior_seqs_f = file "${params.repo}/seqs.fasta"
prior_si_f = file "${params.repo}/seq_info.csv"
prior_records_f = file "${params.repo}/records.txt"
prior_unknowns_f = file "${params.repo}/unknowns.txt"
prior_pubmed_f = file "${params.repo}/pubmed_info.csv"
prior_refseq_info_f = file "${params.repo}/refseq_info.csv"
prior_references_f = file "${params.repo}/references.csv"
prior_refseqinfo_f = file "${params.repo}/refseq_info.csv"
prior_outliers_f = file "${params.repo}/dedup/1200bp/named/filtered/outliers.csv"

// Step 1: Retrieve current accessions with a 16S rRNA
// Archaea
process retrieveAcc_archaea {
    container 'golob/medirect:0.14.0__bcw.0.3.1B'
    label 'io_limited'
    errorStrategy 'finish'

    output:
        file "archaea_acc.txt" into acc_archaea_f

    """
    set -e

    ncbi_get_nt_accessions_for_query \
    --email ${params.email} \
    --ncbi_concurrent_connections ${params.ncbi_concurrent_connections} \
    --retry_max ${params.retry_max} --retry_delay ${params.retry_delay} \
    --query "${params.rRNA16S_arch_search}" \
    --out archaea_acc.txt
    """
}
// Bacteria
process retrieveAcc_bacteria {
    container 'golob/medirect:0.14.0__bcw.0.3.1B'
    label 'io_limited'
    errorStrategy 'finish'

    output:
        file "bacteria_acc.txt" into acc_bacteria_f

    """
    set -e

    ncbi_get_nt_accessions_for_query \
    --email ${params.email} \
    --ncbi_concurrent_connections ${params.ncbi_concurrent_connections} \
    --retry_max ${params.retry_max} --retry_delay ${params.retry_delay} \
    --query "${params.rRNA16S_bact_search}" \
    --out bacteria_acc.txt
    """
}

// From type strain organisms
process retrieveAcc_types {
    container 'golob/medirect:0.14.0__bcw.0.3.1B'
    label 'io_limited'
    errorStrategy 'finish'

    output:
        file "types_acc.txt" into acc_types_f

    """
    set -e

    ncbi_get_nt_accessions_for_query \
    --email ${params.email} \
    --ncbi_concurrent_connections ${params.ncbi_concurrent_connections} \
    --retry_max ${params.retry_max} --retry_delay ${params.retry_delay} \
    --query "${params.rRNA16S_type_search}" \
    --out types_acc.txt
    """
}

// TODO check for modified records by extracting modified dates


// Step 2. Use set adventures to figure out the records we have not looked at
// Can also include here masked records to remove (known empty, etc)
process combineCurrentAccessions {
    container 'golob/medirect:0.14.0__bcw.0.3.1B'
    label 'io_limited'
    errorStrategy 'finish'

    input:
        file acc_types_f
        file acc_archaea_f
        file acc_bacteria_f
    output:
        file "records.current.txt" into records_dl_f
"""
#!/usr/bin/env python
archaea_ver = {
    l.strip()
    for l in open("${acc_archaea_f}", 'rt')
}
bacteria_ver = {
    l.strip()
    for l in open("${acc_bacteria_f}", 'rt')
}
type_ver = {
    l.strip()
    for l in open("${acc_types_f}", 'rt')
}
with open('records.current.txt', 'wt') as out_h:
    for acc in sorted(archaea_ver.union(bacteria_ver).union(type_ver)):
        out_h.write(acc+"\\n")
"""    


}

process accessionsToDownload {
    container 'golob/medirect:0.14.0__bcw.0.3.1B'
    label 'io_limited'
    errorStrategy 'finish'
    cache 'deep'

    input:
        file records_dl_f
        file prior_records_f

    output:
        file "download.txt" into acc_to_download_f
script:
if (params.debug == false)
"""
#!/usr/bin/env python
current_version = {
    l.strip()
    for l in open("${records_dl_f}")
}
prior_version = {
    l.strip()
    for l in open("${prior_records_f}")
}
new_versions = sorted(current_version - prior_version)
with open('download.txt', 'wt') as out_h:
    for acc in new_versions:
        out_h.write(acc+"\\n")
"""
else
"""
#!/usr/bin/env python
current_version = {
    l.strip()
    for l in open("${records_dl_f}")
}
prior_version = {
    l.strip()
    for l in open("${prior_records_f}")
}
new_versions = sorted(current_version - prior_version)[0:1000]
with open('download.txt', 'wt') as out_h:
    for acc in new_versions:
        out_h.write(acc+"\\n")
"""
}


// Step 3. For each new accession, find the 16S rRNA features

process get16SrRNA_feat {
    container 'golob/medirect:0.14.0__bcw.0.3.1B'
    label 'io_limited'
    errorStrategy 'finish'

    input:
        file acc_to_download_f

    output:
        file "new_16s_rrna_feat.csv" into new_16s_rRNA_feat_f

    script:
    if (params.api_key == false)
    """
    set -e

    mefetch -proc ${task.cpus} -max-retry ${params.retry_max} -retry ${params.retry_delay} \
    --email ${params.email} -db nucleotide -mode text -format ft -id ${acc_to_download_f} |
    ftract -feature "rrna:product:16S ribosomal RNA" -min-length ${params.min_len} -on-error continue \
    -out new_16s_rrna_feat.csv
    """
    else
    """
    set -e

    mefetch -proc ${task.cpus} -max-retry ${params.retry_max} -retry ${params.retry_delay} \
    --email ${params.email} -api-key ${params.api_key} -db nucleotide -mode text \
    -format ft -id ${acc_to_download_f} |
    ftract -feature "rrna:product:16S ribosomal RNA" -min-length ${params.min_len} -on-error continue \
    -out new_16s_rrna_feat.csv
    """   
}

// Step 4: Download the new 16S rRNA features into genbank format and extract directly
today = new Date().format('dd-MMM-yyyy')

process get16SrRNA_gb {
    container 'golob/medirect:0.14.0__bcw.0.3.1B'
    label 'io_limited'
    errorStrategy 'finish'


    input:
        file new_16s_rRNA_feat_f
        val today
    
    output:
        file "seqs.fasta" into new_16s_seqs_fasta_f
        file "seq_info.csv" into new_16s_si_f
        file "pubmed_info.csv" into new_16s_pubmed_f
        file "references.csv" into new_16s_refs_f
        file "refseq_info.csv" into new_16s_refseqinfo_f

    script:
    if (params.api_key == false)
    """
    set -e

    mefetch -proc ${task.cpus} -max-retry ${params.retry_max} -retry ${params.retry_delay} \
    --email ${params.email} -db nucleotide -mode text \
    -csv -format gbwithparts -id ${new_16s_rRNA_feat_f} |
    extract_genbank ${today} \
    seqs.fasta \
    seq_info.csv \
    pubmed_info.csv \
    references.csv \
    refseq_info.csv
    """
    else
    """
    set -e

    mefetch -proc ${task.cpus} -max-retry ${params.retry_max} -retry ${params.retry_delay} \
    --email ${params.email} --api-key ${params.api_key} -db nucleotide -mode text \
    -csv -format gbwithparts -id ${new_16s_rRNA_feat_f} |
    extract_genbank ${today} \
    seqs.fasta \
    seq_info.csv \
    pubmed_info.csv \
    references.csv \
    refseq_info.csv
    """
}

// Step 5: Align against RDP type strains with vsearch to validate
process vsearch_rdp_validate {
    container 'golob/ya16sdb:0.2C'
    label 'mem_veryhigh'
    cache 'deep'
    errorStrategy 'finish'

    input:
        file new_16s_seqs_fasta_f
        file new_16s_si_f
        file prior_unknowns_f
        
    
    output:
        file "vsearch/seqs.fasta" into vsearch_seqs_f
        file "vsearch/seq_info.csv" into vsearch_si_f
        file "vsearch/unknown.fasta" into vsearch_unknown_fasta_f
        file "vsearch/unknowns.txt" into vsearch_unknowns_txt_f

    script:
    """
    set -e
    touch ${prior_unknowns_f}

    vsearch --threads ${task.cpus} \
    --db /db/rdp_16s_type_strains.fasta.gz \
    --id 0.70 --iddef 2 --mincols 350 --query_cov 0.70 --strand both  \
    --maxaccepts 1 --maxrejects 32 --top_hits_only \
    --output_no_hits --userfields query+target+qstrand+id+tilo+tihi \
    --usearch_global ${new_16s_seqs_fasta_f} \
    --userout vsearch.tsv
    mkdir -p vsearch
    vsearch.py vsearch.tsv ${new_16s_seqs_fasta_f} ${new_16s_si_f} ${prior_unknowns_f} \
    vsearch/seqs.fasta vsearch/seq_info.csv vsearch/unknown.fasta vsearch/unknowns.txt
    """
}

// Step 6: Build taxonomy DB Check tax IDs in seq_info using taxtastic
process downloadTaxdump {
    container = 'golob/ya16sdb:0.2C'
    label = 'io_limited'
    errorStrategy 'finish'
    publishDir path: "${params.out}/", mode: "copy"

    output: 
        file "taxdmp.zip" into taxdmp_f
    """
    set -e

    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
    """    
}

process buildTaxtasticDB {
    container = 'golob/ya16sdb:0.2C'
    label = 'io_limited'
    publishDir path: "${params.out}/", mode: "copy"
    errorStrategy 'finish'

    input:
        file taxdmp_f

    output:
        file "taxonomy.db" into taxonomy_db_f

    """
    taxit new_database taxonomy.db -z ${taxdmp_f}
    """
}

process filterUnknownTaxa {
    container = 'golob/ya16sdb:0.2C'
    label = 'io_limited'
    errorStrategy 'finish'

    input:
        file taxonomy_db_f
        file vsearch_si_f
    
    output:
        file "filtered/seq_info.csv" into new_filtered_si_f

    """
    mkdir -p filtered/
    taxit update_taxids \
    --unknown-action drop --outfile filtered/seq_info.csv \
    ${vsearch_si_f} ${taxonomy_db_f}
    """
}


// Step 7: Refresh the repo seqs!

// A bit of fakery for debug
if (args.debug != false) {
process debug_records {
    container 'golob/ya16sdb:0.2C'
    label 'io_limited'
    errorStrategy 'finish'
    cache 'deep'
    
    input:
        file "prior_records.txt" from prior_records_f
        file "current_records.txt" from records_dl_f
    output:
        file 'records.txt' into records_current_f
"""
#!/usr/bin/env python3

n = 0
with open('records.txt', 'wt') as out_h:
    for l in open('current_records.txt', 'rt'):
        out_h.write(l.strip()+'\\n')
    for l in open('prior_records.txt', 'rt'):
        out_h.write(l.strip()+'\\n')
        n += 1
        if n > 1000:
            break
"""
}} else {
    records_current_f = records_dl_f
}

process refreshRecords {
    container 'golob/ya16sdb:0.2C'
    label 'io_mem'
    cache 'deep'
    publishDir path: "${params.out}/", mode: "copy"
    errorStrategy 'finish'

    input:
        file "new/records.txt" from records_current_f
        file "new/seqs.fasta" from vsearch_seqs_f
        file "prior/seqs.fasta" from prior_seqs_f
        file "new/seq_info.csv" from new_filtered_si_f
        file "prior/seq_info.csv" from prior_si_f
        file "new/pubmed_info.csv" from new_16s_pubmed_f
        file "prior/pubmed_info.csv" from prior_pubmed_f
        file "new/references.csv" from new_16s_refs_f
        file "prior/references.csv" from prior_references_f
        file "new/refseq_info.csv" from new_16s_refseqinfo_f
        file "prior/refseq_info.csv" from prior_refseqinfo_f
        file vsearch_unknowns_txt_f
        file "prior/records.txt" from prior_records_f

    output:
        file "seqs.fasta" into refresh_seqs_f
        file "seq_info.csv" into refresh_si_unverified_f
        file "pubmed_info.csv" into refresh_pubmed_f
        file "references.csv" into refresh_references_f
        file "refseq_info.csv" into refresh_refseqinfo_f
        file "records.txt" into refresh_records_f

    script:
    """
    set -e

    touch prior/seqs.fasta
    touch prior/seq_info.csv
    touch prior/pubmed_info.csv
    touch prior/references.csv
    touch prior/refseq_info.csv
    touch prior/records.txt
    mkdir -p refresh/

    refresh.py \
    new/records.txt \
    new/seqs.fasta prior/seqs.fasta \
    new/seq_info.csv prior/seq_info.csv \
    new/pubmed_info.csv prior/pubmed_info.csv \
    new/references.csv prior/references.csv \
    new/refseq_info.csv prior/refseq_info.csv \
    ${vsearch_unknowns_txt_f} prior/records.txt \
    seqs.fasta \
    seq_info.csv \
    pubmed_info.csv \
    references.csv \
    refseq_info.csv \
    records.txt
    """
}

// Step 8: Verify the refreshed tax ID and make a taxtable
process refresh_verifyTaxIds {
    container 'golob/ya16sdb:0.2C'
    label 'io_mem'
    errorStrategy 'finish'

    input:
        file taxonomy_db_f
        file "unverified_si.csv" from refresh_si_unverified_f
    
    output:
        file "seq_info.csv" into refresh_si_f
    
    """
    taxit update_taxids \
    --unknown-action drop --outfile seq_info.csv \
    unverified_si.csv ${taxonomy_db_f}
    """
}

process taxonomyTable_refresh {
    container 'golob/ya16sdb:0.2C'
    label 'io_mem'
    publishDir path: "${params.out}/", mode: "copy"

    input:
        file taxonomy_db_f
        file refresh_si_f
    
    output:
        file "taxonomy.csv" into refresh_taxonomy_table_f
    
    """
    taxit -v taxtable --seq-info ${refresh_si_f} --out taxonomy.csv ${taxonomy_db_f}
    """
}

// Step 9: Make a feather file to help sort items around
process buildFeatherSI {
    container 'golob/ya16sdb:0.2C'
    label 'io_mem'
    errorStrategy 'finish'

    input:
        file refresh_seqs_f
        file refresh_si_f
        file refresh_taxonomy_table_f
        file acc_types_f
        file refresh_pubmed_f
        file refresh_refseqinfo_f
        file taxonomy_db_f

    output:
        file "seq_info.feather" into refresh_feather_si

    """
    to_feather.py ${refresh_si_f} seq_info.feather
    taxonomy.py seq_info.feather ${refresh_taxonomy_table_f}
    is_type.py seq_info.feather ${acc_types_f}
    is_published.py seq_info.feather ${refresh_pubmed_f}
    is_refseq.py seq_info.feather ${refresh_refseqinfo_f}
    is_valid.py seq_info.feather sqlite:///${taxonomy_db_f}
    confidence.py seq_info.feather
    seqhash.py seq_info.feather ${refresh_seqs_f}
    sort_values.py seq_info.feather "is_type,is_published,is_refseq,ambig_count,modified_date,download_date,seqhash"
    """
}

// Step 10: Split out our deduplicated 1200bp seqs.
process refresh_dd1200bp {
    container 'golob/ya16sdb:0.2C'
    label 'io_mem'
    publishDir path: "${params.out}/dedup/1200bp/", mode: "copy"
    errorStrategy 'finish'

    input:
        file refresh_feather_si
        file "unfiltered_seqs.fasta" from refresh_seqs_f

    output:
        file "seqs.fasta" into refresh_dd1200_seqs_f
        file "seq_info.csv" into refresh_dd1200_si_f
        file "blast.nhr" into refresh_dd1200_nhr_f
        file "blast.nsq" into refresh_dd1200_nsq_f
        file "blast.nin" into refresh_dd1200_nin_f

    """
    set -e

    partition_refs.py \
    --drop-duplicate-sequences \
    --min-length 1200 \
    --prop-ambig-cutoff 0.01 \
    unfiltered_seqs.fasta ${refresh_feather_si} \
    seqs.fasta seq_info.csv

    makeblastdb -dbtype nucl -in seqs.fasta -out blast
    """
}

process refresh_types {
    container 'golob/ya16sdb:0.2C'
    label 'io_mem'
    publishDir path: "${params.out}/dedup/1200bp/types/", mode: "copy"
    errorStrategy 'finish'

    input:
        file refresh_feather_si
        file "unfiltered_seqs.fasta" from refresh_dd1200_seqs_f
        file taxonomy_db_f

    output:
        file "seqs.fasta" into refresh_dd1200_types_seqs_f
        file "seq_info.csv" into refresh_dd1200_types_si_f
        file "taxonomy.csv" into refresh_dd1200_types_taxonomy_f
        file "lineages.csv" into refresh_dd1200_types_lineages_f
        file "lineages.txt" into refresh_dd1200_types_lineages_mothur_f
        file "blast.nhr" into refresh_dd1200_types_nhr_f
        file "blast.nsq" into refresh_dd1200_types_nsq_f
        file "blast.nin" into refresh_dd1200_types_nin_f

    """
    set -e
    mkdir -p partition
    partition_refs.py \
    --is_species \
    --is_type \
    --is_valid \
    unfiltered_seqs.fasta ${refresh_feather_si} \
    seqs.fasta seq_info.csv
    
    taxit -v taxtable --seq-info seq_info.csv --out taxonomy.csv ${taxonomy_db_f}
    taxit lineage_table --csv-table lineages.csv taxonomy.csv seq_info.csv
    taxit lineage_table --taxonomy-table lineages.txt taxonomy.csv seq_info.csv
    makeblastdb -dbtype nucl -in seqs.fasta -out blast
    """
}

process refresh_named {
    container 'golob/ya16sdb:0.2C'
    label 'io_mem'
    publishDir path: "${params.out}/dedup/1200bp/named/", mode: "copy"
    errorStrategy 'finish'

    input:
        file refresh_feather_si
        file "unfiltered_seqs.fasta" from refresh_dd1200_seqs_f
        file taxonomy_db_f

    output:
        file "seqs.fasta" into refresh_dd1200_named_seqs_f
        file "seq_info.csv" into refresh_dd1200_named_si_f
        file "taxonomy.csv" into refresh_dd1200_named_taxonomy_f
        file "lineages.csv" into refresh_dd1200_named_lineages_f
        file "lineages.txt" into refresh_dd1200_named_lineages_mothur_f
        file "blast.nhr" into refresh_dd1200_named_nhr_f
        file "blast.nsq" into refresh_dd1200_named_nsq_f
        file "blast.nin" into refresh_dd1200_named_nin_f

    """
    set -e

    partition_refs.py \
    --is_species \
    --is_valid \
    --species-cap ${params.species_cap} \
    unfiltered_seqs.fasta ${refresh_feather_si} \
    seqs.fasta seq_info.csv
    
    taxit -v taxtable --seq-info seq_info.csv --out taxonomy.csv ${taxonomy_db_f}
    taxit lineage_table --csv-table lineages.csv taxonomy.csv seq_info.csv
    taxit lineage_table --taxonomy-table lineages.txt taxonomy.csv seq_info.csv
    makeblastdb -dbtype nucl -in seqs.fasta -out blast
    """
}

// Group the named seqs by annotated tax id to create a channel
process groupNamedByTaxa {
    container 'golob/ya16sdb:0.2C'
    label 'io_limited'
    errorStrategy 'finish'

    input:
        file refresh_dd1200_named_si_f
    
    output:
        stdout into named_tax_groups_out

    """
    set -e

    csvcut.py --columns seqname,tax_id --out /dev/stdout ${refresh_dd1200_named_si_f}
    """
}
// convert the stdout into a channel
named_tax_groups_out
    .splitCsv(header: true)
    .map{r-> [r.tax_id, r.seqname]}
    .groupTuple()
    .set { named_tax_group_map }

// Load in the prior outliers.csv
Channel.from(prior_outliers_f)
    .splitCsv(header: true)
    .map { r -> [
        r.seqname,
        r.tax_id,
        r.centroid,
        r.cluster,
        r.dist,
        r.is_out,
        r.species,
        r.x,
        r.y
    ]}
    .tap { prior_outliers_ch }
    .map { r -> [
        r[1],   // tax_id
        r[0]    // seqname
    ]}
    .groupTuple()
    .set {
        prior_outliers_map
    }

// Use the prior outliers map to see what changed....
named_tax_group_map.join(prior_outliers_map, remainder: true)
    .filter { r -> r[1] != null}
    .into {
        tax_group_changed_ch;
        tax_group_unchanged_ch;
    }
tax_group_changed_ch
    .filter { r -> r[1] != r[2]}
    .map{ r -> [r[0], r[1]] }
    .set { 
        tax_group_changed_ch
    }
tax_group_unchanged_ch
    .filter { r -> r[1] == r[2]}
    .map{ r -> r[1] }
    .flatten()
    .set { 
        seqname_unchanged_ch
    }

// Of the changed, separate into 'rare' taxa (less than min reps for filtering) and non-rare
tax_group_changed_ch
    .into{
        tax_group_changed_rare_ch;
        tax_group_changed_tofilter_ch
    }
// for the changed rare we want to ungroup into seqname, tax_id tuples
tax_group_changed_rare_ch
    .filter {
        r -> r[1].size() < params.min_seqs_for_filtering
    }
    .flatMap{
            r -> 
            fl = [];
            r[1].eachWithIndex{ 
                it, i ->  fl.add([
                    r[1][i], // seqname
                    r[0]
                ])
            }
            return fl;
    }
    .toSortedList({a, b -> a[0] <=> b[0]})
    .flatMap()
    .set {
        seqname_changed_rare_ch
    }

tax_group_changed_tofilter_ch
    .filter {
        r -> r[1].size() >= params.min_seqs_for_filtering
    }
    .set {
        tax_group_changed_tofilter_ch
    }    




// Use this channel to subset out the named seqs and si
process taxonGroupFiles {
    container 'golob/ya16sdb:0.2C'
    label 'io_limited'
    errorStrategy 'finish'

    input:
        set tax_id, seq_names from tax_group_changed_tofilter_ch
        file refresh_dd1200_named_seqs_f
    
    output:
        set file('taxon_seqs.fasta'), file('taxon_seq_info.csv') into taxon_group_files_ch

"""
#!/usr/bin/env python3
from Bio import SeqIO
import csv

seq_names = {
    sn.strip() for sn in
    "${seq_names}".replace("[", "").replace("]", "").split(",")
}
num_seqs_out = 0
with open('taxon_seqs.fasta', 'wt') as seqs_out_h:
    for sr in SeqIO.parse(open('${refresh_dd1200_named_seqs_f}', 'rt'), 'fasta'):
        if sr.id in seq_names:
            SeqIO.write(sr, seqs_out_h, 'fasta')
            num_seqs_out += 1
assert (len(seq_names) == num_seqs_out)
with open('taxon_seq_info.csv', 'wt') as taxon_si_h:
    writer = csv.writer(taxon_si_h)
    writer.writerow(['seqname', 'tax_id'])
    for sn in seq_names:
        writer.writerow([sn, '${tax_id}'])
"""
}

// Run this channel through deenurp.
process filterOutliers {
    container 'golob/deenurp:0.2.6'
    label 'multithread'
    errorStrategy 'finish'

    input:
        set file(seqs_f), file(seq_info_f) from taxon_group_files_ch
        file refresh_dd1200_named_taxonomy_f
    
    output:
        set file('outliers.csv'), file('unsorted.fasta'), file('deenurp.log') into filter_outliers_ch 
    
    """
    set -e 

    deenurp -vvv filter_outliers \
    --strategy cluster \
    --cluster-type single \
    --distance-percentile 90.0 \
    --filter-rank species \
    --jobs 1 --threads-per-job ${task.cpus} \
    --max-distance 0.02 \
    --min-distance 0.01 \
    --min-seqs-for-filtering 5 \
    --log deenurp.log \
    --detailed-seqinfo outliers.csv \
    --output-seqs unsorted.fasta \
    ${seqs_f} ${seq_info_f} ${refresh_dd1200_named_taxonomy_f}
    """
}

// Make our new outliers.csv

// Use join to subset the prior outliers to the unchanged
prior_outliers_ch.join(seqname_unchanged_ch)
    .set{
        refresh_outliers_unchanged_ch
    }

// Rare changed set we let pass through as presumed valid
seqname_changed_rare_ch
    .map{ r -> [
        r[0], // seqname
        r[1], // tax_id
        "", // centroid
        "", // cluster
        "", // dist,
        "False", // is_out (false)
        r[1], // species = tax_id
        "", // x 
        "", // y
    ]}
    .set {
        refresh_outliers_changed_rare_ch
    }

filter_outliers_ch
    // read in each outliers.csv. Remove header. 
    .flatMap { r -> 
        lines = r[0].readLines();
        lines.removeAt(0);
        return lines;
    }
    .set {
        refresh_outliers_changed_ch
    }


refresh_outliers_unchanged_ch
    .mix(
        refresh_outliers_changed_rare_ch
    )
    .map { r -> 
        r.join(',')
    }
    .mix(
        refresh_outliers_changed_ch
    )
    .reduce(""){ p,c ->
        p += (c+"\\n");
        return p;
    }
    .set{
        refresh_outlier_details_val
    }


process makeRefreshedOutliers {
    container 'golob/ya16sdb:0.2C'
    label 'io_limited'
    errorStrategy 'finish'
    publishDir path: "${params.out}/dedup/1200bp/named/filtered/", mode: 'copy'

    input:
        val refresh_outlier_details_val
    output:
        file 'outliers.csv' into refresh_outliers_f
    
    """
    printf "seqname,tax_id,centroid,cluster,dist,is_out,species,x,y\\n" > outliers.csv
    printf "${refresh_outlier_details_val}" >> outliers.csv
    """
}

// Integrate in the outlier information to the feather-seq-info

process injectOutlier {
    container 'golob/ya16sdb:0.2C'
    label 'io_mem'
    errorStrategy 'finish'

    input:
        file refresh_feather_si
        file refresh_outliers_f
    
    output:
        file 'seq_info_wOutlier.feather' into refresh_feather_wOutlier_si
    """
    set -e

    cp ${refresh_feather_si} seq_info_wOutlier.feather
    filter_outliers.py seq_info_wOutlier.feather ${refresh_outliers_f}
    """
}

process refresh_filtered {
    container 'golob/ya16sdb:0.2C'
    label 'io_mem'
    errorStrategy 'finish'
    publishDir path: "${params.out}/dedup/1200bp/named/filtered/", mode: "copy"

    input:
        file refresh_feather_wOutlier_si
        file "unfiltered_seqs.fasta" from refresh_dd1200_named_seqs_f
        file taxonomy_db_f

    output:
        file "seqs.fasta" into refresh_dd1200_named_filtered_seqs_f
        file "seq_info.csv" into refresh_dd1200_named_filtered_si_f
        file "taxonomy.csv" into refresh_dd1200_named_filtered_taxonomy_f
        file "lineages.csv" into refresh_dd1200_named_filtered_lineages_f
        file "lineages.txt" into refresh_dd1200_named_filtered_lineages_mothur_f
        file "blast.nhr" into refresh_dd1200_named_filtered_nhr_f
        file "blast.nsq" into refresh_dd1200_named_filtered_nsq_f
        file "blast.nin" into refresh_dd1200_named_filtered_nin_f

    """
    set -e

    partition_refs.py \
    --inliers \
    --is_species \
    --is_valid \
    --species-cap ${params.species_cap} \
    unfiltered_seqs.fasta ${refresh_feather_wOutlier_si} \
    seqs.fasta seq_info.csv
    
    taxit -v taxtable --seq-info seq_info.csv --out taxonomy.csv ${taxonomy_db_f}
    taxit lineage_table --csv-table lineages.csv taxonomy.csv seq_info.csv
    taxit lineage_table --taxonomy-table lineages.txt taxonomy.csv seq_info.csv
    makeblastdb -dbtype nucl -in seqs.fasta -out blast
    """
}

/*


// */