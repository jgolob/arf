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
    7) Deduplicate the sequences to:
        /dedup/
    8) Validate the length and percent of ambiguous bases in the rRNA gene to:
        /dedup/1200bp/
    9) Segregate out 16s rRNA genes from type strains to:
        /dedup/1200bp/
    10) Regex based parsing of taxonomic annotations to find out named, to:
        /dedup/1200bp/named/
    11) Cluster and remove outliers to generate:
        /dedup/1200bp/named/filtered/
    12) Append the new entries to the extant library
    13) Use git to version and commit the changes.
*/


// User params initialization
params.help = false
params.testing = true


params.repo = './output'
params.email = false
params.ncbi_concurrent_connections = 3
params.retry_max = 1
params.retry_delay = 60000
params.min_len = 500
params.api_key = false
params.debug = false

if (params.debug == true){
    params.rRNA16S_bact_search = """16s[All Fields] AND rRNA[Feature Key] AND Bacteria[Organism] AND ${params.min_len} : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism]) AND sequence_from_type[Filter] AND ("2019/06/01"[Publication Date] : "2019/06/02"[Publication Date])"""
}
else {
    params.rRNA16S_bact_search = """16s[All Fields] AND rRNA[Feature Key] AND Bacteria[Organism] AND ${params.min_len} : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism])"""

}
params.rRNA16S_arch_search = """16s[All Fields] AND rRNA[Feature Key] AND Archaea[Organism] AND ${params.min_len} : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism])"""
params.rRNA16S_type_search = """16s[All Fields] AND rRNA[Feature Key] AND (Bacteria[Organism] OR  Archaea[Organism]) AND ${params.min_len} : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism]) AND sequence_from_type[Filter] AND ("2019/06/01"[Publication Date] : "2019/06/02"[Publication Date])"""

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run jgolob/ya16sdb <ARGUMENTS>
    
    Required Arguments:
    --repo                          path to directory holding the extant repo (default = './output')
    --email                         Valid email to use with NCBI
    --ncbi_concurrent_connections   Number of concurrent connections (default = 3)
    --retry_max                     Max retries with NCBI requests (default = 1)
    --retry_delay                   Delay (ms) between retries (default=60000)


    Options:
    --debug                         Cuts down the number of records for testing

    """.stripIndent()
}

// parameter validation

/*
    Repo structure:
    ./seqs.fasta  <- *all* 16S rRNA to date, valid, invalid, etc
    ./seq_info.csv <- Extract seq info, including version / accession
    ./refseq_info.csv <- CSV formatted information about the sequences
    ./pubmed_info.csv <- Pubmed sources for the sequences
    ./references.csv  <- references for the sequences
    
    ./dedup/ <- DIR: Holds deduplicated 16S rRNA
    ./dedup/seqs.fasta
    ./dedup/seq_info.csv
    ./dedup/taxonomy.csv <- raw taxonomic annotations
    
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

// Step 1: Retrieve current accessions with a 16S rRNA
// Archaea
process retrieveAcc_archaea {
    container 'golob/medirect:0.14.0__bcw.0.3.1B'
    label 'io_limited'
    errorStrategy 'retry'

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
    errorStrategy 'retry'

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
    errorStrategy 'retry'

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
    input:
        file acc_types_f
        file acc_archaea_f
        file acc_bacteria_f
    output:
        file "records.current.txt" into records_current_f
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
    for acc in archaea_ver.union(bacteria_ver).union(type_ver):
        out_h.write(acc+"\\n")
"""    


}

process accessionsToDownload {
    container 'golob/medirect:0.14.0__bcw.0.3.1B'
    label 'io_limited'
    //errorStrategy 'retry'

    input:
        file records_current_f
        file prior_records_f

    output:
        file "download.txt" into acc_to_download_f
script:
if (params.debug == false)
"""
#!/usr/bin/env python
current_version = {
    l.strip()
    for l in open("${records_current_f}")
}
prior_version = {
    l.strip()
    for l in open("${prior_records_f}")
}
new_versions = current_version - prior_version
with open('download.txt', 'wt') as out_h:
    for acc in new_versions:
        out_h.write(acc+"\\n")
"""
else
"""
#!/usr/bin/env python
current_version = {
    l.strip()
    for l in open("${records_current_f}")
}
prior_version = {
    l.strip()
    for l in open("${prior_records_f}")
}
new_versions = list(current_version - prior_version)[0:10]
with open('download.txt', 'wt') as out_h:
    for acc in new_versions:
        out_h.write(acc+"\\n")
"""
}


// Step 3. For each new accession, find the 16S rRNA features

process get16SrRNA_feat {
    container 'golob/medirect:0.14.0__bcw.0.3.1B'
    label 'io_limited'
    //errorStrategy 'retry'

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
    //errorStrategy 'retry'


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
    //errorStrategy 'retry'

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
process DlBuildTaxtasticDB {
    container = 'golob/ya16sdb:0.2B'
    label = 'io_limited'
    // errorStrategy = 'retry'

    output:
        file "taxonomy.db" into taxonomy_db_f

    afterScript "rm -rf dl/"

    """
    mkdir -p dl/ && \
    taxit new_database taxonomy.db -u ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip -p dl/
    """
}

process filterUnknownTaxa {
    container = 'golob/ya16sdb:0.2C'
    label = 'io_limited'
    // errorStrategy = 'retry'

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

process refreshRecords {
    container 'golob/ya16sdb:0.2C'
    label 'io_mem'
    //errorStrategy 'retry'

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
        file "refresh/seqs.fasta" into refresh_seqs_f
        file "refresh/seq_info.csv" into refresh_si_f
        file "refresh/pubmed_info.csv" into refresh_pubmed_f
        file "refresh/references.csv" into refresh_references_f
        file "refresh/refseq_info.csv" into refresh_refseqinfo_f
        file "refresh/records.txt" into refresh_records_f

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
    refresh/seqs.fasta \
    refresh/seq_info.csv \
    refresh/pubmed_info.csv \
    refresh/references.csv \
    refresh/refseq_info.csv \
    refresh/records.txt
    """
}

// Step 8: Make a (mothur-style) taxonomy table of the refreshed reads
process taxonomyTable_refresh {
    container 'golob/ya16sdb:0.2C'
    label 'io_mem'

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

/*


// */