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

//params.rRNA16S_bact_search = "16s[All Fields] AND rRNA[Feature Key] AND Bacteria[Organism] AND 500 : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism])"
params.rRNA16S_bact_search = """16s[All Fields] AND rRNA[Feature Key] AND Bacteria[Organism] AND ("2019/06/01"[Publication Date] : "2019/06/02"[Publication Date]) AND ${params.min_len} : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism])"""
//params.rRNA16S_arch_search = "16s[All Fields] AND rRNA[Feature Key] AND Archaea[Organism] AND 500 : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism])"
params.rRNA16S_arch_search = """16s[All Fields] AND rRNA[Feature Key] AND Archaea[Organism] AND ("2019/06/01"[Publication Date] : "2019/06/30"[Publication Date]) AND ${params.min_len} : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism])"""


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

// Step 1a and 1b: Retrieve current accessions with a 16S rRNA (archaea and bacteria)
process retrieveAccArchaea {
    container 'golob/medirect:0.14.0__bcw.0.3.1A'
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

process retrieveAccBacteria {
    container 'golob/medirect:0.14.0__bcw.0.3.1A'
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

current_records_f = file "${params.repo}/records.txt"
process getNewAccessions {
    container 'golob/medirect:0.14.0__bcw.0.3.1A'
    label 'io_limited'
    //errorStrategy 'retry'

    input:
        file acc_archaea_f
        file acc_bacteria_f
        file current_records_f

    output:
        stdout into acc_to_download_out

    """
    #!/usr/bin/env python
    import csv

    current_version = {
        l.strip()
        for l in open("${current_records_f}")
    }
    
    archaea_ver = {
        l.strip()
        for l in open("${acc_archaea_f}", 'rt')
    }
    bacteria_ver = {
        l.strip()
        for l in open("${acc_bacteria_f}", 'rt')
    }
    new_versions = archaea_ver.union(bacteria_ver) - current_version
    for acc in new_versions:
        if acc is not None and acc is not "":
            print(acc)
    """
}
acc_to_download_out
    .splitText()
    .map{ r-> r.strip()}
    .filter{ s -> !s.isEmpty()}
    .first()
    .set { acc_to_download_ch }

// Step 2. For each new accession, collect a feature table

process get16SrRNA_feat {
    container 'golob/medirect:0.14.0__bcw.0.3.1A'
    label 'io_limited'
    //errorStrategy 'retry'
    maxForks 1

    input:
        val id from acc_to_download_ch

    output:
        stdout into ft_16srRNA_out

    """
    set -e

    mefetch -proc ${task.cpus} -max-retry ${params.retry_max} -retry ${params.retry_delay} \
    --email ${params.email} -db nucleotide -mode text -format ft -id ${id} |
    ftract -feature "rrna:product:16S ribosomal RNA" -min-length ${params.min_len} -on-error continue -out /dev/stdout
    """
}

ft_16srRNA_out
    .splitCsv(header: true)
    .map{ r->
        return [r.id, r.seq_start, r.seq_stop, r.strand, (r.id =~ /gb\|(\w+)\.\d+\|/).group(0)]
    }
    .set{ rRNA16s_for_dl_ch }

rRNA16s_for_dl_ch.println()

/*
// Step 3: Download genbank
process get16SrRNA {
    container 'golob/medirect:0.14.0__bcw.0.3.1A'
    label 'io_limited'
    //errorStrategy 'retry'
    maxForks 1

    input:
        val id, val seq_start, val seq_stop, val strand from rRNA16s_for_dl_ch
    
    output:
        val id, val seq_start, val seq_stop, val strand, file ("${}") into rRNA16s_dl_ch


    """
    set -e

    echo id,seq_start,seq_stop,strand
    echo ${id},${seq_start},${seq_stop},${strand}
    """

}

/*


// */