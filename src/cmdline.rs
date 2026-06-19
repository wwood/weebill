use clap::{Args, Parser, Subcommand};
use crate::constants::*;

#[derive(Parser)]
#[clap(author, version, about = "Ultrafast genome ANI queries and taxonomic profiling for metagenomic shotgun samples.\n\n--- Preparing inputs by sketching (indexing)\n## fastq (reads) and fasta (genomes all at once)\n## *.sylsp found in -d; *.syldb given by -o\nsylph sketch -t 5 sample1.fq sample2.fq genome1.fa genome2.fa -o genome1+genome2 -d sample_dir\n\n## paired-end reads\nsylph sketch -1 a_1.fq b_1.fq -2 b_2.fq b_2.fq -d paired_sketches\n\n--- Nearest neighbour containment ANI\nsylph query *.syldb *.sylsp > all-to-all-query.tsv\n\n--- Taxonomic profiling with relative abundances and ANI\nsylph profile *.syldb *.sylsp > all-to-all-profile.tsv", arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    /// Sketch sequences into samples (reads) and databases (genomes). Each sample.fq -> sample.sylsp. All *.fa -> *.syldb. 
    #[clap(display_order = 1)]
    Sketch(SketchArgs),
    /// Coverage-adjusted ANI querying between databases and samples.
    #[clap(display_order = 3)]
    Query(ContainArgs),
    ///Species-level taxonomic profiling with abundances and ANIs. 
    #[clap(display_order = 2)]
    Profile(ContainArgs),
    ///Inspect sketched .syldb and .sylsp files.
    #[clap(arg_required_else_help = true, display_order = 4)]
    Inspect(InspectArgs),
    ///Build a k-mer dereplicated reference database (.sylref) for reference-delta sample compression.
    #[clap(arg_required_else_help = true, display_order = 5)]
    RefBuild(RefBuildArgs),
    ///Compress sample sketches against a reference DB (.sylsp -> .sylspr), or --decompress to reverse.
    #[clap(arg_required_else_help = true, display_order = 6)]
    RefCompress(RefCompressArgs),
    ///Convert a standard database (.syldb) into a two-stage seekable database (.syl2db) for `profile --two-stage`.
    #[clap(arg_required_else_help = true, display_order = 7)]
    DbConvert(DbConvertArgs),
}

#[derive(Args)]
pub struct DbConvertArgs {
    #[clap(multiple=true, help = "Standard genome database sketches (*.syldb) to convert")]
    pub files: Vec<String>,
    #[clap(short='o', long="output", help = "Output two-stage database name (.syl2db appended)")]
    pub output: String,
    #[clap(long="screen-c", default_value_t = SCREEN_C_DEFAULT, help = "Subsampling rate -c of the small in-memory stage-1 SCREEN index (the bincoded sparse hashes). Must be >= the database -c. A coarser (larger) value gives a smaller/faster screen index. The dense per-genome blocks always keep every k-mer at the database -c.")]
    pub screen_c: usize,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(long="trace", help = "Trace output")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output")]
    pub debug: bool,
}

#[derive(Args)]
pub struct RefBuildArgs {
    #[clap(multiple=true, help = "Genome database sketches (*.syldb) to build the reference from")]
    pub files: Vec<String>,
    #[clap(short='T', long="taxonomy", help_heading = "INPUT", help = "TSV with one line per genome: <genome_file_name><TAB><species><TAB><rep|strain>. The genome name matches the sketched path or its basename. Genomes absent from the file are treated as their own single-genome species representative. Strains of a species are placed contiguously, representatives first.")]
    pub taxonomy: Option<String>,
    #[clap(short='o', long="output", help = "Output reference database name (.sylref appended)")]
    pub output: String,
    #[clap(long="sparse-subsample", default_value_t = 16, help = "Stage-1 sparse index subsampling divisor: 1/N of each genome's distinctive k-mers are kept uncompressed for fast hit detection. Larger N gives a smaller/faster stage-1 index but coarser detection (1 keeps all).")]
    pub sparse_div: u64,
    #[clap(long="pool-min-genomes", default_value_t = 3, help = "Minimum number of same-tier genomes required before a k-mer is placed in the shared pool. With 3, k-mers shared by exactly two reps/strains are assigned to the first such genome instead of the pool.")]
    pub pool_min_genomes: u32,
    #[clap(long="max-ram", help = "Approximate peak RAM target (GB) for building. Sizes the number of on-disk partitions the build streams through; a soft target, not a hard limit.")]
    pub max_ram: Option<usize>,
    #[clap(long="tmp-dir", help = "Directory for build scratch files (needs roughly the input database size of free space). Default: alongside the output.")]
    pub tmp_dir: Option<String>,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(long="trace", help = "Trace output")]
    pub trace: bool,
}

#[derive(Args)]
pub struct RefCompressArgs {
    #[clap(multiple=true, help = "Sample sketches (*.sylsp) to compress, or (*.sylspr) with --decompress")]
    pub files: Vec<String>,
    #[clap(short='r', long="reference", help = "Reference database (*.sylref) produced by `sylph ref-build`")]
    pub ref_db: Option<String>,
    #[clap(long="decompress", help = "Reverse the operation: reconstruct *.sylsp from *.sylspr")]
    pub decompress: bool,
    #[clap(long="inspect", help = "Inspect reference-delta sketches (*.sylspr) and report metadata plus encoded section sizes")]
    pub inspect: bool,
    #[clap(long="verify", help = "Verify existing *.sylspr inputs by decompressing them and requiring exact equality to the original sketch path stored in each file")]
    pub verify: bool,
    #[clap(short='d', long="output-directory", default_value = "./", help = "Output directory")]
    pub output_dir: String,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(long="trace", help = "Trace output")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output")]
    pub debug: bool,
}


#[derive(Args, Default)]
pub struct SketchArgs {
    #[clap(multiple=true, help_heading = "INPUT", help = "fasta/fastq files; gzip optional. Default: fastq file produces a sample sketch (*.sylsp) while fasta files are combined into a database (*.syldb).")]
    pub files: Vec<String>,
    #[clap(short='o',long="out-name-db", default_value = "database", help_heading = "OUTPUT", help = "Output name for database sketch (with .syldb appended)")]
    pub db_out_name: String,
    #[clap(short='d',long="sample-output-directory", default_value = "./", help_heading = "OUTPUT", help = "Output directory for sample sketches")]
    pub sample_output_dir: String,
    #[clap(long="compressed-output", help_heading = "OUTPUT", help = "Like -o, but writes a compressed database sketch (with .syldbc appended). Compressed sketches are smaller on disk and readable by query/profile/inspect")]
    pub compressed_db_out_name: Option<String>,
    #[clap(long="compressed-database", help_heading = "OUTPUT", help = "Like -d, but writes compressed sample sketches (with .sylspc appended). Compressed sketches are smaller on disk and readable by query/profile/inspect")]
    pub compressed_sample_output_dir: Option<String>,
    #[clap(short,long="individual-records", help_heading = "GENOME INPUT", help = "Use individual records (contigs) for database construction")]
    pub individual: bool,
    #[clap(multiple=true,short,long="reads", help_heading = "SINGLE-END INPUT", help = "Single-end fasta/fastq reads")]
    pub reads: Option<Vec<String>>,
    #[clap(multiple=true,short='g', long="genomes", help_heading = "GENOME INPUT", help = "Genomes in fasta format")]
    pub genomes: Option<Vec<String>>,
    #[clap(short,long="list", help_heading = "INPUT", help = "Newline delimited file with inputs; fastas -> database, fastq -> sample")]
    pub list_sequence: Option<String>,
    #[clap(long="rl", hidden=true, help_heading = "SINGLE-END INPUT", help = "Newline delimited file; inputs assumed reads")]
    pub list_reads: Option<String>,
    #[clap(long="gl", help_heading = "GENOME INPUT", help = "Newline delimited file; inputs assumed genomes")]
    pub list_genomes: Option<String>,
    #[clap(long="l1", help_heading = "PAIRED-END INPUT", help = "Newline delimited file; inputs are first pair of PE reads")]
    pub list_first_pair: Option<String>,
    #[clap(long="l2", help_heading = "PAIRED-END INPUT", help = "Newline delimited file; inputs are second pair of PE reads")]
    pub list_second_pair: Option<String>,
    #[clap(long="lS", help_heading = "INPUT", help = "Newline delimited file; read sketches are renamed to given sample names")]
    pub list_sample_names: Option<String>,
    #[clap(multiple=true, short='S', long="sample-names", help_heading = "INPUT", help = "Read sketches are renamed to given sample names")]
    pub sample_names: Option<Vec<String>>,

    #[clap(short, default_value_t = 31,help_heading = "ALGORITHM", help ="Value of k. Only k = 21, 31 are currently supported")]
    pub k: usize,
    #[clap(short, default_value_t = 200, help_heading = "ALGORITHM", help = "Subsampling rate")]
    pub c: usize,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(long="ram-barrier", help = "Stop multi-threaded read sketching when (virtual) RAM is past this value (in GB). Does NOT guarantee max RAM limit", hidden=true)]
    pub max_ram: Option<usize>,
    #[clap(long="trace", help = "Trace output (caution: very verbose)")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output")]
    pub debug: bool,


    #[clap(long="no-dedup", help_heading = "ALGORITHM", help = "Disable read deduplication procedure. Reduces memory; not recommended for illumina data")]
    pub no_dedup: bool,
    #[clap(long="disable-profiling", help_heading = "ALGORITHM", help = "Disable sylph profile usage for databases; may decrease size and make sylph query slightly faster", hidden=true)]
    pub no_pseudotax: bool,
    #[clap(long="min-spacing", default_value_t = 30, help_heading = "ALGORITHM", help = "Minimum spacing between selected k-mers on the genomes")]
    pub min_spacing_kmer: usize,
    #[clap(long="fpr", default_value_t = DEFAULT_FPR, help_heading = "ALGORITHM", help = "False positive rate for read deduplicate hashing; valid values in [0,1).")]
    pub fpr: f64,
    #[clap(short='1',long="first-pairs", multiple=true, help_heading = "PAIRED-END INPUT", help = "First pairs for paired end reads")]
    pub first_pair: Vec<String>,
    #[clap(short='2',long="second-pairs", multiple=true, help_heading = "PAIRED-END INPUT", help = "Second pairs for paired end reads")]
    pub second_pair: Vec<String>,
}

#[derive(Args, Clone)]
pub struct ContainArgs {
    #[clap(multiple=true, help = "Pre-sketched *.syldb/*.sylsp files. Raw single-end fastq/fasta are allowed and will be automatically sketched to .sylsp/.syldb")]
    pub files: Vec<String>,

    
    #[clap(short='l',long="list", help = "Newline delimited file of file inputs",help_heading = "INPUT/OUTPUT")]
    pub file_list: Option<String>,

    #[clap(long,default_value_t = 3., help_heading = "ALGORITHM", help = "Minimum k-mer multiplicity needed for coverage correction. Higher values gives more precision but lower sensitivity")]
    pub min_count_correct: f64,
    #[clap(short='M',long,default_value_t = 50., help_heading = "ALGORITHM", help = "Exclude genomes with less than this number of sampled k-mers")]
    pub min_number_kmers: f64,
    #[clap(short, long="minimum-ani", help_heading = "ALGORITHM", help = "Minimum adjusted ANI to consider (0-100). Default is 90 for query and 95 for profile. Smaller than 95 for profile will give inaccurate results." )]
    pub minimum_ani: Option<f64>,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(short='s', long="sample-threads", help = "Number of samples to be processed concurrently. Default: (# of total threads / 3) + 1 for profile, 1 for query")]
    pub sample_threads: Option<usize>,
    #[clap(long="trace", help = "Trace output (caution: very verbose)")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output")]
    pub debug: bool,

    #[clap(long="estimate-read-counts", help_heading = "ALGORITHM", help = "Very roughly estimate read counts in the 'Sequence_abundance' column instead of relative abundance. This forces `-u`, which may have caveats for long reads and complex environments.")]
    pub estimate_read_counts: bool,

    #[clap(short='u', long="estimate-unknown", help_heading = "ALGORITHM", help = "Estimate true coverage and scale sequence abundance in `profile` by estimated unknown sequence percentage" )]
    pub estimate_unknown: bool,
    
    #[clap(short='I',long="read-seq-id", help_heading = "ALGORITHM", help = "Sequence identity (%) of reads. Only used in -u option and overrides automatic detection. ")]
    pub seq_id: Option<f64>,

    //#[clap(short='l', long="read-length", help_heading = "ALGORITHM", help = "Read length (single-end length for pairs). Only necessary for short-read coverages when using --estimate-unknown. Not needed for long-reads" )]
    //pub read_length: Option<usize>,

    #[clap(short='R', long="redundancy-threshold", help_heading = "ALGORITHM", help = "Removes redundant genomes up to a rough ANI percentile when profiling", default_value_t = 99.0, hidden=true)]
    pub redundant_ani: f64,

    #[clap(short='r',long="reads", multiple=true, help = "Single-end raw reads (fastx/gzip)", display_order = 1, help_heading = "SKETCHING")]
    pub reads: Vec<String>,

    #[clap(short='1', long="first-pairs", multiple=true, help = "First pairs for raw paired-end reads (fastx/gzip)", help_heading = "SKETCHING")]
    pub first_pair: Vec<String>,

    #[clap(short='2', long="second-pairs", multiple=true, help = "Second pairs for raw paired-end reads (fastx/gzip)", help_heading = "SKETCHING")]
    pub second_pair: Vec<String>,

    #[clap(short, default_value_t = 200, help_heading = "SKETCHING", help = "Subsampling rate. Does nothing for pre-sketched files")]
    pub c: usize,
    #[clap(short, default_value_t = 31, help_heading = "SKETCHING", help = "Value of k. Only k = 21, 31 are currently supported. Does nothing for pre-sketched files")]
    pub k: usize,
    #[clap(short,long="individual-records", help_heading = "SKETCHING", help = "Use individual records (e.g. contigs) for database construction instead. Does nothing for pre-sketched files")]
    pub individual: bool,
    #[clap(long="min-spacing", default_value_t = 30, help_heading = "SKETCHING", help = "Minimum spacing between selected k-mers on the database genomes. Does nothing for pre-sketched files")]
    pub min_spacing_kmer: usize,

    #[clap(short='o',long="output-file", help = "Output to this file (TSV format). [default: stdout]", help_heading="INPUT/OUTPUT")]
    pub out_file_name: Option<String>,
    #[clap(long="reference", help_heading="INPUT/OUTPUT", help = "Reference database (*.sylref from `sylph ref-build`) used to decode reference-delta compressed samples (*.sylspr). Required when any input is a *.sylspr file.")]
    pub reference: Option<String>,
    #[clap(long="log-reassignments", help = "Output information for how k-mers for genomes are reassigned during `profile`. Caution: can be verbose and slows down computation.")]
    pub log_reassignments: bool,

    #[clap(long="two-stage", help_heading = "TWO-STAGE PROFILING", help = "Two-stage profiling (profile only): cheaply SCREEN the sample against the (sparse) database, then densely profile ONLY the genomes that pass the screen. Lets a sparse pre-built database (e.g. -c 200 GTDB) deliver dense -c profiling without ever building/loading a dense full database.")]
    pub two_stage: bool,
    #[clap(long="dense-c", default_value_t = DENSE_C_DEFAULT, help_heading = "TWO-STAGE PROFILING", help = "Subsampling rate -c for the dense second stage. Genomes passing the screen are (re)sketched at this rate from their source fasta if the database is sparser than this. The sample sketch must have -c <= this value.")]
    pub dense_c: usize,
    #[clap(long="screen-c", help_heading = "TWO-STAGE PROFILING", help = "Subsampling rate -c for the cheap first-stage screen. Default: the database's own -c. Must be >= the database -c (a sketch can only be made sparser, never denser).")]
    pub screen_c: Option<usize>,
    #[clap(long="screen-ani", default_value_t = SCREEN_MIN_ANI_DEFAULT, help_heading = "TWO-STAGE PROFILING", help = "Minimum adjusted ANI (0-100) for a genome to pass the first-stage screen. Deliberately permissive; the dense stage recovers specificity.")]
    pub screen_ani: f64,
    #[clap(long="screen-min-matches", default_value_t = 1, help_heading = "TWO-STAGE PROFILING", help = "Minimum number of matched stage-1 screen k-mers for a genome to pass the screen and be densely decoded. Default 1 keeps the same results as single-stage; raising it (e.g. with a permissive --screen-ani) cheaply prunes genomes that pass on a handful of chance-shared k-mers, cutting wasted dense decodes at a small sensitivity cost for very-low-coverage genomes.")]
    pub screen_min_matches: usize,
    #[clap(long="dense-cache", help_heading = "TWO-STAGE PROFILING", help = "Directory of cached per-genome dense sketches (*.sylgn). Genomes (re)sketched for the dense stage are stored here and reused across samples/runs, so a dense database is grown lazily only for genomes that actually appear.")]
    pub dense_cache: Option<String>,
    #[clap(long="screen-dump", hidden=true, help_heading = "TWO-STAGE PROFILING", help = "Debug: write a TSV of every stage-1 screen survivor (genome, matched/total screen k-mers, naive/adjusted ANI, median coverage) to this file.")]
    pub screen_dump: Option<String>,


    //Hidden options that are embedded in the args but no longer used... 
    #[clap(short, hidden=true, long="pseudotax", help_heading = "ALGORITHM", help = "Pseudo taxonomic classification mode. This removes shared k-mers between species by assigning k-mers to the highest ANI species. Requires sketches with --enable-pseudotax option" )]
    pub pseudotax: bool,
    #[clap(long="ratio", hidden=true)]
    pub ratio: bool,
    #[clap(long="mme", hidden=true)]
    pub mme: bool,
    #[clap(long="mle", hidden=true)]
    pub mle: bool,
    #[clap(long="nb", hidden=true)]
    pub nb: bool,
    #[clap(long="no-ci", help = "Do not output confidence intervals", hidden=true)]
    pub no_ci: bool,
    #[clap(long="no-adjust", hidden=true)]
    pub no_adj: bool,
    #[clap(long="mean-coverage", help_heading = "ALGORITHM", help = "Use the robust mean coverage estimator instead of median estimator", hidden=true )]
    pub mean_coverage: bool,

}

#[derive(Args)]
pub struct InspectArgs {
    #[clap(multiple=true, help = "Pre-sketched *.syldb/*.sylsp files.")]
    pub files: Vec<String>,
    #[clap(short='o',long="output-file", help = "Output to this file (YAML format). [default: stdout]")]
    pub out_file_name: Option<String>,

}
    
