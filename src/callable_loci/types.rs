use serde::{Deserialize, Serialize};

/// The `CalledState` enumeration represents various states or classifications
/// for a callable entity, typically used in analysis pipelines, genomic contexts,
/// or scenarios requiring classification of coverage or mapping quality.
///
/// # Variants
///
/// - `REF_N`:
///   Indicates a reference 'N' state, typically used to represent unknown or ambiguous bases.
///
/// - `CALLABLE`:
///   Represents a state where the entity is callable, meeting all necessary criteria for analysis.
///
/// - `NO_COVERAGE`:
///   Indicates that there is no coverage for the entity, implying that no supporting data was observed.
///
/// - `LOW_COVERAGE`:
///   Represents a state where the coverage is low, potentially limiting confidence or accuracy of analysis.
///
/// - `EXCESSIVE_COVERAGE`:
///   Indicates that the coverage far exceeds expected thresholds, which might be indicative of artifacts or repetitive regions.
///
/// - `POOR_MAPPING_QUALITY`:
///   Represents a state where the entity has poor mapping quality, suggesting unreliable alignments or ambiguity in mapping.
///
/// # Derive Attributes
///
/// - `Debug`: Allows formatting the enum for debugging purposes.
/// - `Clone`: Enables the creation of copies of enum instances.
/// - `Copy`: Indicates that the enum can be copied without consuming the original instance.
/// - `PartialEq`: Enables comparison for equality between instances.
/// - `Eq`: Confirms the capability for strict equality comparisons.
#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalledState {
    REF_N,
    CALLABLE,
    NO_COVERAGE,
    LOW_COVERAGE,
    EXCESSIVE_COVERAGE,
    POOR_MAPPING_QUALITY,
}

/// A structure representing a range of code coverage, including its start and end
/// positions, as well as the state indicating whether the range was covered.
///
/// # Fields
/// * `start` - The starting point (inclusive) of the coverage range.
/// * `end` - The ending point (inclusive) of the coverage range.
/// * `state` - The state that indicates whether the range was covered, represented
///             by the `CalledState` enum.
/// ```
#[derive(Debug)]
pub(crate) struct CoverageRange {
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) state: CalledState,
}

/// Represents a BAM analysis report containing metadata, statistics, and detailed analysis per contig.
///
/// This struct is part of the BAM (Binary Alignment/Map) file analysis pipeline and is used to
/// encapsulate the results of analyzing a BAM file. It contains metadata associated with the
/// analysis, high-level statistics, and a list of detailed analyses for each contig.
///
/// Attributes:
/// - `metadata` (`BamMetadata`): Metadata information about the analysis, such as the file name,
///    reference genome, or date of analysis.
/// - `bam_stats` (`BamStatistics`): Aggregated statistics and metrics computed from the analyzed BAM file.
/// - `contig_analyses` (`Vec<ContigAnalysis>`): A collection of analyses that provide detailed
///    statistics and information for each individual contig in the BAM file.
///
/// This struct derives `Debug`, `Serialize`, and `Deserialize` for easier debugging and
/// compatibility with serialization frameworks like JSON or YAML.
#[derive(Debug, Serialize, Deserialize)]
pub struct BamAnalysisReport {
    pub metadata: BamMetadata,
    pub bam_stats: BamStatistics,
    pub contig_analyses: Vec<ContigAnalysis>,
}

/// Represents metadata information about a BAM file.
///
/// This struct contains details about the reference genome build, the aligner
/// used for generating the BAM file, and the number of samples included.
///
/// # Fields
/// - `reference_build` (`String`): The name of the reference genome build used
///   for aligning the sequences (e.g., "GRCh38", "hg19").
/// - `aligner` (`String`): The name of the aligner software used to process the data
///   (e.g., "BWA", "Bowtie").
/// - `sample_count` (`usize`): The number of samples included in the BAM file.
///
/// # Traits
/// - `Debug`: Enables formatting and debugging the struct.
/// - `Serialize`, `Deserialize`: Allows serialization and deserialization of this struct,
///   typically used for converting to or from formats like JSON or YAML.
/// ```
#[derive(Debug, Serialize, Deserialize)]
pub struct BamMetadata {
    pub reference_build: String,
    pub aligner: String,
    pub sample_count: usize,
}

/// A structure to represent statistics derived from BAM (Binary Alignment/Map) file analysis.
///
/// This struct contains the following fields:
///
/// - `average_read_length` (`f64`):
///   The average length of reads found in the BAM file. This value is calculated as the mean of the lengths of all reads.
///
/// - `paired_percentage` (`f64`):
///   The percentage of reads in the BAM file that are part of a paired-end sequencing experiment.
///   This represents the proportion of reads that have a mate mapped to the same sequencing context.
///
/// - `average_insert_size` (`f64`):
///   The average insert size of paired-end reads, expressed as the mean of the distances between paired read alignments.
///
/// This struct can be serialized to or deserialized from formats such as JSON using Serde.
/// It can also be easily debugged through the derived `Debug` trait.
#[derive(Debug, Serialize, Deserialize)]
pub struct BamStatistics {
    pub average_read_length: f64,
    pub paired_percentage: f64,
    pub average_insert_size: f64,
}

/// Represents the analysis information of a genomic contig.
///
/// # Fields
/// - `name`:
///   The name or identifier of the contig.
///
/// - `length`:
///   The total length of the contig in base pairs.
///
/// - `coverage_stats`:
///   Statistics related to the coverage of the contig, such as depth or distribution.
///   This is represented by an instance of `ContigCoverageStats`.
///
/// - `quality_stats`:
///   Statistics related to the quality of the contig, such as read accuracy or error rates.
///   This is represented by an instance of `ContigQualityStats`.
///
/// - `coverage_plot`:
///   An optional field containing the file path or URL to a visual representation of
///   the contig's coverage plot. If no plot is available, this field will be `None`.
///
/// # Traits
/// - `Serialize` and `Deserialize`:
///   Supports serialization and deserialization for ease of data storage or transmission.
///
/// - `Debug`:
///   Enables debugging output for the struct when needed.
#[derive(Debug, Serialize, Deserialize)]
pub struct ContigAnalysis {
    pub name: String,
    pub length: usize,
    pub coverage_stats: ContigCoverageStats,
    pub quality_stats: ContigQualityStats,
    pub coverage_plot: Option<String>,
}

/// Represents coverage statistics for a genomic contig.
///
/// This struct is used to store detailed data about the coverage of a
/// specific contig in a genomic analysis. The coverage statistics include
/// the number of unique reads, state counts, percentage coverage, average
/// depth, covered bases, and total bases.
///
/// # Fields
/// - `unique_reads` (`u64`): The number of unique reads mapped to the contig.
/// - `state_counts` (`ContigStateCounts`): A collection of counts representing
///   different states for the contig (e.g., mapped, unmapped, or error states).
/// - `coverage_percent` (`f64`): The percentage of bases in the contig that
///   are covered by at least one read.
/// - `average_depth` (`f64`): The average sequencing depth of the contig, calculated
///   as the total number of bases aligned to the contig divided by its length.
/// - `covered_bases` (`u64`): The number of bases in the contig that are covered
///   by one or more reads.
/// - `total_bases` (`u64`): The total number of bases in the contig.
///
/// # Trait Implementations
/// The struct derives the following traits:
/// - `Debug`: Enables formatting the structure for debugging purposes.
/// - `Serialize`: Allows the structure to be serialized into formats such as JSON.
/// - `Deserialize`: Allows the structure to be deserialized from formats such as JSON.
#[derive(Debug, Serialize, Deserialize)]
pub struct ContigCoverageStats {
    pub unique_reads: u64,
    pub state_counts: ContigStateCounts,
    pub coverage_percent: f64,
    pub average_depth: f64,
    pub covered_bases: u64,
    pub total_bases: u64,
}

/// A structure that represents quality statistics for a contig in genomic data analysis.
///
/// This struct provides information about mapping quality, base quality, and the percentage of bases
/// with a quality score of at least 30 (Q30), which are commonly used quality metrics in sequencing analysis.
///
/// # Fields
///
/// * `average_mapq` (`f64`):
///   The average mapping quality (MAPQ) score of the reads aligned to a contig.
///   Mapping quality represents the confidence of a read's alignment to the reference genome.
///
/// * `average_baseq` (`f64`):
///   The average base quality score of the reads aligned to a contig.
///   Base quality reflects the probability of an error in the base call at each position in the read.
///
/// * `q30_percentage` (`f64`):
///   The percentage of bases in the contig with a Phred quality score of 30 or higher (Q30).
///   A Q30 score indicates a 99.9% accuracy of the base call and is a key metric in sequencing data quality assessment.
///
/// # Derive Attributes
///
/// The `ContigQualityStats` structure derives the following traits:
///
/// * `Debug`: Allows for formatting the struct using the `{:?}` formatter, useful for debugging purposes.
/// * `Serialize` and `Deserialize`: Enables serialization and deserialization of the struct,
///   allowing it to be converted to and from data formats such as JSON, YAML, or binary formats.
///
/// This example demonstrates creating an instance of `ContigQualityStats` and printing it in debug format.
#[derive(Debug, Serialize, Deserialize)]
pub struct ContigQualityStats {
    pub average_mapq: f64,
    pub average_baseq: f64,
    pub q30_percentage: f64,
}

/// Represents a collection of state counts for a contig (a contiguous sequence of DNA).
///
/// This structure is typically used to record various statistics related to
/// the sequencing and analysis of genomic data for a specified contig.
///
/// # Fields
///
/// * `reference_n` (`u64`):
///   The count of reference 'N' nucleotides in the contig. 'N' is used
///   to represent unknown or ambiguous bases in genomic sequences.
///
/// * `no_coverage` (`u64`):
///   The count of positions in the contig where there is no sequencing coverage.
///
/// * `low_coverage` (`u64`):
///   The count of positions in the contig classified as having low sequencing coverage.
///
/// * `excessive_coverage` (`u64`):
///   The count of positions in the contig classified as having excessively high sequencing coverage,
///   potentially indicating systematic errors or repetitive sequences.
///
/// * `poor_mapping_quality` (`u64`):
///   The count of positions in the contig characterized by poor mapping quality, suggesting
///   unreliable read mapping to the reference genome.
///
/// * `callable` (`u64`):
///   The count of positions in the contig that are deemed callable, meaning they have sufficient
///   coverage and quality metrics to support reliable variant calling or confidence in base determination.
///
/// # Traits
///
/// This struct derives the following traits:
/// - `Debug`: Enables formatting of the structure for debugging purposes.
/// - `Serialize`: Allows the structure to be serialized for use with formats like JSON or other data representations.
/// - `Deserialize`: Allows the structure to be deserialized from representations like JSON, enabling data loading.
///
/// # Usage
///
/// This struct is ideal for tracking various statistics related to genomic sequencing and
/// can be a component in broader bioinformatics analyses.
#[derive(Debug, Serialize, Deserialize)]
pub struct ContigStateCounts {
    pub reference_n: u64,
    pub no_coverage: u64,
    pub low_coverage: u64,
    pub excessive_coverage: u64,
    pub poor_mapping_quality: u64,
    pub callable: u64,
}
