use crate::callable_loci::types::ContigQualityStats;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

///
/// A struct representing the export of coverage data, including a summary,
/// contigs information, and quality metrics. This structure is designed to
/// be serializable and deserializable using Serde.
///
/// # Fields
///
/// * `summary` (`CoverageSummary`) - Provides a summary of the coverage information.
/// * `contigs` (`Arc<Vec<ContigExport>>`) - A reference-counted pointer to a vector
///   of contigs export data. Serialized and deserialized with custom Serde support
///   using the `arc_serde` module.
/// * `quality_metrics` (`QualityMetrics`) - Contains quality metrics related
///   to the coverage data.
///
/// # Attributes
///
/// * `#[derive(Debug, Serialize, Deserialize)]` - Automatically implements the `Debug`,
///   `Serialize`, and `Deserialize` traits for debugging and (de)serialization.
/// * `#[serde(with = "arc_serde")]` - Ensures the proper serialization and deserialization
///   of the `Arc<Vec<ContigExport>>` field.
/// ```
#[derive(Debug, Serialize, Deserialize)]
pub struct CoverageExport {
    pub summary: CoverageSummary,
    #[serde(with = "arc_serde")]
    pub contigs: Arc<Vec<ContigExport>>,
    pub quality_metrics: QualityMetrics,
}

// Implementation of Clone using Arc for the non-cloneable Vec<ContigExport>
impl Clone for CoverageExport {
    fn clone(&self) -> Self {
        CoverageExport {
            summary: self.summary.clone(),
            contigs: Arc::clone(&self.contigs),
            quality_metrics: self.quality_metrics.clone(),
        }
    }
}

// Helper module for serializing Arc<Vec<T>>
mod arc_serde {
    use serde::{Deserialize, Deserializer, Serialize, Serializer};
    use std::sync::Arc;

    pub fn serialize<S, T>(val: &Arc<Vec<T>>, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
        T: Serialize,
    {
        (**val).serialize(serializer)
    }

    pub fn deserialize<'de, D, T>(deserializer: D) -> Result<Arc<Vec<T>>, D::Error>
    where
        D: Deserializer<'de>,
        T: Deserialize<'de>,
    {
        Vec::deserialize(deserializer).map(Arc::new)
    }
}

/// The `CoverageSummary` struct provides a summary of sequencing coverage statistics.
/// It is designed to represent coverage-related metrics commonly used in bioinformatics analysis,
/// particularly in the context of sequencing reads and alignment data.
///
/// # Fields
///
/// - `total_bases` (`u64`):
///   The total number of sequenced bases. This includes all bases in the analyzed dataset,
///   regardless of their quality or usability.
///
/// - `callable_bases` (`u64`):
///   The number of bases considered "callable" after applying quality and coverage filters.
///   Callable bases are typically those that meet thresholds for sequencing depth and base quality,
///   making them suitable for variant calling or downstream analysis.
///
/// - `callable_percentage` (`f64`):
///   The percentage of callable bases relative to the total bases.
///   This is calculated as: `(callable_bases / total_bases) * 100.0`.
///   It represents the proportion of usable bases in the dataset.
///
/// - `average_depth` (`f64`):
///   The average sequencing depth (coverage) over the analyzed regions or dataset.
///   Depth is usually defined as the average number of sequencing reads overlapping a base position.
///
/// - `contigs_analyzed` (`usize`):
///   The number of contigs (genomic regions or sequences) analyzed during the calculation of coverage statistics.
///   This field represents the breadth of the dataset or analysis scope.
///
/// # Derive Attributes
///
/// - `Debug`:
///   Enables formatted printing of `CoverageSummary` instances for debugging purposes.
///
/// - `Serialize` and `Deserialize`:
///   Makes the struct compatible with serialization formats (e.g., JSON, YAML, etc.).
///   This is useful for saving and loading `CoverageSummary` data.
///
/// - `Clone`:
///   Allows for creating exact copies of `CoverageSummary` instances.
///
/// This struct is useful for summarizing and reporting coverage-related results in genomics pipelines.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct CoverageSummary {
    pub aligner: String,
    pub reference_build: String,
    pub sequencing_platform: String,
    pub read_length: usize,
    pub total_bases: u64,
    pub callable_bases: u64,
    pub callable_percentage: f64,
    pub average_depth: f64,
    pub contigs_analyzed: usize,
}

/// A structure representing the exportable data of a contig in a sequencing or bioinformatics analysis.
///
/// This struct contains metadata and statistics about a contig, such as its name, length, read coverage,
/// and quality metrics, as well as distribution and optional visualization data.
///
/// # Fields
/// * `name` - The name or identifier of the contig.
/// * `length` - The length of the contig in base pairs (bp).
/// * `unique_reads` - The number of reads mapped uniquely to the contig.
/// * `coverage_percent` - The percentage of the contig covered by sequencing reads.
/// * `average_depth` - The average sequencing depth (coverage) across the contig.
/// * `covered_bases` - The number of bases in the contig that are covered by sequencing reads.
/// * `total_bases` - The total number of base pairs in the contig (same as `length`).
/// * `quality_stats` - A `ContigQualityStats` struct, providing detailed quality metrics for the contig.
/// * `state_distribution` - A `StateDistribution` struct, describing the distribution of states
///   (e.g., different types of bases or genomic features) across the contig.
/// * `coverage_histogram` - An optional field containing plot data for rendering coverage histograms
///   in a graphical user interface (GUI). This is a `Vec` of `HistogramBin` and is serialized only if present.
///
/// # Notes
/// - The use of the `#[serde(skip_serializing_if = "Option::is_none")]` attribute ensures that the
///   `coverage_histogram` field is omitted from the serialized output if it is `None`.
/// - The structure is derived with `Debug`, `Serialize`, and `Deserialize` traits for easier debugging and
///   seamless integration with serialization/deserialization libraries like `serde`.
#[derive(Debug, Serialize, Deserialize)]
pub struct ContigExport {
    pub name: String,
    pub length: usize,
    pub unique_reads: u64,
    pub coverage_percent: f64,
    pub average_depth: f64,
    pub covered_bases: u64,
    pub total_bases: u64,

    pub quality_stats: ContigQualityStats,
    pub state_distribution: StateDistribution,

    // Optional: include plot data for GUI rendering
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coverage_histogram: Option<Vec<HistogramBin>>,
}

/// Represents the distribution of various states in a given context with quantitative measures.
///
/// This structure is serializable and deserializable, supports debugging,
/// and can be cloned for enhanced usability. Each field corresponds to a specific
/// state category with `u64` values representing their respective counts.
///
/// # Fields
///
/// * `ref_n` - The count of references categorized under "N", typically used to represent undefined or unclassified bases.
/// * `callable` - The count of callable states, indicating the regions that meet the criteria for being callable.
/// * `no_coverage` - The count of regions or states with no observed coverage.
/// * `low_coverage` - The count of regions or states with lower-than-threshold coverage.
/// * `excessive_coverage` - The count of regions or states with coverage far exceeding the expected threshold.
/// * `poor_mapping_quality` - The count of regions or states with subpar mapping quality, signaling unreliable data or mapping mismatches.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct StateDistribution {
    pub ref_n: u64,
    pub callable: u64,
    pub no_coverage: u64,
    pub low_coverage: u64,
    pub excessive_coverage: u64,
    pub poor_mapping_quality: u64,
}

/// Represents a single bin in a histogram, defined by a range and an associated state.
///
/// This struct is typically used for grouping data points within a specified range (`start` to `end`)
/// and associating that range with a descriptive state or label.
///
/// # Fields
/// - `start` (u32): The starting value of the range for this histogram bin (inclusive).
/// - `end` (u32): The ending value of the range for this histogram bin (exclusive).
/// - `state` (String): A string representing the state or label associated with this bin,
///   providing additional context or categorization.
///
/// # Traits
/// This struct derives the following traits:
/// - `Debug`: Enables formatting with `{:?}` for debugging purposes.
/// - `Serialize`: Allows serialization of `HistogramBin` instances, which is useful for encoding the data in formats such as JSON.
/// - `Deserialize`: Allows deserialization of data into `HistogramBin` instances, enabling decoding from formats like JSON.
/// - `Clone`: Allows creating a duplicate of a `HistogramBin` instance.
///
/// ```
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct HistogramBin {
    pub start: u32,
    pub end: u32,
    pub state: String,
}

/// The `QualityMetrics` structure represents a set of quality control metrics
/// for sequencing or genomic data. This structure is serializable and
/// deserializable utilizing Serde and can be cloned or debugged.
///
/// # Fields
///
/// * `average_mapq` (`f64`):
///   The average mapping quality of the sequencing reads. Mapping quality
///   is typically used to assess the alignment confidence of the reads
///   to the reference genome.
///
/// * `average_baseq` (`f64`):
///   The average base quality score across all bases in the reads.
///   Base quality scores indicate the probability of a base call being
///   incorrect, with higher values representing better quality.
///
/// * `q30_percentage` (`f64`):
///   The percentage of bases with a Phred base quality score of 30 or higher.
///   A Phred score of 30 corresponds to a 1 in 1,000 chance of an incorrect
///   base call, which is often used as a standard high-quality threshold in
///   sequencing data.
///
/// # Traits Implemented
///
/// * `Debug`: Enables formatting the structure using the `{:?}` formatter for debugging purposes.
/// * `Serialize` and `Deserialize`: Allows for easy serialization and deserialization
///   of the structure to/from formats like JSON or TOML.
/// * `Clone`: Allows the structure to be duplicated.
/// ```
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct QualityMetrics {
    pub average_mapq: f64,
    pub average_baseq: f64,
    pub q30_percentage: f64,
}