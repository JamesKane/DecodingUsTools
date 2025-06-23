# DecodingUs Tools

A suite of future tools intended to be used with the main [DecodingUs](https://github.com/JamesKane/decodingus) project.

Disclaimer: This project is barely ALPHA release quality.  Use at your own risk!

## Installation
```shell
git clone https://github.com/JamesKane/DecodingUsTools.git
cd DecodingUsTools
cargo build --release
```

After building, the binary will be in target/release/decodingus-tools.  Add this directory to your PATH or 
copy the binary to a location already in your PATH.

## Usage

### Coverage Analysis

```shell
decodingus-tools coverage \
  --reference <REFERENCE_FILE> \
  -o <callable_regions.bed> \
  -s <summary.html> \
  [-L  ]... \
  <BAM_FILE> 
```

Results in a callable region bed report and an interactive HTML summary for each contig in the BAM. Requires the original reference.
Intended to replace `samtools coverage` and `gatk CallableLoci` with a single tool.

Options:
- `-L, --contig <CONTIG>`: Limit analysis to specific contigs. Can be specified multiple times. (Example: `-L chr1 -L chr2`)
- `-r, --reference <FILE>`: Reference FASTA file
- `-o, --output <FILE>`: Output BED file (default: callable_regions.bed)
- `-s, --summary <FILE>`: Output HTML summary (default: summary.html)


Example output: (callable_regions.bed)
```text
chr1	10000	10006	LOW_COVERAGE
chr1	10007	10131	POOR_MAPPING_QUALITY
chr1	10132	10170	LOW_COVERAGE
chr1	10171	10228	POOR_MAPPING_QUALITY
chr1	10229	10402	LOW_COVERAGE
chr1	10403	10405	POOR_MAPPING_QUALITY
chr1	10406	10406	LOW_COVERAGE
chr1	10407	10407	CALLABLE
chr1	10408	10410	LOW_COVERAGE
chr1	10411	10412	CALLABLE
```

The HTML summary provides an interactive report with:
- Overall BAM statistics (read length, paired percentage, insert size)
- Per-contig statistics including:
    - Coverage metrics (callable bases, low coverage regions, etc.)
    - Average depth and quality scores
    - Coverage distribution plots
- Interactive contig selection via searchable dropdown
- SVG coverage plots for visual analysis

The summary can be viewed in any modern web browser.

### Y-DNA Branch Finding

```shell
decodingus-tools find-y-branch \
  --reference <REFERENCE_FILE> \
  --provider <PROVIDER> \
[--show-snps] \
  <BAM_FILE> <OUTPUT_FILE>
```

Find the closest YDNA branch. Available providers:
- `ftdna` (default): Uses FamilyTreeDNA's haplotree
- `decodingus`: Uses DecodingUs project haplotree

Example output:

|Haplogroup|Score|Matching_SNPs|Mismatching_SNPs|Ancestral_Matches|No_Calls|Total_SNPs|Cumulative_SNPs|Depth|
|----------|------|-------------|----------------|-----------------|--------|----------|---------------|-----|
|R-FGC29071|2.75|3|0|2|1|9|1917|53|
|R-FGC29067|3.08|7|0|0|2|4|1911|52|

### mtDNA Branch Finding

```shell
decodingus-tools find-mt-branch \
--reference <REFERENCE_FILE> \
--provider <PROVIDER> \
[--show-snps] \
 <BAM_FILE> <OUTPUT_FILE>
```

Find the closest mtDNA branch. Available providers:
- `ftdna` (default): Uses FamilyTreeDNA's haplotree

Example output:

|Haplogroup| Score |Matching_SNPs|Mismatching_SNPs|Ancestral_Matches|No_Calls|Total_SNPs|Cumulative_SNPs|Depth|
|----------|-------|-------------|----------------|-----------------|--------|----------|---------------|-----|
|U5a1b1g|1.65|1|0|0|0|1|55|15|

### Fix Surjected BAM

```shell
decodingus-tools fix-surjected-bam \
  --reference <REFERNECE_FILE> \ 
  -o <OUTPUT_BAM> \ 
  <SURJECTED_BAM>
```

When a GAM file is surjected back via `vg surject` to a linear BAM, the results are mixed up compared to a traditional linear reference. This automates:
- Reheadering the @SQ details with a known reference to match the order
- Removing the PanSN-spec prefix from the individual reads
- Invoking samtools to sort the final result

### Generate Sequence Fingerprint

```shell
decodingus-tools fingerprint \
[-r <REFERENCE_FILE>] \
[-o <OUTPUT_FILE>] \
[-R  ] \
[--ksize <K-MER_SIZE>] \
[--scaled <SCALED_FACTOR>] \
<INPUT_FILE>
```
#### The tool supports various input formats:
- BAM/CRAM files (requires reference for CRAM)
- FASTQ files (.fastq, .fq, optionally gzipped)
- GAM files (Graph Alignment/Map format from vg)

### Optimized MinHash Variant for Sequence Fingerprinting
This fingerprinting method utilizes a performance-enhanced variant of the MinHash algorithm, designed to efficiently identify the original source library of sequencing reads processed through downstream pipelines. While effective for certain applications, its use is context-dependent and is best suited for specific scenarios.
#### Key Use Cases:
- **Source Library Identification**: Ideal for determining whether raw sequencing reads originated from the same source library, even if they were processed using different alignment tools or references.
- **Pipeline Validation**: Enables cross-comparisons of sequencing data processed via diverse downstream bioinformatics workflows, providing insights into pipeline consistency.

#### Limitations:
- **Impact of Wet Lab and Sequencer Variability**: When sequencing reads are generated using different instrumentation (e.g., sequencing platforms) or wet lab protocols (e.g., library prep kits), Jaccard similarity scores between datasets from the same donor can appear as distant as those from completely unrelated individuals.
- **Population and Family Genetics**: This method is not designed to reliably identify individual donors, familial relationships, or relatedness in population genomics for human-sized genomes, due to the loss of subtle genetic variations in the fingerprinting process.

#### Jaccard Similarity Context:
- **Effective Use Cases**: Jaccard similarity scores are meaningful when comparing fingerprints derived from the same raw sequencing datasets, such as raw FASTQ files or aligned outputs from varying pipelines (e.g., different alignment tools or genome references).
- **Observed Failure Scenarios**: Similarity calculations fail when datasets are derived using differing sequencing technologies or library preparation methods, with results resembling comparisons between unrelated samples rather than same-origin datasets. This behavior is inherent to MinHash-based algorithms and is not unique to this implementation.

Notably, a similar behavior is observed in tools like **sourmash**, which is frequently used in the literature for assessing relatedness in human samples. Therefore, this is a characteristic of the algorithm itself, rather than a limitation of this particular implementation.
#### Performance Characteristics:
This implementation strikes a balance between speed and resource efficiency:
- **Whole Genome Fingerprinting**: Processes 30x WGS datasets (~90 billion bases) in ~20–30 minutes on modern processors.
- **Region-Specific Analysis** (e.g., Y Chromosome or Mitochondrial DNA): Completes in under 2 minutes with minimal memory usage.

#### Output and Interpretation:
- **Similarity Analysis**: Produces compact k-mer files that allow rapid Jaccard similarity comparisons between datasets, provided the datasets share the same raw source or aligned set of reads. Users must interpret results cautiously when fingerprints are generated using different sequencing methods or workflows.
- **File Format**: Includes metadata (k-mer size, scaling factor, and region specificity) to ensure output reproducibility.

#### Practical Considerations:
- **Parameter Tuning**: Adjust k-mer size (`--ksize`) and scaling (`--scaled`) based on genome size and research focus. Larger datasets may benefit from fine-tuning these parameters for improved performance.
- **Raw Dataset Comparability**: This method works best when applied to fingerprints derived from the same set of source reads, even if processed differently downstream.
- **Runtime Efficiency**: While whole-genome analysis may take longer and require overnight processing for batch comparisons, region-specific analyses (e.g., specific chromosomes) typically complete within minutes.

This MinHash-based implementation is especially suited for identifying and validating sequencing library origins across bioinformatics pipelines. Its behavior in defining sample similarity mirrors that of other tools like **sourmash** and should be considered inherent to the algorithm itself rather than a tool-specific limitation. Users should employ this method with care, particularly in contexts involving diverse wet lab preparations or sequencing technologies.

### Configuration

The tool supports configuration via a TOML file located at:
* Linux: ~/.config/decodingus-tools/config.toml
* macOS: ~/Library/Application Support/com.decodingus.decodingus-tools/config.toml
* Windows: %APPDATA%\decodingus-tools\config.toml

Example configuration file:
```toml
# Time in seconds to wait for haplotree downloads
download_timeout = 120
```

Available configuration options:
* `download_timeout`: Time in seconds to wait for haplotree downloads

You can create or modify this file manually to adjust the tool's behavior. The configuration file is optional - the tool will use default values if the file doesn't exist.


### Caching
The branch finding algorithms store files in these locations:
* Linux: ~/.cache/decodingus-tools
* macOS: ~/Library/Caches/com.decodingus.decodingus-tools
* Windows: %LOCALAPPDATA%\decodingus-tools\Cache

Each tree will only be retrieved once per week if the file exists.  Remove the weekly file if you wish to force a reload.

## [![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

This project is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.

