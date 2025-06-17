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
  [--ksize <K-MER_SIZE>] \
  [--scaled <SCALED_FACTOR>] \
  <INPUT_FILE>
```
Generate a unique MinHash-based fingerprint for sequencing data files (BAM, CRAM, or FASTQ). Useful for:
- Quick sequence similarity comparison
- Data provenance tracking
- Identifying duplicate or related samples

Options:
- : Reference FASTA file (required for CRAM files) `-r, --reference <FILE>`
- `--ksize <INT>`: K-mer size for MinHash sketching (default: 31)
- `--scaled <INT>`: MinHash scaled factor (default: 1000)

The tool processes the input file in parallel using available CPU cores, providing progress updates every few seconds.


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

