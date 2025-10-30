use crate::callable_loci::detect_aligner;
use crate::types::ReferenceGenome;
use crate::utils::bam_reader::BamReaderFactory;
use crate::utils::progress_manager::ProgressManager;
use rust_htslib::bam::Read;
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SequencingPlatform {
    Illumina,
    PacBio,
    Nanopore,
    MGI,
    Unknown,
}

pub struct BamStats {
    read_count: usize,
    total_read_length: usize,
    paired_reads: usize,
    total_insert_size: i64,
    paired_count: usize,
    pub(crate) max_samples: usize,
    length_distribution: HashMap<usize, usize>,
    insert_size_distribution: HashMap<i64, usize>,
    aligner: String,
    reference_build: String,
    flow_cells: HashMap<String, usize>,
    instruments: HashMap<String, usize>,
    platform_counts: HashMap<SequencingPlatform, usize>,
}

impl BamStats {
    pub fn new(max_samples: usize) -> Self {
        BamStats {
            read_count: 0,
            total_read_length: 0,
            paired_reads: 0,
            total_insert_size: 0,
            paired_count: 0,
            max_samples,
            length_distribution: HashMap::new(),
            insert_size_distribution: HashMap::new(),
            aligner: String::new(),
            reference_build: String::new(),
            flow_cells: HashMap::new(),
            instruments: HashMap::new(),
            platform_counts: HashMap::new(),
        }
    }

    pub fn collect_stats(
        &mut self,
        bam_path: &str,
        reference_path: Option<&str>,
        progress_mgr: &ProgressManager,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut bam = BamReaderFactory::open(bam_path, reference_path)?;
        let header = bam.header().clone();

        self.aligner = detect_aligner(&header);
        self.reference_build = ReferenceGenome::from_header(&header)
            .map(|g| g.name().to_string())
            .unwrap_or_else(|| "Unknown".to_string());

        let progress = progress_mgr.add_spinner("Collecting BAM statistics...");

        for (i, record_result) in bam.records().enumerate() {
            if i >= self.max_samples {
                break;
            }

            let record = record_result?;

            // Only count primary alignments
            if !record.is_secondary() && !record.is_supplementary() {
                let seq = record.seq();
                let seq_len = seq.len();

                *self.length_distribution.entry(seq_len).or_insert(0) += 1;

                self.read_count += 1;
                self.total_read_length += seq_len;

                // Extract flow cell and instrument from read name
                if let Ok(qname) = std::str::from_utf8(record.qname()) {
                    let platform = Self::detect_platform_from_qname(qname);
                    *self.platform_counts.entry(platform.clone()).or_insert(0) += 1;

                    match platform {
                        SequencingPlatform::Illumina => {
                            if let Some((instrument, flow_cell)) =
                                Self::parse_illumina_read_name(qname)
                            {
                                *self.instruments.entry(instrument.to_string()).or_insert(0) += 1;
                                *self.flow_cells.entry(flow_cell.to_string()).or_insert(0) += 1;
                            }
                        }
                        SequencingPlatform::PacBio => {
                            if let Some(instrument) = Self::parse_pacbio_read_name(qname) {
                                *self.instruments.entry(instrument.to_string()).or_insert(0) += 1;
                            }
                        }
                        SequencingPlatform::Nanopore => {
                            if let Some(instrument) = Self::parse_nanopore_read_name(qname) {
                                *self.instruments.entry(instrument.to_string()).or_insert(0) += 1;
                            }
                        }
                        SequencingPlatform::MGI => {
                            if let Some((instrument, flow_cell)) = Self::parse_mgi_read_name(qname)
                            {
                                *self.instruments.entry(instrument.to_string()).or_insert(0) += 1;
                                *self.flow_cells.entry(flow_cell.to_string()).or_insert(0) += 1;
                            }
                        }
                        SequencingPlatform::Unknown => {}
                    }
                }

                if record.is_paired() {
                    self.paired_reads += 1;

                    if record.is_proper_pair() && record.is_first_in_template() {
                        let insert_size = record.insert_size().abs();
                        if insert_size > 0 {
                            *self
                                .insert_size_distribution
                                .entry(insert_size)
                                .or_insert(0) += 1;
                            self.total_insert_size += insert_size;
                            self.paired_count += 1;
                        }
                    }
                }
            }

            if i % 10000 == 0 {
                progress.set_message(format!("Processed {} reads...", i));
            }
        }

        progress.finish_with_message("BAM Statistics collected");
        Ok(())
    }

    /// Detect sequencing platform from read name format
    fn detect_platform_from_qname(qname: &str) -> SequencingPlatform {
        // Oxford Nanopore format: typically UUID-like or has "runid" in the name
        // Example: 0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4c5d or read_id with channel info
        // Nanopore reads often have very long alphanumeric IDs with dashes or underscores
        if qname.len() > 30 && (qname.contains('-') || qname.contains('_')) {
            // Check for UUID pattern (8-4-4-4-12 hex digits)
            let parts: Vec<&str> = qname.split('-').collect();
            if parts.len() == 5 {
                let is_uuid = parts[0].len() == 8
                    && parts[1].len() == 4
                    && parts[2].len() == 4
                    && parts[3].len() == 4
                    && parts[4].len() >= 12;
                if is_uuid
                    && parts
                        .iter()
                        .all(|p| p.chars().all(|c| c.is_ascii_hexdigit()))
                {
                    return SequencingPlatform::Nanopore;
                }
            }
            // Alternative Nanopore format with channel/read info
            if qname.contains("ch") && qname.contains("read") {
                return SequencingPlatform::Nanopore;
            }
        }

        // PacBio format: m<instrument>_<date>_<time>/<zmw>/<subread_type>
        // Example: m64023e_230414_133043/1/ccs
        if qname.starts_with('m') && qname.contains('/') {
            let parts: Vec<&str> = qname.split('/').collect();
            if parts.len() >= 2 && parts[0].contains('_') {
                return SequencingPlatform::PacBio;
            }
        }

        // MGI/BGI format: Similar to Illumina but with specific prefixes
        // Example: V300012345L1C001R00100000001 or CL100012345L1C001R001_1
        // MGI instruments start with specific prefixes: V300, E100, CL100, G99, G400, etc.
        if qname.len() > 15 {
            let prefix = &qname[..5].to_uppercase();
            if prefix.starts_with("V300")
                || prefix.starts_with("E100")
                || prefix.starts_with("CL100")
                || prefix.starts_with("G400")
                || prefix.starts_with("G99")
            {
                return SequencingPlatform::MGI;
            }
            // Another MGI format with colons but different structure than Illumina
            if qname.matches(':').count() >= 6 {
                // Check if it looks like MGI by examining the instrument field
                let parts: Vec<&str> = qname.split(':').collect();
                if parts[0].starts_with('V')
                    || parts[0].starts_with('E')
                    || parts[0].starts_with("CL")
                    || parts[0].starts_with('G')
                {
                    // Further check: MGI flow cell IDs often have 'L' followed by digits
                    if parts.len() >= 3 && parts[2].starts_with('L') {
                        return SequencingPlatform::MGI;
                    }
                }
            }
        }

        // Illumina format: has multiple colons
        // Example: A00123:123:HXXXYDRXX:1:1101:1000:1000
        if qname.matches(':').count() >= 6 {
            return SequencingPlatform::Illumina;
        }

        SequencingPlatform::Unknown
    }

    /// Parse Illumina-style read names to extract instrument ID and flow cell ID
    /// Format: @<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y> <read>:<filtered>:<control>:<index>
    /// Example: @A00123:123:HXXXYDRXX:1:1101:1000:1000 1:N:0:ATCACG
    fn parse_illumina_read_name(qname: &str) -> Option<(&str, &str)> {
        let parts: Vec<&str> = qname.split(':').collect();
        if parts.len() >= 3 {
            let instrument = parts[0];
            let flow_cell = parts[2];
            Some((instrument, flow_cell))
        } else {
            None
        }
    }

    /// Parse PacBio read names to extract instrument ID
    /// Format: m<instrument>_<date>_<time>/<zmw>/<subread_type>
    /// Example: m64023e_230414_133043/1/ccs
    /// Returns the instrument ID (e.g., m64023e)
    fn parse_pacbio_read_name(qname: &str) -> Option<&str> {
        if let Some(slash_pos) = qname.find('/') {
            let movie_name = &qname[..slash_pos];
            if movie_name.starts_with('m') {
                // Extract just the instrument ID (e.g., m64023e)
                if let Some(first_underscore) = movie_name.find('_') {
                    return Some(&movie_name[..first_underscore]);
                }
            }
        }
        None
    }

    /// Parse Oxford Nanopore read names to extract device/run ID
    /// Format varies: UUID-based or descriptive format
    /// Example: 0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4c5d or similar
    /// Returns a simplified identifier
    fn parse_nanopore_read_name(qname: &str) -> Option<&str> {
        // For Nanopore, we'll use a portion of the read name as the "instrument"
        // This could be the run ID or device ID depending on the naming scheme

        // If it's a UUID format, take the first segment
        if qname.len() > 30 && qname.contains('-') {
            let parts: Vec<&str> = qname.split('-').collect();
            if parts.len() >= 5 {
                // Return a combination of first segments to represent the run/device
                return Some(
                    qname
                        .split('_')
                        .next()
                        .unwrap_or(qname)
                        .split('-')
                        .next()
                        .unwrap_or("nanopore"),
                );
            }
        }

        // For other formats, try to extract a meaningful prefix
        if let Some(underscore_pos) = qname.find('_') {
            return Some(&qname[..underscore_pos]);
        }

        Some("nanopore")
    }

    /// Parse MGI-Seq read names to extract instrument ID and flow cell ID
    /// Format 1: V300012345L1C001R00100000001 (no colons)
    /// Format 2: V300012345:L1:C001:R001 (with colons, similar to Illumina)
    /// Returns (instrument_id, flow_cell_id)
    fn parse_mgi_read_name(qname: &str) -> Option<(&str, &str)> {
        // Format with colons (similar to Illumina)
        if qname.matches(':').count() >= 3 {
            let parts: Vec<&str> = qname.split(':').collect();
            if parts.len() >= 3 {
                let instrument = parts[0];
                let flow_cell = parts[1]; // Lane is often the flow cell identifier in MGI
                return Some((instrument, flow_cell));
            }
        }

        // Format without colons: V300012345L1C001R00100000001
        // Try to extract instrument prefix and lane info
        if qname.len() > 10 {
            if let Some(l_pos) = qname.find('L') {
                let instrument = &qname[..l_pos];
                // Extract lane+chip as flow cell identifier
                if let Some(c_pos) = qname[l_pos..].find('C') {
                    let end_pos = qname[l_pos..].find('R').unwrap_or(qname[l_pos..].len());
                    let flow_cell = &qname[l_pos..l_pos + end_pos];
                    return Some((instrument, flow_cell));
                }
            }
        }

        None
    }

    pub fn get_stats(&self) -> HashMap<String, f64> {
        let mut stats = HashMap::new();

        if self.read_count > 0 {
            let modal_length = self
                .length_distribution
                .iter()
                .max_by_key(|(_, &count)| count)
                .map(|(length, _)| *length)
                .unwrap_or(0);

            stats.insert("average_read_length".to_string(), modal_length as f64);
            stats.insert(
                "paired_percentage".to_string(),
                (self.paired_reads as f64 / self.read_count as f64) * 100.0,
            );
        }

        if self.paired_count > 0 {
            let modal_insert_size = self
                .insert_size_distribution
                .iter()
                .max_by_key(|(_, &count)| count)
                .map(|(size, _)| *size)
                .unwrap_or(0);

            stats.insert("average_insert_size".to_string(), modal_insert_size as f64);
            stats.insert(
                "proper_pair_percentage".to_string(),
                (self.paired_count as f64 * 2.0 / self.paired_reads as f64) * 100.0,
            );
        }

        stats
    }

    pub fn aligner(&self) -> &str {
        &self.aligner
    }

    pub fn reference_build(&self) -> &str {
        &self.reference_build
    }

    pub fn average_read_length(&self) -> usize {
        if self.read_count > 0 {
            self.total_read_length / self.read_count
        } else {
            0
        }
    }

    pub fn modal_read_length(&self) -> usize {
        if self.read_count > 0 {
            self.length_distribution
                .iter()
                .max_by_key(|(_, &count)| count)
                .map(|(length, _)| *length)
                .unwrap_or(0)
        } else {
            0
        }
    }

    pub fn flow_cells(&self) -> &HashMap<String, usize> {
        &self.flow_cells
    }

    pub fn instruments(&self) -> &HashMap<String, usize> {
        &self.instruments
    }

    pub fn get_primary_platform(&self) -> SequencingPlatform {
        self.platform_counts
            .iter()
            .max_by_key(|(_, &count)| count)
            .map(|(platform, _)| platform.clone())
            .unwrap_or(SequencingPlatform::Unknown)
    }

    pub fn infer_platform(&self) -> String {
        let primary_platform = self.get_primary_platform();

        match primary_platform {
            SequencingPlatform::PacBio => {
                // PacBio instrument prefixes:
                // m64xxx = Sequel II/IIe
                // m54xxx = Sequel I
                // m84xxx = Revio
                if let Some((instrument_id, _)) =
                    self.instruments.iter().max_by_key(|(_, &count)| count)
                {
                    if instrument_id.starts_with("m84") {
                        return "PacBio Revio".to_string();
                    } else if instrument_id.starts_with("m64") {
                        return "PacBio Sequel II/IIe".to_string();
                    } else if instrument_id.starts_with("m54") {
                        return "PacBio Sequel".to_string();
                    } else {
                        return "PacBio".to_string();
                    }
                }
                "PacBio".to_string()
            }
            SequencingPlatform::Nanopore => {
                // Oxford Nanopore platforms
                // Could infer MinION, GridION, PromethION based on read characteristics
                // For now, return generic Nanopore
                "Oxford Nanopore".to_string()
            }
            SequencingPlatform::MGI => {
                // MGI/BGI instrument prefixes:
                // V300 = DNBseq/MGISEQ-2000
                // E100 = MGISEQ-200
                // CL100 = MGISEQ-T7
                // G400 = DNBSEQ-G400
                // G99 = MGISEQ-T1
                if let Some((instrument_id, _)) =
                    self.instruments.iter().max_by_key(|(_, &count)| count)
                {
                    if instrument_id.starts_with("V300") {
                        return "MGI DNBSEQ/MGISEQ-2000".to_string();
                    } else if instrument_id.starts_with("E100") {
                        return "MGI MGISEQ-200".to_string();
                    } else if instrument_id.starts_with("CL100") {
                        return "MGI MGISEQ-T7".to_string();
                    } else if instrument_id.starts_with("G400") {
                        return "MGI DNBSEQ-G400".to_string();
                    } else if instrument_id.starts_with("G99") {
                        return "MGI MGISEQ-T1".to_string();
                    } else {
                        return "MGI DNBseq".to_string();
                    }
                }
                "MGI DNBseq".to_string()
            }
            SequencingPlatform::Illumina => {
                // Most common instrument ID
                if let Some((instrument_id, _)) =
                    self.instruments.iter().max_by_key(|(_, &count)| count)
                {
                    // Infer platform based on instrument ID prefix
                    // See: https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm
                    let prefix = instrument_id.chars().next().unwrap_or(' ');
                    match prefix {
                        'A' | 'a' => "NovaSeq".to_string(),
                        'D' | 'd' => "HiSeq 2500".to_string(),
                        'J' | 'j' => "HiSeq 3000".to_string(),
                        'K' | 'k' => "HiSeq 4000".to_string(),
                        'E' | 'e' => "HiSeq X".to_string(),
                        'N' | 'n' => "NextSeq".to_string(),
                        'M' | 'm' => "MiSeq".to_string(),
                        'V' | 'v' => "NovaSeq X".to_string(),
                        'F' | 'f' => "iSeq".to_string(),
                        _ => "Unknown Illumina".to_string(),
                    }
                } else {
                    "Unknown Illumina".to_string()
                }
            }
            SequencingPlatform::Unknown => "Unknown".to_string(),
        }
    }

    pub fn get_instrument_summary(&self) -> Vec<(String, usize, f64)> {
        let total_reads = self.read_count as f64;
        let mut summary: Vec<_> = self
            .instruments
            .iter()
            .map(|(id, &count)| {
                let percentage = (count as f64 / total_reads) * 100.0;
                (id.clone(), count, percentage)
            })
            .collect();

        summary.sort_by(|a, b| b.1.cmp(&a.1));
        summary
    }

    pub fn get_flow_cell_summary(&self) -> Vec<(String, usize, f64)> {
        let total_reads = self.read_count as f64;
        let mut summary: Vec<_> = self
            .flow_cells
            .iter()
            .map(|(id, &count)| {
                let percentage = (count as f64 / total_reads) * 100.0;
                (id.clone(), count, percentage)
            })
            .collect();

        summary.sort_by(|a, b| b.1.cmp(&a.1));
        summary
    }

    pub fn get_platform_summary(&self) -> Vec<(String, usize, f64)> {
        let total_reads = self.read_count as f64;
        let mut summary: Vec<_> = self
            .platform_counts
            .iter()
            .map(|(platform, &count)| {
                let platform_name = match platform {
                    SequencingPlatform::Illumina => "Illumina",
                    SequencingPlatform::PacBio => "PacBio",
                    SequencingPlatform::Nanopore => "Oxford Nanopore",
                    SequencingPlatform::MGI => "MGI DNBseq",
                    SequencingPlatform::Unknown => "Unknown",
                }
                .to_string();
                let percentage = (count as f64 / total_reads) * 100.0;
                (platform_name, count, percentage)
            })
            .collect();

        summary.sort_by(|a, b| b.1.cmp(&a.1));
        summary
    }
}
