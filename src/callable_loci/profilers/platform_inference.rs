use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SequencingPlatform {
    Illumina,
    PacBio,
    Nanopore,
    MGI,
    Unknown,
}

/// Helper struct for inferring sequencing platforms from BAM file metadata
pub struct PlatformInference;

impl PlatformInference {
    /// Detect sequencing platform from read name format
    pub fn detect_platform_from_qname(qname: &str) -> SequencingPlatform {
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
    pub fn parse_illumina_read_name(qname: &str) -> Option<(&str, &str)> {
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
    pub fn parse_pacbio_read_name(qname: &str) -> Option<&str> {
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
    pub fn parse_nanopore_read_name(qname: &str) -> Option<&str> {
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
    pub fn parse_mgi_read_name(qname: &str) -> Option<(&str, &str)> {
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

    /// Infers the sequencing platform based on the primary platform and the associated instrument IDs.
    ///
    /// This method uses specific rules and conventions to deduce the platform name from the available data.
    /// It categorizes platforms into PacBio, Oxford Nanopore, MGI, Illumina, or Unknown, with further distinctions
    /// for PacBio, MGI, and Illumina based on instrument-specific prefixes.
    ///
    /// # Returns
    /// A `String` representing the inferred sequencing platform. Possible outputs include:
    /// - `"PacBio Revio"`, `"PacBio Sequel II/IIe"`, `"PacBio Sequel"`, or `"PacBio"` (for PacBio platforms)
    /// - `"Oxford Nanopore"` (for Nanopore platforms)
    /// - `"MGI DNBSEQ/MGISEQ-2000"`, `"MGI MGISEQ-200"`, `"MGI MGISEQ-T7"`, `"MGI DNBSEQ-G400"`, `"MGI MGISEQ-T1"`, or `"MGI DNBseq"` (for MGI platforms)
    /// - Specific model names such as `"NovaSeq"`, `"HiSeq 2500"`, `"HiSeq 3000"`, `"HiSeq 4000"`, `"HiSeq X"`, `"NextSeq"`, `"MiSeq"`, `"NovaSeq X"`, `"iSeq"`
    ///   or `"Unknown Illumina"` (for Illumina platforms)
    /// - `"Unknown"` if the platform can't be determined.
    ///
    /// # Method
    /// - For **PacBio** platforms, uses instrument prefixes like `m84xxx` for Revio, `m64xxx` for Sequel II/IIe, and `m54xxx` for Sequel I.
    /// - For **Oxford Nanopore**, returns a generic "Oxford Nanopore" since detailed inference isn't implemented.
    /// - For **MGI**, uses prefixes such as `V300`, `E100`, `CL100`, `G400`, `G99` for specific platform types.
    /// - For **Illumina**, it examines the first character of the instrument ID to infer the platform (e.g., `A` for NovaSeq, `D` for HiSeq 2500, etc.).
    ///
    /// # Notes
    /// - If no instrument data is available for a specific platform, a default generic name for the platform is returned, e.g., "PacBio", "MGI DNBseq",
    ///   or "Unknown Illumina".
    /// - `Unknown` is returned when the sequencing platform cannot be identified.
    pub fn infer_specific_platform(
        primary_platform: &SequencingPlatform,
        instruments: &HashMap<String, usize>,
    ) -> String {
        match primary_platform {
            SequencingPlatform::PacBio => {
                // PacBio instrument prefixes:
                // m64xxx = Sequel II/IIe
                // m54xxx = Sequel I
                // m84xxx = Revio
                if let Some((instrument_id, _)) = instruments.iter().max_by_key(|(_, &count)| count)
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
                if let Some((instrument_id, _)) = instruments.iter().max_by_key(|(_, &count)| count)
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
                if let Some((instrument_id, _)) = instruments.iter().max_by_key(|(_, &count)| count)
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
}
