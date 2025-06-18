use crate::cli::Region;

impl Region {
    pub fn to_chromosome_names(&self) -> Vec<String> {
        match self {
            Region::Full => vec![],  // Empty vec means process all chromosomes
            Region::Ychr => vec!["chrY".to_string(), "Y".to_string()],  // Try both names
            Region::Mtchr => vec!["chrM".to_string(), "MT".to_string(), "M".to_string()],
        }
    }

    pub fn to_output_name(&self) -> &'static str {
        match self {
            Region::Full => "full",
            Region::Ychr => "chrY",  // Changed from "Ychr" to match BAM conventions
            Region::Mtchr => "chrM",
        }
    }
}
