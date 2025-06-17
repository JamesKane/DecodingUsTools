use crate::cli::Region;

impl Region {
    pub(crate) fn to_chromosome_names(&self) -> Vec<String> {
        match self {
            Region::Full => vec![], // Empty vec means process all chromosomes
            Region::Ychr => vec![
                "chrY".to_string(),
                "Y".to_string(),
                "NC_000024.10".to_string(),
                "CM000686.2".to_string(),
                "CP086569.2".to_string(),
                "NC_060948.1".to_string(),
            ],
            Region::Mtchr => vec![
                "chrM".to_string(),
                "MT".to_string(),
                "M".to_string(),
                "NC_012920.1".to_string(),
                "J01415.2".to_string(),
            ],
        }
    }
}