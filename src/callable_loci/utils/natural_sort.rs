// Add a helper function for natural sorting of chromosome names
pub(crate) fn natural_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    // Special handling for chrM to always put it at the end
    if a == "chrM" {
        return std::cmp::Ordering::Greater;
    }
    if b == "chrM" {
        return std::cmp::Ordering::Less;
    }

    // Special handling for chrX and chrY
    if a == "chrX" && b == "chrY" {
        return std::cmp::Ordering::Less;
    }
    if a == "chrY" && b == "chrX" {
        return std::cmp::Ordering::Greater;
    }
    if a == "chrX" || a == "chrY" {
        return std::cmp::Ordering::Greater;
    }
    if b == "chrX" || b == "chrY" {
        return std::cmp::Ordering::Less;
    }

    // For other chromosomes, split into prefix and number
    let (a_prefix, a_num) = a.split_at(a.chars().take_while(|c| !c.is_ascii_digit()).count());
    let (b_prefix, b_num) = b.split_at(b.chars().take_while(|c| !c.is_ascii_digit()).count());

    // Compare prefixes first
    match a_prefix.cmp(b_prefix) {
        std::cmp::Ordering::Equal => {
            // If prefixes are equal, compare numbers
            let a_num = a_num.parse::<u32>().unwrap_or(0);
            let b_num = b_num.parse::<u32>().unwrap_or(0);
            a_num.cmp(&b_num)
        }
        other => other,
    }
}