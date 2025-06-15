use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;
use tempfile::Builder;
use crate::callable_loci::types::{CalledState, CoverageRange};

struct HistogramPlotter {
    min_cutoff: u32,
    max_cutoff: u32,
    stride_len: u32,
    bar_height: u32,
    canvas_background: u32,
}

struct SvgTag {
    name: &'static str,
    attributes: HashMap<&'static str, String>,
}

impl SvgTag {
    fn new(name: &'static str) -> Self {
        Self {
            name,
            attributes: HashMap::new(),
        }
    }

    fn attr(mut self, key: &'static str, value: impl ToString) -> Self {
        self.attributes.insert(key, value.to_string());
        self
    }

    fn render(&self, self_closing: bool) -> String {
        let attrs: String = self.attributes
            .iter()
            .map(|(k, v)| format!("{}=\"{}\"", k, escape_xml_attr(v)))
            .collect::<Vec<_>>()
            .join(" ");

        if self_closing {
            format!("<{} {}/>", self.name, attrs)
        } else {
            format!("<{} {}>", self.name, attrs)
        }
    }
}

fn escape_xml_attr(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('\'', "&apos;")
        .replace('"', "&quot;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}


impl HistogramPlotter {
    pub fn new(min_cutoff: u32, max_cutoff: u32, stride_len: u32, bar_height: u32, canvas_background: u32) -> Self {
        Self {
            min_cutoff,
            max_cutoff,
            stride_len,
            bar_height,
            canvas_background,
        }
    }

    fn process_coverage_ranges(&self, ranges: &[CoverageRange]) -> (Vec<u32>, Vec<u32>) {
        let array_size = ((self.max_cutoff - self.min_cutoff) / self.stride_len) as usize + 1;
        let mut callable_depths = vec![0; array_size];
        let mut low_qual_depths = vec![0; array_size];

        for range in ranges {
            // Convert genomic coordinates to array indices
            let start_idx = ((range.start - self.min_cutoff) / self.stride_len) as usize;
            let end_idx = ((range.end - self.min_cutoff) / self.stride_len) as usize;

            // Record the range based on its state
            match range.state {
                CalledState::CALLABLE => {
                    for idx in start_idx..=end_idx {
                        if idx < array_size {
                            callable_depths[idx] = 100; // Using 100 as full height
                        }
                    }
                },
                CalledState::POOR_MAPPING_QUALITY => {
                    for idx in start_idx..=end_idx {
                        if idx < array_size {
                            low_qual_depths[idx] = 100; // Using 100 as full height
                        }
                    }
                },
                // Other states (NO_COVERAGE, LOW_COVERAGE, etc.) are not displayed in the histogram
                _ => {}
            }
        }

        (callable_depths, low_qual_depths)
    }


    fn generate_svg(&self, callable_depths: &[u32], low_qual_depths: &[u32], contig_name: &str) -> String {
        let svg_width = ((self.max_cutoff - self.min_cutoff) / self.stride_len) as u32;
        let svg_height = self.bar_height;
        // Add 30px padding at the top for the label
        let total_height = svg_height + 30;

        let mut svg = String::from("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");

        // Root SVG tag with increased height
        svg.push_str(&SvgTag::new("svg")
            .attr("xmlns", "http://www.w3.org/2000/svg")
            .attr("width", svg_width)
            .attr("height", total_height)
            .attr("style", format!("background:#{:06x}", self.canvas_background))
            .render(false));
        svg.push('\n');

        // Add contig name text
        svg.push_str(&SvgTag::new("text")
            .attr("x", svg_width / 2)
            .attr("y", 20)
            .attr("text-anchor", "middle")
            .attr("font-family", "Arial")
            .attr("font-size", "16")
            .attr("fill", "#000000")
            .render(false));
        svg.push_str(contig_name);
        svg.push_str("</text>\n");

        // Create a group for the histogram and translate it down by 30px
        svg.push_str(&SvgTag::new("g")
            .attr("transform", "translate(0,30)")
            .render(false));
        svg.push('\n');

        // Plot histograms
        for x in (self.min_cutoff..self.max_cutoff).step_by(self.stride_len as usize) {
            let idx = ((x - self.min_cutoff) / self.stride_len) as usize;
            let x_pos = idx;

            // Plot callable depth
            if callable_depths[idx] > 0 {
                let height = (callable_depths[idx] as f32 / 100.0 * svg_height as f32) as u32;
                svg.push_str(&SvgTag::new("rect")
                    .attr("x", x_pos)
                    .attr("y", svg_height - height)
                    .attr("width", 1)
                    .attr("height", height)
                    .attr("fill", "#007700")
                    .render(true));
                svg.push('\n');
            }

            // Plot low quality depth
            if low_qual_depths[idx] > 0 {
                let height = (low_qual_depths[idx] as f32 / 100.0 * svg_height as f32) as u32;
                let y_pos = svg_height - height - (callable_depths[idx] as f32 / 100.0 * svg_height as f32) as u32;
                svg.push_str(&SvgTag::new("rect")
                    .attr("x", x_pos)
                    .attr("y", y_pos)
                    .attr("width", 1)
                    .attr("height", height)
                    .attr("fill", "#770000")
                    .render(true));
                svg.push('\n');
            }
        }
        
        // Add bar crawls (grid lines) with labels
        let bar_stride = svg_width / (self.max_cutoff / 5_000_000);
        for x in (0..svg_width).step_by(bar_stride as usize) {
            // Add position label first with background
            let pos = (x * self.stride_len + self.min_cutoff) / 1_000_000;

            // Add a white background rectangle for the text
            svg.push_str(&SvgTag::new("rect")
                .attr("x", x - bar_stride/2)
                .attr("y", 0)
                .attr("width", bar_stride)
                .attr("height", 25)
                .attr("fill", format!("#{:06x}", self.canvas_background))
                .render(true));
            svg.push('\n');

            // Add the position label with larger font
            svg.push_str(&SvgTag::new("text")
                .attr("x", x)
                .attr("y", 18)
                .attr("text-anchor", "middle")
                .attr("font-family", "Arial")
                .attr("font-size", "16")
                .attr("font-weight", "bold")
                .attr("fill", "#800080")
                .render(false));
            svg.push_str(&format!("{}Mb", pos));
            svg.push_str("</text>\n");

            // Start bar crawls below the labels
            for bar in (2..21).step_by(2) {
                let bar_height = svg_height / 20;
                let y = bar * bar_height;
                for dx in -4..=4 {
                    svg.push_str(&SvgTag::new("line")
                        .attr("x1", x as i32 + dx)
                        .attr("y1", y)
                        .attr("x2", x as i32 + dx)
                        .attr("y2", y + 25)
                        .attr("stroke", "#800080")
                        .attr("stroke-width", 1)
                        .attr("stroke-opacity", "0.5")
                        .render(true));
                    svg.push('\n');
                }
            }
        }
        
        // Add gradients definitions
        svg.push_str("<defs>\n");
        svg.push_str(&SvgTag::new("linearGradient")
            .attr("id", "callableGradient")
            .attr("x1", "0%")
            .attr("y1", "0%")
            .attr("x2", "100%")
            .attr("y2", "0%")
            .attr("fill", "url(#callableGradient)")
            .render(false));
        svg.push_str("<stop offset=\"0%\" style=\"stop-color:#007700;stop-opacity:0.8\"/>\n");
        svg.push_str("<stop offset=\"100%\" style=\"stop-color:#00aa00;stop-opacity:0.8\"/>\n");
        svg.push_str("</linearGradient>\n");

        svg.push_str(&SvgTag::new("linearGradient")
            .attr("id", "lowQualGradient")
            .attr("x1", "0%")
            .attr("y1", "0%")
            .attr("x2", "100%")
            .attr("y2", "0%")
            .attr("fill", "url(#lowQualGradient)")
            .render(false));
        svg.push_str("<stop offset=\"0%\" style=\"stop-color:#770000;stop-opacity:0.8\"/>\n");
        svg.push_str("<stop offset=\"100%\" style=\"stop-color:#aa0000;stop-opacity:0.8\"/>\n");
        svg.push_str("</linearGradient>\n");
        svg.push_str("</defs>\n");

        // Add legend
        let legend_y = svg_height - 40;
        // Callable coverage legend
        svg.push_str(&SvgTag::new("rect")
            .attr("x", 10)
            .attr("y", legend_y)
            .attr("width", 20)
            .attr("height", 10)
            .attr("fill", "url(#callableGradient)")
            .render(true));
        svg.push_str(&SvgTag::new("text")
            .attr("x", 35)
            .attr("y", legend_y + 8)
            .attr("font-family", "Arial")
            .attr("font-size", "12")
            .attr("fill", "#000000")
            .render(false));
        svg.push_str("Callable Coverage");
        svg.push_str("</text>\n");

        // Low quality coverage legend
        svg.push_str(&SvgTag::new("rect")
            .attr("x", 150)
            .attr("y", legend_y)
            .attr("width", 20)
            .attr("height", 10)
            .attr("fill", "url(#lowQualGradient)")
            .render(true));
        svg.push_str(&SvgTag::new("text")
            .attr("x", 175)
            .attr("y", legend_y + 8)
            .attr("font-family", "Arial")
            .attr("font-size", "12")
            .attr("fill", "#000000")
            .render(false));
        svg.push_str("Low Quality Coverage");
        svg.push_str("</text>\n");


        // Close the group tag
        svg.push_str("</g>\n");
        svg.push_str("</svg>\n");
        svg
    }
}

pub(crate) fn generate_histogram(
    ranges: Vec<CoverageRange>,
    file_name_prefix: &str,
    contig_name: &str,
    contig_length: u32,
    largest_contig_length: u32,
) -> std::io::Result<PathBuf> {
    const MAX_WIDTH: u32 = 1000; // Maximum width for the largest contig
    const MIN_WIDTH: u32 = 200;  // Minimum width for chrM
    const CHRM_LENGTH: u32 = 16569; // Human mitochondrial genome length

    // Special handling for chrM
    let stride_len = if contig_name == "chrM" {
        // Use a smaller stride for chrM to ensure visibility
        // This will make it MIN_WIDTH pixels wide
        (CHRM_LENGTH + MIN_WIDTH - 1) / MIN_WIDTH
    } else {
        // Normal stride calculation for other chromosomes
        (largest_contig_length + MAX_WIDTH - 1) / MAX_WIDTH
    };

    let plotter = HistogramPlotter::new(
        0,              // min_cutoff
        contig_length,  // max_cutoff uses actual contig length
        stride_len,     // adjusted stride
        400,           // bar_height
        0xFFFFFF       // canvas_background (white)
    );

    let (callable_depths, low_qual_depths) = plotter.process_coverage_ranges(&ranges);

    let svg_content = plotter.generate_svg(&callable_depths, &low_qual_depths, contig_name);

    let temp_file = Builder::new()
        .prefix(file_name_prefix)
        .suffix(".svg")
        .tempfile()?;

    let mut file = temp_file.as_file();
    file.write_all(svg_content.as_bytes())?;

    Ok(temp_file.into_temp_path().keep()?.to_path_buf())
}
