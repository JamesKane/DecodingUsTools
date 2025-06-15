use crate::callable_loci::types::{CalledState, CoverageRange};
use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;
use tempfile::Builder;

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
        let attrs: String = self
            .attributes
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
    pub fn new(
        min_cutoff: u32,
        max_cutoff: u32,
        stride_len: u32,
        bar_height: u32,
        canvas_background: u32,
    ) -> Self {
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

            // Record the range based on its state by incrementing counts
            match range.state {
                CalledState::CALLABLE => {
                    for idx in start_idx..=end_idx {
                        if idx < array_size {
                            callable_depths[idx] += 1;
                        }
                    }
                }
                CalledState::POOR_MAPPING_QUALITY => {
                    for idx in start_idx..=end_idx {
                        if idx < array_size {
                            low_qual_depths[idx] += 1;
                        }
                    }
                }
                // Other states (NO_COVERAGE, LOW_COVERAGE, etc.) are not displayed in the histogram
                _ => {}
            }
        }

        (callable_depths, low_qual_depths)
    }

    fn generate_svg(
        &self,
        callable_depths: &[u32],
        low_qual_depths: &[u32],
        contig_name: &str,
    ) -> String {
        let svg_width = ((self.max_cutoff - self.min_cutoff) / self.stride_len) as u32;
        let histogram_height = self.bar_height;
        let title_height = 30;
        let header_padding = 15; // Padding between title and header content
        let label_height = 25; // Height for Mb labels
        let header_separator = 10; // Space between header and histogram
        let legend_height = 50;
        let total_header_height = title_height + header_padding + label_height + header_separator;
        let total_height = total_header_height + histogram_height + legend_height;
        let notch_height = 10;

        let mut svg =
            String::from("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");

        // Root SVG tag
        svg.push_str(
            &SvgTag::new("svg")
                .attr("xmlns", "http://www.w3.org/2000/svg")
                .attr("width", svg_width)
                .attr("height", total_height)
                .attr(
                    "style",
                    format!("background:#{:06x}", self.canvas_background),
                )
                .render(false),
        );
        svg.push('\n');

        // Title (contig name at the top)
        svg.push_str(
            &SvgTag::new("text")
                .attr("x", svg_width / 2)
                .attr("y", 20)
                .attr("text-anchor", "middle")
                .attr("font-family", "Arial")
                .attr("font-size", "16")
                .attr("font-weight", "bold")
                .attr("fill", "#000000")
                .render(false),
        );
        svg.push_str(contig_name);
        svg.push_str("</text>\n");

        // Header separator line
        svg.push_str(
            &SvgTag::new("line")
                .attr("x1", 0)
                .attr("y1", total_header_height - header_separator / 2)
                .attr("x2", svg_width)
                .attr("y2", total_header_height - header_separator / 2)
                .attr("stroke", "#808080")
                .attr("stroke-width", 1)
                .render(true),
        );
        svg.push('\n');

        // Header background
        svg.push_str(
            &SvgTag::new("rect")
                .attr("x", 0)
                .attr("y", 0)
                .attr("width", svg_width)
                .attr("height", total_header_height - header_separator)
                .attr("fill", "#F8F8F8")
                .attr("opacity", "0.8")
                .render(true),
        );
        svg.push('\n');

        // Draw position labels and notches
        let crawl_interval = 10_000_000; // 10Mb
        let positions: Vec<u32> = (0..=self.max_cutoff)
            .step_by(crawl_interval as usize)
            .collect();

        let label_y = title_height + header_padding + label_height - 5;

        for pos in positions {
            let x = ((pos - self.min_cutoff) / self.stride_len) as i32;

            if x >= svg_width as i32 || x < 0 {
                continue;
            }

            let is_label_fully_visible = x >= 20 && x <= (svg_width as i32 - 20);

            if is_label_fully_visible {
                let mb_pos = pos / 1_000_000;
                // Position label
                svg.push_str(
                    &SvgTag::new("text")
                        .attr("x", x)
                        .attr("y", label_y)
                        .attr("text-anchor", "middle")
                        .attr("font-family", "Arial")
                        .attr("font-size", "16")
                        .attr("font-weight", "bold")
                        .attr("fill", "#800080")
                        .render(false),
                );
                svg.push_str(&format!("{}Mb", mb_pos));
                svg.push_str("</text>\n");
            }

            // Top notch (at histogram top)
            svg.push_str(
                &SvgTag::new("line")
                    .attr("x1", x)
                    .attr("y1", total_header_height)
                    .attr("x2", x)
                    .attr("y2", total_header_height + notch_height)
                    .attr("stroke", "#800080")
                    .attr("stroke-width", 2)
                    .render(true),
            );
            svg.push('\n');

            // Bottom notch
            svg.push_str(
                &SvgTag::new("line")
                    .attr("x1", x)
                    .attr("y1", total_header_height + histogram_height - notch_height)
                    .attr("x2", x)
                    .attr("y2", total_header_height + histogram_height)
                    .attr("stroke", "#800080")
                    .attr("stroke-width", 2)
                    .render(true),
            );
            svg.push('\n');
        }

        // Plot histogram bars - adjusted for total_header_height
        for x in (self.min_cutoff..self.max_cutoff).step_by(self.stride_len as usize) {
            let idx = ((x - self.min_cutoff) / self.stride_len) as usize;
            let x_pos = idx;

            if callable_depths[idx] > 0 {
                let height = (callable_depths[idx] as f32 / 100.0 * histogram_height as f32) as u32;
                svg.push_str(
                    &SvgTag::new("rect")
                        .attr("x", x_pos)
                        .attr("y", total_header_height + histogram_height - height)
                        .attr("width", 1)
                        .attr("height", height)
                        .attr("fill", "#007700")
                        .render(true),
                );
                svg.push('\n');
            }

            if low_qual_depths[idx] > 0 {
                let height = (low_qual_depths[idx] as f32 / 100.0 * histogram_height as f32) as u32;
                let y_pos = total_header_height + histogram_height
                    - height
                    - (callable_depths[idx] as f32 / 100.0 * histogram_height as f32) as u32;
                svg.push_str(
                    &SvgTag::new("rect")
                        .attr("x", x_pos)
                        .attr("y", y_pos)
                        .attr("width", 1)
                        .attr("height", height)
                        .attr("fill", "#770000")
                        .render(true),
                );
                svg.push('\n');
            }
        }

        // Legend position adjusted
        let legend_y = total_header_height + histogram_height + 10;

        let legend_total_width = 300;
        let legend_start_x = (svg_width - legend_total_width) / 2;

        // Add gradients definitions before the legend
        svg.push_str("<defs>\n");
        svg.push_str(
            &SvgTag::new("linearGradient")
                .attr("id", "callableGradient")
                .attr("x1", "0%")
                .attr("y1", "0%")
                .attr("x2", "100%")
                .attr("y2", "0%")
                .attr("fill", "url(#callableGradient)")
                .render(false),
        );
        svg.push_str("<stop offset=\"0%\" style=\"stop-color:#007700;stop-opacity:0.8\"/>\n");
        svg.push_str("<stop offset=\"100%\" style=\"stop-color:#00aa00;stop-opacity:0.8\"/>\n");
        svg.push_str("</linearGradient>\n");

        svg.push_str(
            &SvgTag::new("linearGradient")
                .attr("id", "lowQualGradient")
                .attr("x1", "0%")
                .attr("y1", "0%")
                .attr("x2", "100%")
                .attr("y2", "0%")
                .attr("fill", "url(#lowQualGradient)")
                .render(false),
        );
        svg.push_str("<stop offset=\"0%\" style=\"stop-color:#770000;stop-opacity:0.8\"/>\n");
        svg.push_str("<stop offset=\"100%\" style=\"stop-color:#aa0000;stop-opacity:0.8\"/>\n");
        svg.push_str("</linearGradient>\n");
        svg.push_str("</defs>\n");

        // Callable coverage legend
        svg.push_str(
            &SvgTag::new("rect")
                .attr("x", legend_start_x)
                .attr("y", legend_y)
                .attr("width", 20)
                .attr("height", 10)
                .attr("fill", "url(#callableGradient)")
                .render(true),
        );
        svg.push_str(
            &SvgTag::new("text")
                .attr("x", legend_start_x + 25)
                .attr("y", legend_y + 8)
                .attr("font-family", "Arial")
                .attr("font-size", "12")
                .attr("fill", "#000000")
                .render(false),
        );
        svg.push_str("Callable Coverage");
        svg.push_str("</text>\n");

        // Low quality coverage legend
        svg.push_str(
            &SvgTag::new("rect")
                .attr("x", legend_start_x + 150)
                .attr("y", legend_y)
                .attr("width", 20)
                .attr("height", 10)
                .attr("fill", "url(#lowQualGradient)")
                .render(true),
        );
        svg.push_str(
            &SvgTag::new("text")
                .attr("x", legend_start_x + 175)
                .attr("y", legend_y + 8)
                .attr("font-family", "Arial")
                .attr("font-size", "12")
                .attr("fill", "#000000")
                .render(false),
        );
        svg.push_str("Low Quality Coverage");
        svg.push_str("</text>\n");

        // Close the group tag
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
    const MAX_WIDTH: u32 = 2000; // Maximum width for the largest contig
    const MIN_WIDTH: u32 = 200; // Minimum width for chrM
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
        0,             // min_cutoff
        contig_length, // max_cutoff uses actual contig length
        stride_len,    // adjusted stride
        100,           // bar_height
        0xFFFFFF,      // canvas_background (white)
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
