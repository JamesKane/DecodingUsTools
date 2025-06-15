use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;
use tempfile::Builder;


#[derive(Debug)]
pub(crate) struct CoverageRange {
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) depth: u32,
    pub(crate) is_low_mapq: bool,
}

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
    fn new(min_cutoff: u32, max_cutoff: u32, stride_len: u32, bar_height: u32, canvas_background: u32) -> Self {
        Self {
            min_cutoff,
            max_cutoff,
            stride_len,
            bar_height,
            canvas_background,
        }
    }

    fn process_coverage_ranges(&self, ranges: &[CoverageRange]) -> (Vec<u32>, Vec<u32>) {
        let array_size = ((self.max_cutoff - self.min_cutoff) / self.stride_len) as usize;
        let mut callable_depths = vec![0; array_size];
        let mut low_qual_depths = vec![0; array_size];

        for range in ranges {
            let start_idx = ((range.start - self.min_cutoff) / self.stride_len) as usize;
            let end_idx = ((range.end - self.min_cutoff) / self.stride_len) as usize;

            for idx in start_idx..=end_idx {
                if idx < array_size {
                    if range.is_low_mapq {
                        low_qual_depths[idx] = range.depth;
                    } else {
                        callable_depths[idx] = range.depth;
                    }
                }
            }
        }

        (callable_depths, low_qual_depths)
    }

    fn generate_svg(&self, callable_depths: &[u32], low_qual_depths: &[u32], contig_name: &str) -> String {
        let svg_width = ((self.max_cutoff - self.min_cutoff) / self.stride_len) as u32;
        let svg_height = self.bar_height;

        let mut svg = String::from("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");

        // Root SVG tag
        svg.push_str(&SvgTag::new("svg")
            .attr("xmlns", "http://www.w3.org/2000/svg")
            .attr("width", svg_width)
            .attr("height", svg_height)
            .attr("style", format!("background:#{:06x}", self.canvas_background))
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

        // Add bar crawls (grid lines)
        let bar_stride = svg_width / (self.max_cutoff / 5_000_000);
        for x in (0..svg_width).step_by(bar_stride as usize) {
            for bar in (0..21).step_by(2) {
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
                        .render(true));
                    svg.push('\n');
                }
            }
        }

        svg.push_str("</svg>\n");
        svg
    }
}

pub(crate) fn generate_histogram(
    ranges: Vec<CoverageRange>,
    file_name_prefix: &str,
    contig_name: &str,
) -> std::io::Result<PathBuf> {
    let plotter = HistogramPlotter::new(
        0,      // min_cutoff
        60_000_000, // max_cutoff
        1000,   // stride_len
        400,    // bar_height
        0xFFFFFF // canvas_background (white)
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
