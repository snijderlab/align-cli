use bio::alignment::{Alignment, AlignmentOperation};
use colored::{ColoredString, Colorize};
use imgt_germlines::Allele;
use rustyms::{align::MatchType, MassTolerance};
use std::fmt::Write;

use crate::stats::*;

pub fn show_alignment(
    alignment: &Alignment,
    sequence_x: &[u8],
    sequence_y: &[u8],
    semi_global: bool,
    line_width: usize,
) {
    let (identical, similar, gaps, length) = score_stats(alignment, sequence_x, sequence_y);

    println!(
        "Identity: {} {}, Similarity: {} {}, Gaps: {:} {}, Score: {}{}\n",
        format!("{:.3}", identical as f64 / length as f64).blue(),
        format!("({}/{})", identical, length).dimmed(),
        format!("{:.3}", similar as f64 / length as f64).cyan(),
        format!("({}/{})", similar, length).dimmed(),
        format!("{:.3}", gaps as f64 / length as f64).green(),
        format!("({}/{})", gaps, length).dimmed(),
        format!("{}", alignment.score).yellow(),
        if semi_global {
            format!("\nCIGAR: {}", alignment.cigar(false))
        } else {
            String::new()
        }
    );

    macro_rules! line {
        ($lines:ident, $x:expr, $y: expr, $c:expr, $colour:ident) => {
            write!(&mut $lines.0, "{}", $x.$colour()).unwrap();
            write!(&mut $lines.1, "{}", $y.$colour()).unwrap();
            write!(&mut $lines.2, "{}", $c.$colour()).unwrap();
        };
    }

    let mut lines = (String::new(), String::new(), String::new());
    let mut numbers = String::new();
    let mut x = alignment.xstart;
    let mut y = alignment.ystart;
    // Similar: · Gap: ◯ Identical: ∘ ⏺ ∘◯◦ ⚪ ⚫ ⬤ ⭘ 🞄 ∘ ○ ● ◦ ◯ ⴰ ⨉⨯+-
    for (index, step) in alignment.operations.iter().enumerate() {
        match step {
            AlignmentOperation::Del => {
                line!(
                    lines,
                    "-",
                    String::from_utf8_lossy(&[sequence_y[y]]),
                    "+",
                    yellow
                );
                y += 1;
            }
            AlignmentOperation::Ins => {
                line!(
                    lines,
                    String::from_utf8_lossy(&[sequence_x[x]]),
                    "-",
                    "+",
                    yellow
                );
                x += 1;
            }
            AlignmentOperation::Subst => {
                if SIMILAR.contains(&(sequence_x[x], sequence_y[y])) {
                    line!(
                        lines,
                        String::from_utf8_lossy(&[sequence_x[x]]),
                        String::from_utf8_lossy(&[sequence_y[y]]),
                        "-",
                        green
                    );
                } else {
                    line!(
                        lines,
                        String::from_utf8_lossy(&[sequence_x[x]]),
                        String::from_utf8_lossy(&[sequence_y[y]]),
                        "⨯",
                        red
                    );
                }
                x += 1;
                y += 1;
            }
            AlignmentOperation::Match => {
                line!(
                    lines,
                    String::from_utf8_lossy(&[sequence_x[x]]),
                    String::from_utf8_lossy(&[sequence_y[y]]),
                    " ",
                    normal
                );
                x += 1;
                y += 1;
            }
            AlignmentOperation::Xclip(_) => todo!(),
            AlignmentOperation::Yclip(_) => todo!(),
        }
        write!(&mut numbers, " ").unwrap();
        if (index + 1) % 10 == 0 {
            numbers.truncate(numbers.len() - number_length(index + 1));
            write!(&mut numbers, "{}", index + 1).unwrap();
        }
        if (index + 1) % line_width == 0 {
            println!("{}", numbers.dimmed());
            println!("{}", lines.0);
            println!("{}", lines.1);
            println!("{}", lines.2);
            lines = (String::new(), String::new(), String::new());
            numbers = String::new();
        }
    }
    println!("{}", numbers.dimmed());
    println!("{}", lines.0);
    println!("{}", lines.1);
    println!("{}", lines.2);
}

pub fn show_mass_alignment(
    alignment: &rustyms::align::Alignment,
    line_width: usize,
    tolerance: MassTolerance,
) {
    let (identical, similar, gap, length) = alignment.stats();

    println!(
        "Identity: {} {}, Similarity: {} {}, Gaps: {} {}, Score: {} (Normalised: {}), Mass difference: {} Da {} ppm, {}\nPath: {}\n",
        format!("{:.3}", identical as f64 / length as f64).bright_blue(),
        format!("({}/{})", identical, length).dimmed(),
        format!("{:.3}", similar as f64 / length as f64).blue(),
        format!("({}/{})", similar, length).dimmed(),
        format!("{:.3}", gap as f64 / length as f64).cyan(),
        format!("({}/{})", gap, length).dimmed(),
        format!("{}", alignment.absolute_score).green(),
        format!("{:.3}", alignment.normalised_score).green(),
        alignment
            .mass_difference()
            .map_or("??".to_string(), |m| format!("{:.2}", m.value))
            .yellow(),
        alignment
            .ppm()
            .map_or("??".to_string(), |m| format!("{:.2}", m))
            .yellow(),
        format!("Tolerance: {}", tolerance).dimmed(),
        alignment.short(),
    );

    let mut lines = (String::new(), String::new(), String::new());
    let mut numbers = String::new();
    let mut a = alignment.start_a;
    let mut b = alignment.start_b;

    #[derive(PartialEq, Eq)]
    enum StepType {
        Insertion,
        Deletion,
        Match,
        Mismatch,
        Special,
    }

    for (index, step) in alignment.path.iter().enumerate() {
        let ty = match (step.match_type, step.step_a, step.step_b) {
            (MatchType::Isobaric, _, _) => StepType::Special, // Catch any 1/1 isobaric sets before they are counted as Match/Mismatch
            (MatchType::FullIdentity, _, _) => StepType::Match,
            (MatchType::IdentityMassMismatch, _, _) => StepType::Match,
            (MatchType::Mismatch, _, _) => StepType::Mismatch,
            (_, 0, 1) => StepType::Insertion,
            (_, 1, 0) => StepType::Deletion,
            _ => StepType::Special,
        };
        let (colour, ch) = match ty {
            StepType::Insertion => ("yellow", "+"),
            StepType::Deletion => ("yellow", "+"),
            StepType::Match => ("white", " "),
            StepType::Mismatch => ("red", "⨯"),
            StepType::Special => ("yellow", "-"), // ⇤⇥ ⤚---⤙ ├─┤ ║ ⤚⤙ l╴r╶
        };
        let len = step.step_a.max(step.step_b) as usize;
        write!(
            &mut lines.0,
            "{}",
            if step.step_a == 0 {
                "-".repeat(len)
            } else {
                format!(
                    "{:·<width$}",
                    alignment.seq_a.sequence[a..a + step.step_a as usize]
                        .iter()
                        .map(|a| a.aminoacid.char())
                        .collect::<String>(),
                    width = len
                )
            }
            .color(
                if alignment.seq_a.sequence[a..a + step.step_a as usize]
                    .iter()
                    .any(|a| !a.modifications.is_empty())
                {
                    "blue"
                } else {
                    colour
                }
            )
        )
        .unwrap();
        write!(
            &mut lines.1,
            "{}",
            if step.step_b == 0 {
                "-".repeat(len)
            } else {
                format!(
                    "{:·<width$}",
                    alignment.seq_b.sequence[b..b + step.step_b as usize]
                        .iter()
                        .map(|a| a.aminoacid.char())
                        .collect::<String>(),
                    width = len
                )
            }
            .color(
                if alignment.seq_b.sequence[b..b + step.step_b as usize]
                    .iter()
                    .any(|a| !a.modifications.is_empty())
                {
                    "blue"
                } else {
                    colour
                }
            )
        )
        .unwrap();
        let bottom = if ty == StepType::Special {
            match len {
                1 => "─".to_string(),
                2 => "╶╴".to_string(),
                n => format!("╶{}╴", "─".repeat(n - 2)),
            }
        } else {
            ch.repeat(len)
        };
        write!(&mut lines.2, "{}", bottom.color(colour)).unwrap();

        a += step.step_a as usize;
        b += step.step_b as usize;

        write!(&mut numbers, " ").unwrap();
        if (index + 1) % 10 == 0 {
            numbers.truncate(numbers.len() - number_length(index + 1));
            write!(&mut numbers, "{}", index + 1).unwrap();
        }
        if (index + 1) % line_width == 0 {
            println!("{}", numbers.dimmed());
            println!("{}", lines.0);
            println!("{}", lines.1);
            println!("{}", lines.2);
            lines = (String::new(), String::new(), String::new());
            numbers = String::new();
        }
    }
    println!("{}", numbers.dimmed());
    println!("{}", lines.0);
    println!("{}", lines.1);
    println!("{}", lines.2);
}

trait Legend {
    fn fg_color(&self) -> Option<&'static str>;
    fn bg_color(&self) -> Option<&'static str>;
}

impl Legend for imgt_germlines::Annotation {
    fn fg_color(&self) -> Option<&'static str> {
        Some(match self {
            Self::Cysteine1 => "yellow",
            Self::Cysteine2 => "yellow",
            Self::Tryptophan => "yellow",
            Self::Phenylalanine => "yellow",
        })
    }
    fn bg_color(&self) -> Option<&'static str> {
        None
    }
}

impl Legend for imgt_germlines::Region {
    fn fg_color(&self) -> Option<&'static str> {
        match self {
            Self::CDR1 | Self::CDR2 | Self::CDR3 => Some("black"),
            _ => None,
        }
    }
    fn bg_color(&self) -> Option<&'static str> {
        match self {
            Self::CDR1 | Self::CDR2 | Self::CDR3 => Some("red"),
            _ => None,
        }
    }
}

trait ExtendedColorize {
    fn style(self, style: &str) -> Self;
    fn on_color_e(self, color: Option<&str>) -> Self;
    fn color_e(self, color: Option<&str>) -> Self;
}

impl ExtendedColorize for ColoredString {
    fn style(self, style: &str) -> Self {
        match style {
            "underline" => self.underline(),
            _ => self,
        }
    }
    fn on_color_e(self, color: Option<&str>) -> Self {
        match color {
            Some(clr) => self.on_color(clr),
            None => self,
        }
    }
    fn color_e(self, color: Option<&str>) -> Self {
        match color {
            Some(clr) => self.color(clr),
            None => self,
        }
    }
}

pub fn show_annotated_mass_alignment(
    alignment: &rustyms::align::Alignment,
    line_width: usize,
    tolerance: MassTolerance,
    imgt: &Allele,
) {
    let (identical, similar, gap, length) = alignment.stats();

    println!(
        "IMGT: {}\nIdentity: {} {}, Similarity: {} {}, Gaps: {} {}, Score: {} (Normalised: {}), Mass difference: {} Da {} ppm, {}\nPath: {}\n",
        imgt.name().purple(),
        format!("{:.3}", identical as f64 / length as f64).bright_blue(),
        format!("({}/{})", identical, length).dimmed(),
        format!("{:.3}", similar as f64 / length as f64).blue(),
        format!("({}/{})", similar, length).dimmed(),
        format!("{:.3}", gap as f64 / length as f64).cyan(),
        format!("({}/{})", gap, length).dimmed(),
        format!("{}", alignment.absolute_score).green(),
        format!("{:.3}", alignment.normalised_score).green(),
        alignment
            .mass_difference()
            .map_or("??".to_string(), |m| format!("{:.2}", m.value))
            .yellow(),
        alignment
            .ppm()
            .map_or("??".to_string(), |m| format!("{:.2}", m))
            .yellow(),
        format!("Tolerance: {}", tolerance).dimmed(),
        alignment.short(),
    );

    let mut lines = (String::new(), String::new(), String::new(), String::new());
    let mut a = alignment.start_a;
    let mut b = alignment.start_b;

    #[derive(PartialEq, Eq)]
    enum StepType {
        Insertion,
        Deletion,
        Match,
        Mismatch,
        Special,
    }

    let mut number_tail = String::new();
    for (index, step) in alignment.path.iter().enumerate() {
        let ty = match (step.match_type, step.step_a, step.step_b) {
            (MatchType::Isobaric, _, _) => StepType::Special, // Catch any 1/1 isobaric sets before they are counted as Match/Mismatch
            (MatchType::FullIdentity, _, _) => StepType::Match,
            (MatchType::IdentityMassMismatch, _, _) => StepType::Match,
            (MatchType::Mismatch, _, _) => StepType::Mismatch,
            (_, 0, 1) => StepType::Insertion,
            (_, 1, 0) => StepType::Deletion,
            _ => StepType::Special,
        };
        let (colour, ch) = match ty {
            StepType::Insertion => ("yellow", "+"),
            StepType::Deletion => ("yellow", "+"),
            StepType::Match => ("white", " "),
            StepType::Mismatch => ("red", "⨯"),
            StepType::Special => ("yellow", "-"), // ⇤⇥ ⤚---⤙ ├─┤ ║ ⤚⤙ l╴r╶
        };
        // TODO: it uses step_a steps, while these regions and stuff want a single char... (also fixes the weird alignment if it is done based on a single char)
        let region = imgt.region(a + step.step_a as usize).unwrap();
        let annotations = imgt.annotations(a + step.step_a as usize).next(); // Only interested in the first annotation
        let len = step.step_a.max(step.step_b) as usize;

        if a % 10 == 0 && number_tail.is_empty() {
            number_tail = a.to_string().chars().rev().collect();
        }
        if region.1 {
            number_tail = region.0.to_string().chars().rev().collect();

            // Compress the region from CDR3 to 3 or C3 depending on how much room is left
            let len = alignment.path[index..].iter().map(|a| a.step_a).sum::<u8>();
            if len <= 1 {
                number_tail = number_tail.chars().next().unwrap().to_string();
            } else if len <= 4 {
                number_tail = format!(
                    "{}{}",
                    number_tail.chars().next().unwrap(),
                    number_tail.chars().last().unwrap(),
                );
            }
        }
        write!(
            &mut lines.0,
            "{}",
            number_tail
                .pop()
                .unwrap_or(' ')
                .to_string()
                .normal()
                .color_e(region.0.fg_color())
                .on_color_e(region.0.bg_color())
        )
        .unwrap();
        write!(
            &mut lines.1,
            "{}",
            if step.step_a == 0 {
                "-".repeat(len)
            } else {
                format!(
                    "{:·<width$}",
                    alignment.seq_a.sequence[a..a + step.step_a as usize]
                        .iter()
                        .map(|a| a.aminoacid.char())
                        .collect::<String>(),
                    width = len
                )
            }
            .normal()
            .style(
                if alignment.seq_a.sequence[a..a + step.step_a as usize]
                    .iter()
                    .any(|a| !a.modifications.is_empty())
                {
                    "underline"
                } else {
                    "normal"
                }
            )
            .on_color_e(region.0.bg_color())
            .color(
                annotations
                    .and_then(|a| a.fg_color())
                    .or(region.0.fg_color())
                    .unwrap_or(colour)
            )
        )
        .unwrap();
        write!(
            &mut lines.2,
            "{}",
            if step.step_b == 0 {
                "-".repeat(len)
            } else {
                format!(
                    "{:·<width$}",
                    alignment.seq_b.sequence[b..b + step.step_b as usize]
                        .iter()
                        .map(|a| a.aminoacid.char())
                        .collect::<String>(),
                    width = len
                )
            }
            .normal()
            .style(
                if alignment.seq_b.sequence[b..b + step.step_b as usize]
                    .iter()
                    .any(|a| !a.modifications.is_empty())
                {
                    "underline"
                } else {
                    "normal"
                }
            )
            .on_color_e(region.0.bg_color())
            .color(region.0.fg_color().unwrap_or(colour))
        )
        .unwrap();
        let bottom = if ty == StepType::Special {
            match len {
                1 => "─".to_string(),
                2 => "╶╴".to_string(),
                n => format!("╶{}╴", "─".repeat(n - 2)),
            }
        } else {
            ch.repeat(len)
        };
        write!(
            &mut lines.3,
            "{}",
            bottom
                .normal()
                .on_color_e(region.0.bg_color())
                .color(region.0.fg_color().unwrap_or(colour))
        )
        .unwrap();

        a += step.step_a as usize;
        b += step.step_b as usize;

        if (index + 1) % line_width == 0 {
            println!("{}", lines.0);
            println!("{}", lines.1);
            println!("{}", lines.2);
            println!("{}", lines.3);
            lines = (String::new(), String::new(), String::new(), String::new());
        }
    }
    println!("{}", lines.0);
    println!("{}", lines.1);
    println!("{}", lines.2);
    println!("{}", lines.3);
}

pub fn table(data: &[(String, String, String, String, String, String, String)]) {
    let sizes = data.iter().fold(
        (0, 0, 0, 0, 0, 0, 0),
        |(aa, ab, ac, ad, ae, af, ag), (a, b, c, d, e, f, g)| {
            (
                aa.max(a.len()),
                ab.max(b.len()),
                ac.max(c.len()),
                ad.max(d.len()),
                ae.max(e.len()),
                af.max(f.len()),
                ag.max(g.len()),
            )
        },
    );
    println!(
        "┌{}┬{}┬{}┬{}┬{}┬{}┬{}┐",
        "─".repeat(sizes.0),
        "─".repeat(sizes.1),
        "─".repeat(sizes.2),
        "─".repeat(sizes.3),
        "─".repeat(sizes.4),
        "─".repeat(sizes.5),
        "─".repeat(sizes.6),
    );
    for (a, b, c, d, e, f, g) in data {
        println!(
            "│{:w0$}│{:w1$}│{:w2$}│{:w3$}│{:w4$}│{:w5$}│{:w6$}│",
            a,
            b,
            c,
            d,
            e,
            f,
            g,
            w0 = sizes.0,
            w1 = sizes.1,
            w2 = sizes.2,
            w3 = sizes.3,
            w4 = sizes.4,
            w5 = sizes.5,
            w6 = sizes.6,
        );
    }
    println!(
        "└{}┴{}┴{}┴{}┴{}┴{}┴{}┘",
        "─".repeat(sizes.0),
        "─".repeat(sizes.1),
        "─".repeat(sizes.2),
        "─".repeat(sizes.3),
        "─".repeat(sizes.4),
        "─".repeat(sizes.5),
        "─".repeat(sizes.6),
    );
}
