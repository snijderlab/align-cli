use bio::alignment::{Alignment, AlignmentOperation};
use colored::{Color, ColoredString, Colorize};
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

trait Legend {
    fn fg_color(&self) -> Option<Color>;
    fn bg_color(&self) -> Option<Color>;
}

impl Legend for imgt_germlines::Annotation {
    fn fg_color(&self) -> Option<Color> {
        Some(match self {
            Self::Cysteine1 => Color::Yellow,
            Self::Cysteine2 => Color::Yellow,
            Self::Tryptophan => Color::Yellow,
            Self::Phenylalanine => Color::Yellow,
        })
    }
    fn bg_color(&self) -> Option<Color> {
        None
    }
}

impl Legend for imgt_germlines::Region {
    fn fg_color(&self) -> Option<Color> {
        match self {
            Self::CDR1 | Self::CDR2 | Self::CDR3 => Some(Color::Black),
            _ => None,
        }
    }
    fn bg_color(&self) -> Option<Color> {
        match self {
            Self::CDR1 => Some(Color::Red),
            Self::CDR2 => Some(Color::Green),
            Self::CDR3 => Some(Color::Blue),
            _ => None,
        }
    }
}

trait ExtendedColorize {
    fn style(self, style: &str) -> Self;
    fn on_color_e(self, color: Option<Color>) -> Self;
    fn color_e(self, color: Option<Color>) -> Self;
}

impl ExtendedColorize for ColoredString {
    fn style(self, style: &str) -> Self {
        match style {
            "underline" => self.underline(),
            _ => self,
        }
    }
    fn on_color_e(self, color: Option<Color>) -> Self {
        match color {
            Some(clr) => self.on_color(clr),
            None => self,
        }
    }
    fn color_e(self, color: Option<Color>) -> Self {
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
    imgt: Option<&Allele>,
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
        alignment.short().dimmed(),
    );

    let mut lines = (
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        0,
    );
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

    let mut number_shift_back = 1;
    let mut number_tail = String::new();
    let mut last_region = None;
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
            StepType::Insertion => (Some(Color::Yellow), "+"),
            StepType::Deletion => (Some(Color::Yellow), "+"),
            StepType::Match => (None, " "),
            StepType::Mismatch => (Some(Color::Red), "⨯"),
            StepType::Special => (Some(Color::Yellow), "-"), // ⇤⇥ ⤚---⤙ ├─┤ ║ ⤚⤙ l╴r╶
        };

        let region = imgt.and_then(|imgt| imgt.region(a + step.step_a as usize));
        let len = step.step_a.max(step.step_b) as usize;

        const NUMBER_GAP: usize = 10;
        if let Some(region) = region {
            if region.1 && number_tail.is_empty() && last_region != Some(region.0) {
                number_tail = format!("{} ", region.0).chars().rev().collect();

                // Compress the region from CDR3 to 3 or C3 depending on how much room is left
                let len = alignment.path[index..]
                    .iter()
                    .map(|a| a.step_a as usize)
                    .sum::<usize>();
                if len <= 1 {
                    number_tail = format!(" {}", number_tail.chars().nth(1).unwrap());
                } else if len <= 4 {
                    number_tail = format!(
                        " {}{}",
                        number_tail.chars().nth(1).unwrap(),
                        number_tail.chars().last().unwrap(),
                    );
                }
                last_region = Some(region.0);
            }
        }
        if (a + number_shift_back) % NUMBER_GAP == 0 && number_tail.is_empty() {
            let mut state = 0;
            let full_width = alignment
                .path
                .iter()
                .skip(index)
                .take_while(|s| {
                    state += s.step_a as usize;
                    state < number_shift_back
                })
                .map(|s| s.step_a.max(s.step_b) as usize)
                .sum::<usize>();
            number_tail = (a + number_shift_back).to_string();
            number_tail = format!(
                "{}{number_tail}",
                " ".repeat(full_width.saturating_sub(number_tail.len() - 1))
            )
            .chars()
            .rev()
            .collect();
            number_shift_back = (a + NUMBER_GAP + number_shift_back).to_string().len();
        }

        let a_str = if step.step_a == 0 {
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
        .chars()
        .collect::<Vec<_>>();
        let b_str = if step.step_b == 0 {
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
        .chars()
        .collect::<Vec<_>>();
        let bottom = if ty == StepType::Special {
            match len {
                1 => "─".to_string(),
                2 => "╶╴".to_string(),
                n => format!("╶{}╴", "─".repeat(n - 2)),
            }
        } else {
            ch.repeat(len)
        }
        .chars()
        .collect::<Vec<_>>();

        // Now write to the buffers one character at a time
        for s in 0..len {
            writelines(
                &mut lines,
                line_width,
                (
                    region.and_then(|r| r.0.fg_color()),
                    colour,
                    region.and_then(|r| r.0.bg_color()),
                ),
                number_tail.pop().unwrap_or(' '),
                (
                    a_str[s],
                    imgt.and_then(|imgt| {
                        imgt.annotations(a + s.max(step.step_a as usize))
                            .next()
                            .and_then(|a| a.fg_color())
                    }),
                    if alignment.seq_a.sequence[a..a + step.step_a as usize]
                        .iter()
                        .any(|a| !a.modifications.is_empty())
                    {
                        "underline"
                    } else {
                        "normal"
                    },
                ),
                (
                    b_str[s],
                    if alignment.seq_b.sequence[b..b + step.step_b as usize]
                        .iter()
                        .any(|a| !a.modifications.is_empty())
                    {
                        "underline"
                    } else {
                        "normal"
                    },
                ),
                bottom[s],
            )
        }
        a += step.step_a as usize;
        b += step.step_b as usize;
    }
    // Print any tail
    println!("{}", lines.0);
    println!("{}", lines.1);
    println!("{}", lines.2);
    println!("{}", lines.3);
}

fn writelines(
    lines: &mut (String, String, String, String, usize),
    line_width: usize,
    color: (Option<Color>, Option<Color>, Option<Color>),
    n: char,
    a: (char, Option<Color>, &str),
    b: (char, &str),
    c: char,
) {
    let color_fg = color.0.or(color.1);
    write!(
        &mut lines.0,
        "{}",
        n.to_string().normal().color_e(color.0).on_color_e(color.2)
    )
    .unwrap();
    write!(
        &mut lines.1,
        "{}",
        a.0.to_string()
            .normal()
            .color_e(a.1.or(color_fg))
            .on_color_e(color.2)
            .style(a.2)
    )
    .unwrap();
    write!(
        &mut lines.2,
        "{}",
        b.0.to_string()
            .normal()
            .color_e(color_fg)
            .on_color_e(color.2)
            .style(b.1)
    )
    .unwrap();
    write!(
        &mut lines.3,
        "{}",
        c.to_string().normal().color_e(color_fg).on_color_e(color.2)
    )
    .unwrap();
    lines.4 += 1;

    if lines.4 % line_width == 0 {
        println!("{}", lines.0);
        println!("{}", lines.1);
        println!("{}", lines.2);
        println!("{}", lines.3);
        lines.0.clear();
        lines.1.clear();
        lines.2.clear();
        lines.3.clear();
        lines.4 = 0;
    }
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
