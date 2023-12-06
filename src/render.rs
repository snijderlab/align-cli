use bio::alignment::{Alignment, AlignmentOperation};
use colored::{Color, Colorize, Styles};
use imgt_germlines::Allele;
use rustyms::{align::MatchType, MassTolerance};
use std::fmt::Write;

use crate::legend::*;
use crate::{stats::*, styling::*};

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

pub fn show_annotated_mass_alignment(
    alignment: &rustyms::align::Alignment,
    tolerance: MassTolerance,
    imgt: Option<&Allele>,
    line_width: usize,
    context: bool,
    only_display_a: bool,
) {
    let (identical, similar, gap, length) = alignment.stats();
    if !only_display_a {
        println!(
            "Identity: {} {}, Similarity: {} {}, Gaps: {} {}, Score: {} (Normalised: {}), Mass difference: {} Da {} ppm, {}\nStart A: {}, Start B: {}, Path: {}\n",
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
            alignment.start_a.to_string().magenta(),
            alignment.start_b.to_string().magenta(),
            alignment.short().dimmed(),
        );
    }

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

    const NUMBER_GAP: usize = 10;
    let mut writer = CombinedLines::new(line_width, only_display_a);
    let mut number_shift_back = 1;
    let mut number_tail = String::new();
    let mut is_number = false;
    let mut last_region = None;

    let mut header = |a: usize,
                      len: usize,
                      full_width: usize,
                      step: usize,
                      mut number_tail: String,
                      mut is_number,
                      mut number_shift_back: usize| {
        let region = imgt.and_then(|imgt| imgt.region(a + step));
        if let Some(region) = region {
            if region.1 && number_tail.is_empty() && last_region != Some(region.0) {
                number_tail = format!("{} ", region.0).chars().rev().collect();

                // Compress the region from CDR3 to 3 or C3 depending on how much room is left
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
                is_number = false;
            }
        }
        if (a + number_shift_back) % NUMBER_GAP == 0 && number_tail.is_empty() {
            number_tail = (a + number_shift_back).to_string();
            number_tail = format!(
                "{}{number_tail}",
                " ".repeat(full_width.saturating_sub(number_tail.len() - 1))
            )
            .chars()
            .rev()
            .collect();
            number_shift_back = (a + NUMBER_GAP + number_shift_back).to_string().len();
            is_number = true;
        }
        (number_tail, is_number, number_shift_back)
    };
    if context {
        let prefix = alignment.start_a.max(alignment.start_b);
        let shift_a = prefix - alignment.start_a;
        let shift_b = prefix - alignment.start_b;

        for index in 0..prefix {
            let a_index = (index >= shift_a).then_some(index.saturating_sub(shift_a));
            let b_index = (index >= shift_b).then_some(index.saturating_sub(shift_b));
            (number_tail, is_number, number_shift_back) = header(
                a_index.unwrap_or_default(),
                prefix - index,
                0,
                1,
                number_tail,
                is_number,
                number_shift_back,
            );

            writer.add_column(
                None,
                None,
                None,
                (
                    number_tail.pop().unwrap_or(' '),
                    Styling::with_style(Styles::Dimmed),
                ),
                (
                    a_index.map_or(' ', |a| alignment.seq_a.sequence[a].aminoacid.char()),
                    Styling::with_style(Styles::Dimmed).maybe_style(a_index.and_then(|a| {
                        alignment.seq_a.sequence[a..a + 1]
                            .iter()
                            .any(|a| !a.modifications.is_empty())
                            .then_some(Styles::Underline)
                    })),
                ),
                (
                    b_index.map_or(' ', |b| alignment.seq_b.sequence[b].aminoacid.char()),
                    Styling::with_style(Styles::Dimmed).maybe_style(b_index.and_then(|b| {
                        alignment.seq_b.sequence[b..b + 1]
                            .iter()
                            .any(|a| !a.modifications.is_empty())
                            .then_some(Styles::Underline)
                    })),
                ),
                ' ',
            );
        }
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
            StepType::Insertion => (Some(Color::Yellow), "+"),
            StepType::Deletion => (Some(Color::Yellow), "+"),
            StepType::Match => (None, " "),
            StepType::Mismatch => (Some(Color::Red), "⨯"),
            StepType::Special => (Some(Color::Yellow), "-"), // ⇤⇥ ⤚---⤙ ├─┤ ║ ⤚⤙ l╴r╶
        };

        let region = imgt.and_then(|imgt| imgt.region(a + step.step_a as usize));
        let len = step.step_a.max(step.step_b) as usize;

        let mut state = 0;
        (number_tail, is_number, number_shift_back) = header(
            a,
            alignment.path[index..]
                .iter()
                .map(|a| a.step_a as usize)
                .sum::<usize>(),
            alignment
                .path
                .iter()
                .skip(index)
                .take_while(|s| {
                    state += s.step_a as usize;
                    state < number_shift_back
                })
                .map(|s| s.step_a.max(s.step_b) as usize)
                .sum::<usize>(),
            step.step_a as usize,
            number_tail,
            is_number,
            number_shift_back,
        );

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
            writer.add_column(
                region.and_then(|r| r.0.fg_color()),
                colour,
                region.and_then(|r| r.0.bg_color()),
                (
                    number_tail.pop().unwrap_or(' '),
                    Styling::none().maybe_style(is_number.then_some(Styles::Dimmed)),
                ),
                (
                    a_str[s],
                    Styling::none()
                        .fg(imgt
                            .and_then(|imgt| imgt.annotations(a).next().and_then(|a| a.fg_color())))
                        .maybe_style(
                            (alignment.seq_a.sequence[a..a + step.step_a as usize]
                                .iter()
                                .any(|a| !a.modifications.is_empty()))
                            .then_some(Styles::Underline),
                        ),
                ),
                (
                    b_str[s],
                    Styling::none().maybe_style(
                        (alignment.seq_b.sequence[b..b + step.step_b as usize]
                            .iter()
                            .any(|a| !a.modifications.is_empty()))
                        .then_some(Styles::Underline),
                    ),
                ),
                bottom[s],
            )
        }
        a += step.step_a as usize;
        b += step.step_b as usize;
    }
    if context {
        let len = (alignment.seq_a.len() - a).max(alignment.seq_b.len() - b);

        for index in 0..len {
            let a_index = (a + index < alignment.seq_a.len()).then_some(a + index);
            let b_index = (b + index < alignment.seq_b.len()).then_some(b + index);
            (number_tail, is_number, number_shift_back) = header(
                a_index.unwrap_or_default(),
                len - index,
                0,
                1,
                number_tail,
                is_number,
                number_shift_back,
            );

            writer.add_column(
                None,
                None,
                None,
                (
                    number_tail.pop().unwrap_or(' '),
                    Styling::with_style(Styles::Dimmed),
                ),
                (
                    a_index.map_or(' ', |a| alignment.seq_a.sequence[a].aminoacid.char()),
                    Styling::with_style(Styles::Dimmed).maybe_style(a_index.and_then(|a| {
                        alignment.seq_a.sequence[a..a + 1]
                            .iter()
                            .any(|a| !a.modifications.is_empty())
                            .then_some(Styles::Underline)
                    })),
                ),
                (
                    b_index.map_or(' ', |b| alignment.seq_b.sequence[b].aminoacid.char()),
                    Styling::with_style(Styles::Dimmed).maybe_style(b_index.and_then(|b| {
                        alignment.seq_b.sequence[b..b + 1]
                            .iter()
                            .any(|a| !a.modifications.is_empty())
                            .then_some(Styles::Underline)
                    })),
                ),
                ' ',
            );
        }
    }

    // Print any tail
    writer.flush()
}

struct CombinedLines {
    numbers: String,
    numbers_content: bool,
    a: String,
    a_content: bool,
    b: String,
    b_content: bool,
    marker: String,
    marker_content: bool,
    chars: usize,
    line_width: usize,
    only_display_a: bool,
}

impl CombinedLines {
    fn new(line_width: usize, only_display_a: bool) -> Self {
        Self {
            numbers: String::with_capacity(line_width),
            numbers_content: false,
            a: String::with_capacity(line_width),
            a_content: false,
            b: String::with_capacity(line_width),
            b_content: false,
            marker: String::with_capacity(line_width),
            marker_content: false,
            chars: 0,
            line_width,
            only_display_a,
        }
    }

    fn add_column(
        &mut self,
        region_colour: Option<Color>,
        type_colour: Option<Color>,
        background_colour: Option<Color>,
        n: (char, Styling),
        a: (char, Styling),
        b: (char, Styling),
        c: char,
    ) {
        // Determine the foreground colour for the a/b/marker lines
        let color_fg = region_colour.or(type_colour);

        write!(
            &mut self.numbers,
            "{}",
            n.0.apply(&n.1.clone().fg(region_colour).bg(background_colour))
        )
        .unwrap();
        self.numbers_content |= !n.0.is_whitespace();

        write!(
            &mut self.a,
            "{}",
            a.0.apply(&a.1.clone().or_fg(color_fg).bg(background_colour))
        )
        .unwrap();
        self.a_content |= !a.0.is_whitespace();

        write!(
            &mut self.b,
            "{}",
            b.0.apply(&b.1.clone().fg(color_fg).bg(background_colour))
        )
        .unwrap();
        self.b_content |= !b.0.is_whitespace();

        write!(
            &mut self.marker,
            "{}",
            c.color_e(color_fg).on_color_e(background_colour)
        )
        .unwrap();
        self.marker_content |= !c.is_whitespace();

        // Flush if the maximal number of chars is reached
        self.chars += 1;
        if self.chars % self.line_width == 0 {
            self.flush()
        }
    }

    fn flush(&mut self) {
        // Only print a line if is has content
        if self.numbers_content {
            println!("{}", self.numbers);
        }
        if self.a_content {
            println!("{}", self.a);
        }
        if !self.only_display_a && self.b_content {
            println!("{}", self.b);
        }
        if !self.only_display_a && self.marker_content {
            println!("{}", self.marker);
        }
        // Reset all internal state
        self.numbers.clear();
        self.a.clear();
        self.b.clear();
        self.marker.clear();
        self.numbers_content = false;
        self.a_content = false;
        self.b_content = false;
        self.marker_content = false;
        self.chars = 0;
    }
}

pub fn table<const N: usize>(data: &[[String; N]], header: bool, styling: &[Styling; N]) {
    let sizes = data.iter().fold([0; N], |mut sizes, row| {
        for i in 0..N {
            sizes[i] = sizes[i].max(row[i].len());
        }
        sizes
    });
    print!("╭");
    for size in sizes.iter().take(N - 1).copied() {
        print!("{}┬", "─".repeat(size));
    }
    println!("{}╮", "─".repeat(sizes[N - 1]));
    if header {
        print!("│");
        #[allow(clippy::needless_range_loop)]
        for i in 0..N {
            print!("{:w$}│", data[0][i].blue(), w = sizes[i]);
        }
        println!();
        print!("├");
        for size in sizes.iter().take(N - 1).copied() {
            print!("{}┼", "─".repeat(size));
        }
        println!("{}┤", "─".repeat(sizes[N - 1]));
    }
    for row in data.iter().skip(usize::from(header)) {
        print!("│");
        for i in 0..N {
            print!("{:w$}│", row[i].apply(&styling[i]), w = sizes[i]);
        }
        println!();
    }
    print!("╰");
    for size in sizes.iter().take(N - 1).copied() {
        print!("{}┴", "─".repeat(size));
    }
    println!("{}╯", "─".repeat(sizes[N - 1]));
}
