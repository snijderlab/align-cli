use colored::{Color, Colorize, Styles};
use imgt_germlines::{Allele, Region};
use itertools::Itertools;
use rustyms::{align::MatchType, MassTolerance};
use std::collections::HashSet;
use std::fmt::Display;
use std::fmt::Write;

use crate::legend::*;
use crate::styling::*;

#[derive(PartialEq, Eq)]
enum StepType {
    Insertion,
    Deletion,
    Match,
    Mismatch,
    Special,
}

pub fn show_annotated_mass_alignment(
    alignment: &rustyms::align::Alignment,
    tolerance: MassTolerance,
    imgt: Option<&Allele>,
    line_width: usize,
    context: bool,
    only_display_a: bool,
    line_names: (
        impl Into<String> + Display + Clone,
        impl Into<String> + Display + Clone,
    ),
) {
    if !only_display_a {
        show_alignment_header(alignment, tolerance, line_names.clone(), None);
    }
    let mut writer = CombinedLines::new(line_width, only_display_a, line_names.1);
    show_alignment_inner(
        &mut writer,
        alignment,
        imgt,
        context,
        None,
        false,
        String::new(),
        line_names.0.into(),
    );
    writer.flush();
}

pub fn show_chained_annotated_mass_alignment(
    alignment1: &(&Allele, rustyms::align::Alignment),
    alignment2: &(Allele, rustyms::align::Alignment),
    tolerance: MassTolerance,
    line_width: usize,
    context: bool,
) {
    println!(
        "{}: {}",
        "V gene".blue(),
        format!("{} / {}", alignment1.0.name(), alignment1.0.fancy_name(),).purple(),
    );
    show_alignment_header(
        &alignment1.1,
        tolerance,
        (alignment1.0.name(), "Query"),
        None,
    );
    println!(
        "{}: {}",
        "J gene".blue(),
        format!("{} / {}", alignment2.0.name(), alignment2.0.fancy_name(),).purple(),
    );
    show_alignment_header(
        &alignment2.1,
        tolerance,
        (alignment2.0.name(), "Query"),
        Some(alignment1.1.start_b + alignment1.1.len_b()),
    );

    let mut writer = CombinedLines::new(line_width, false, "Query");
    let number_tail = show_alignment_inner(
        &mut writer,
        &alignment1.1,
        Some(alignment1.0),
        context,
        None,
        true,
        String::new(),
        alignment1.0.name(),
    );
    show_alignment_inner(
        &mut writer,
        &alignment2.1,
        Some(&alignment2.0),
        context,
        Some(Region::CDR3),
        false,
        number_tail,
        alignment2.0.name(),
    );
    writer.flush();
}

fn show_alignment_inner(
    writer: &mut CombinedLines,
    alignment: &rustyms::align::Alignment,
    imgt: Option<&Allele>,
    context: bool,
    start_context_override: Option<Region>,
    room_on_end: bool,
    number_tail: String,
    a_name: String,
) -> String {
    let mut a = alignment.start_a;
    let mut b = alignment.start_b;
    const NUMBER_GAP: usize = 10;
    let mut number_shift_back = 1;
    let mut number_tail = number_tail;
    let mut is_number = false;
    let mut last_region = start_context_override;

    let mut header = |a: usize,
                      len: usize,
                      full_width: usize,
                      step: usize,
                      mut number_tail: String,
                      mut is_number,
                      mut number_shift_back: usize,
                      room_on_end: bool,
                      skip_first: bool| {
        let region = imgt.and_then(|imgt| imgt.region(a + step));
        if let Some(region) = region {
            if region.1
                && number_tail.is_empty()
                && last_region != Some(region.0)
                && !(skip_first && last_region == start_context_override)
            {
                number_tail = format!("{} ", region.0).chars().rev().collect();

                // Compress the region from CDR3 to 3 or C3 depending on how much room is left
                if !room_on_end && len == full_width {
                    if len <= 1 {
                        number_tail = format!(" {}", number_tail.chars().nth(1).unwrap());
                    } else if len <= 4 {
                        number_tail = format!(
                            " {}{}",
                            number_tail.chars().nth(1).unwrap(),
                            number_tail.chars().last().unwrap(),
                        );
                    }
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
    if context || start_context_override.is_some() {
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
                room_on_end,
                start_context_override.is_some(),
            );
            let base_style = start_context_override
                .map(|_| Styling::none())
                .unwrap_or(Styling::with_style(Styles::Dimmed));

            writer.add_column(
                "",
                start_context_override.and_then(|r| r.fg_color()),
                None,
                start_context_override.and_then(|r| r.bg_color()),
                (number_tail.pop().unwrap_or(' '), base_style.clone()),
                (
                    a_index.map_or(' ', |a| alignment.seq_a.sequence[a].aminoacid.char()),
                    base_style.clone().maybe_style(a_index.and_then(|a| {
                        alignment.seq_a.sequence[a..a + 1]
                            .iter()
                            .any(|a| !a.modifications.is_empty())
                            .then_some(Styles::Underline)
                    })),
                ),
                (
                    b_index.map_or(' ', |b| alignment.seq_b.sequence[b].aminoacid.char()),
                    base_style.clone().maybe_style(b_index.and_then(|b| {
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
            room_on_end,
            false,
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
                &a_name,
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
                room_on_end,
                false,
            );

            writer.add_column(
                "",
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
    number_tail
}

pub fn show_alignment_header(
    alignment: &rustyms::align::Alignment,
    tolerance: MassTolerance,
    names: (impl Display, impl Display),
    additional_b_start: Option<usize>,
) {
    let (identical, similar, gap, length) = alignment.stats();
    println!(
        "Identity: {} {}, Similarity: {} {}, Gaps: {} {}, Score: {} (Normalised: {}), Mass difference: {} Da {} ppm, {}\nStart: {} {} {} {}, Path: {}\n",
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
        format!("Tolerance: {}, Alignment: {:?}", tolerance, alignment.ty).dimmed(),
        names.0,
        alignment.start_a.to_string().magenta(),
        names.1,
        (additional_b_start.unwrap_or_default() + alignment.start_b).to_string().magenta(),
        alignment.short().dimmed(),
    );
}

struct CombinedLines {
    numbers: String,
    a: String,
    a_content: bool,
    b: String,
    b_content: bool,
    marker: String,
    marker_content: bool,
    chars: usize,
    lines: usize,
    line_width: usize,
    only_display_a: bool,
    a_names: HashSet<String>,
    b_name: String,
}

impl CombinedLines {
    fn new(line_width: usize, only_display_a: bool, b_name: impl Into<String>) -> Self {
        Self {
            numbers: String::with_capacity(line_width),
            a: String::with_capacity(line_width),
            a_content: false,
            b: String::with_capacity(line_width),
            b_content: false,
            marker: String::with_capacity(line_width),
            marker_content: false,
            chars: 0,
            lines: 0,
            line_width,
            only_display_a,
            a_names: HashSet::new(),
            b_name: b_name.into(),
        }
    }

    fn add_column(
        &mut self,
        a_name: &str,
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
        if !a_name.is_empty() {
            self.a_names.insert(a_name.to_string());
        }

        write!(
            &mut self.numbers,
            "{}",
            n.0.apply(&n.1.clone().fg(region_colour).bg(background_colour))
        )
        .unwrap();

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
        println!("{}", self.numbers);
        let padding = if self.lines > 0 {
            " ".repeat(self.line_width - self.chars)
        } else {
            String::new()
        };
        if self.a_content {
            println!(
                "{}{} {}",
                self.a,
                padding,
                self.a_names.iter().join(" / ").dimmed(),
            );
        }
        if !self.only_display_a && self.b_content {
            println!("{}{} {}", self.b, padding, self.b_name.dimmed(),);
        }
        if !self.only_display_a && self.marker_content {
            println!("{}", self.marker);
        }
        // Reset all internal state
        self.numbers.clear();
        self.a.clear();
        self.b.clear();
        self.marker.clear();
        self.a_content = false;
        self.b_content = false;
        self.marker_content = false;
        self.chars = 0;
        self.lines += 1;
        self.a_names.clear();
    }
}

pub fn table<const N: usize>(data: &[[String; N]], header: bool, styling: &[Styling; N]) {
    let sizes = data.iter().fold([0; N], |mut sizes, row| {
        for i in 0..N {
            sizes[i] = sizes[i].max(row[i].chars().count());
        }
        sizes
    });
    let line = |start, middle, end| {
        print!("{start}");
        for size in sizes.iter().take(N - 1).copied() {
            print!("{}{middle}", "─".repeat(size));
        }
        println!("{}{end}", "─".repeat(sizes[N - 1]));
    };
    line("╭", "┬", "╮");
    if header {
        print!("│");
        #[allow(clippy::needless_range_loop)]
        for i in 0..N {
            print!("{:^w$}│", data[0][i].blue(), w = sizes[i]);
        }
        println!();
        line("├", "┼", "┤");
    }
    for row in data.iter().skip(usize::from(header)) {
        print!("│");
        for i in 0..N {
            print!("{:w$}│", row[i].apply(&styling[i]), w = sizes[i]);
        }
        println!();
    }
    line("╰", "┴", "╯");
}
