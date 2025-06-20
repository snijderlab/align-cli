use colored::{Color, Colorize, Styles};
use itertools::Itertools;
use rustyms::{
    align::{Alignment, MatchType, Piece},
    imgt::Allele,
    prelude::{AminoAcid, IsAminoAcid, Peptidoform},
    quantities::Tolerance,
    sequence::{AnnotatedPeptide, Annotation, AtMax, Linear, Region},
    system::Mass,
};
use std::cmp::Ordering;
use std::collections::HashSet;
use std::fmt::Display;
use std::fmt::Write;

use crate::{legend::*, Cli};
use crate::{styling::*, NUMBER_PRECISION};

#[derive(PartialEq, Eq)]
enum StepType {
    Insertion,
    Deletion,
    Match,
    Mismatch,
    Special,
    MassMismatch,
}

pub fn show_annotated_mass_alignment<
    A: AtMax<Linear>,
    B: AtMax<Linear>,
    Annotated: AnnotatedPeptide,
>(
    alignment: &Alignment<'_, A, B>,
    imgt: Option<&Annotated>,
    only_display_a: bool,
    omit_headers: bool,
    line_names: (
        impl Into<String> + Display + Clone,
        impl Into<String> + Display + Clone,
    ),
    args: &Cli,
) {
    if !only_display_a {
        show_alignment_header(
            alignment,
            args.tolerance,
            line_names.clone(),
            None,
            args.full_number,
        );
    }
    let mut writer =
        CombinedLines::new(args.line_width, only_display_a, omit_headers, line_names.1);
    show_alignment_inner(
        &mut writer,
        alignment,
        imgt,
        args.context,
        None,
        false,
        String::new(),
        line_names.0.into(),
    );
    writer.flush();
}

pub fn show_chained_annotated_mass_alignment<A: AtMax<Linear>, B: AtMax<Linear>>(
    alignments: &[(Allele, Alignment<'_, A, B>)],
    tolerance: Tolerance<Mass>,
    line_width: usize,
    context: bool,
    full_number: bool,
    generate_annotation: bool,
) {
    let mut start = 0;
    for alignment in alignments {
        println!(
            "{}",
            format!("{} / {}", alignment.0.name(), alignment.0.fancy_name(),).purple(),
        );
        show_alignment_header(
            &alignment.1,
            tolerance,
            (alignment.0.name(), "Query"),
            Some(start),
            full_number,
        );
        start += alignment.1.len_b() + alignment.1.start_b();
    }

    let mut writer = CombinedLines::new(line_width, false, false, "Query");
    let mut number_tail = String::new();
    let mut last_context = None;
    for (index, alignment) in alignments.iter().enumerate() {
        (number_tail, last_context) = show_alignment_inner(
            &mut writer,
            &alignment.1,
            Some(&alignment.0),
            index == alignments.len() - 1 && context,
            last_context, // Original overwrite J with CDR3
            index != alignments.len() - 1,
            number_tail,
            alignment.0.name(),
        );
    }
    writer.flush();

    if generate_annotation {
        // Show annotation and regions for fasta
        // let mut annotations = Vec::new();
        let mut regions = Vec::new();
        let mut a_regions = alignments
            .iter()
            .map(|(a, al)| (a, (al.start_b() != 0).then_some((None, al.start_a()))))
            .flat_map(|(a, start)| {
                start
                    .into_iter()
                    .chain(a.regions.iter().map(|(r, l)| (Some(r.clone()), *l)))
            })
            .collect_vec(); // TODO: this misses unmatched regions between alignments
        a_regions.reverse();

        let mut len_a = 0;
        let mut len_b = 0;
        let mut last_region = None;
        for path in alignments
            .iter()
            .map(|(_, al)| {
                (
                    al,
                    (al.start_b() != 0).then_some(Piece {
                        score: 0,
                        local_score: 0,
                        match_type: MatchType::FullIdentity,
                        step_a: al.start_a() as u16,
                        step_b: al.start_b() as u16,
                    }),
                )
            })
            .flat_map(|(al, a)| a.into_iter().chain(al.path().iter().cloned()))
        {
            len_a += path.step_a as usize;
            len_b += path.step_b as usize;
            if let Some((r, l)) = a_regions.last().cloned() {
                if l <= len_a {
                    let region = r
                        .clone()
                        .or(last_region)
                        .unwrap_or(Region::Other("Unknown".to_string()));
                    if regions.last().is_some_and(|(r, _)| *r == region) {
                        regions.last_mut().unwrap().1 += len_b;
                    } else {
                        regions.push((region, len_b));
                    }
                    last_region = r.clone();
                    a_regions.pop();
                    len_a -= l;
                    len_b = 0;
                }
            }
        }
        // Map the remaining piece to the last element
        if let Some((r, _)) = a_regions.last().cloned() {
            let region = r
                .clone()
                .or(last_region)
                .unwrap_or(Region::Other("Unknown".to_string()));
            if regions.last().is_some_and(|(r, _)| *r == region) {
                regions.last_mut().unwrap().1 += len_b;
            } else {
                regions.push((region, len_b));
            }
        }

        println!(
            "REGIONS={}",
            regions.iter().map(|(r, l)| format!("{r}:{l}")).join(";")
        );
    }
}

fn show_alignment_inner<A, B, Annotated: AnnotatedPeptide>(
    writer: &mut CombinedLines,
    alignment: &Alignment<'_, A, B>,
    imgt: Option<&Annotated>,
    context: bool,
    start_context_override: Option<Region>,
    room_on_end: bool,
    number_tail: String,
    a_name: String,
) -> (String, Option<Region>) {
    let (mut a, mut b) = alignment.start();
    let a_glycan = find_possible_n_glycan_locations(alignment.seq_a());
    let b_glycan = find_possible_n_glycan_locations(alignment.seq_b());
    const NUMBER_GAP: usize = 10;
    let mut number_shift_back = 1;
    let mut number_tail = number_tail;
    let mut is_number = false;
    let mut last_region = start_context_override.as_ref();

    let mut header = |a: usize,
                      len: usize,
                      full_width: usize,
                      step: usize,
                      mut number_tail: String,
                      mut is_number,
                      mut number_shift_back: usize,
                      room_on_end: bool,
                      skip_first: bool| {
        let region = imgt.and_then(|imgt| imgt.get_region(a + step));
        if let Some(region) = region {
            if region.1
                && number_tail.is_empty()
                && last_region != Some(region.0)
                && !(skip_first && last_region == start_context_override.as_ref())
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
    // Start context
    if context || start_context_override.is_some() {
        let prefix = alignment.start_a().max(alignment.start_b());
        let shift_a = prefix - alignment.start_a();
        let shift_b = prefix - alignment.start_b();

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
                start_context_override.as_ref().is_some(),
            );
            let base_style = start_context_override
                .as_ref()
                .map(|_| Styling::none())
                .unwrap_or(Styling::with_style(Styles::Dimmed));

            writer.add_column(
                start_context_override
                    .as_ref()
                    .map(|_| "")
                    .unwrap_or(&a_name),
                start_context_override.as_ref().and_then(|r| r.fg_color()),
                None,
                start_context_override.as_ref().and_then(|r| r.bg_color()),
                (number_tail.pop().unwrap_or(' '), base_style.clone()),
                (
                    a_index.map_or(' ', |a| {
                        alignment.seq_a().sequence()[a]
                            .aminoacid
                            .one_letter_code()
                            .unwrap_or('X')
                    }),
                    base_style.clone().maybe_style(a_index.and_then(|a| {
                        alignment.seq_a().sequence()[a..a + 1]
                            .iter()
                            .any(|a| !a.modifications.is_empty())
                            .then_some(Styles::Underline)
                    })),
                ),
                (
                    b_index.map_or(' ', |b| {
                        alignment.seq_b().sequence()[b]
                            .aminoacid
                            .one_letter_code()
                            .unwrap_or('X')
                    }),
                    base_style.clone().maybe_style(b_index.and_then(|b| {
                        alignment.seq_b().sequence()[b..b + 1]
                            .iter()
                            .any(|a| !a.modifications.is_empty())
                            .then_some(Styles::Underline)
                    })),
                ),
                ' ',
            );
        }
    }
    // Actual alignment / middle
    for (index, step) in alignment.path().iter().enumerate() {
        let ty = match (step.match_type, step.step_a, step.step_b) {
            (MatchType::Isobaric, _, _) => StepType::Special, // Catch any 1/1 isobaric sets before they are counted as Match/Mismatch
            (MatchType::FullIdentity, _, _) => StepType::Match,
            (MatchType::IdentityMassMismatch, _, _) => StepType::MassMismatch,
            (MatchType::Mismatch, _, _) => StepType::Mismatch,
            (_, 0, 1) => StepType::Insertion,
            (_, 1, 0) => StepType::Deletion,
            _ => StepType::Special,
        };
        let (colour, ch) = match ty {
            StepType::Insertion => (Some(Color::Yellow), "+"),
            StepType::Deletion => (Some(Color::Yellow), "+"),
            StepType::Match => (None, " "),
            StepType::MassMismatch => (Some(Color::Yellow), "m"),
            StepType::Mismatch => (Some(Color::Red), "⨯"),
            StepType::Special => (Some(Color::Yellow), "-"), // ⇤⇥ ⤚---⤙ ├─┤ ║ ⤚⤙ l╴r╶
        };

        let region = imgt.and_then(|imgt| imgt.get_region(a + step.step_a as usize));
        let len = step.step_a.max(step.step_b) as usize;

        let mut state = 0;
        (number_tail, is_number, number_shift_back) = header(
            a,
            alignment.path()[index..]
                .iter()
                .map(|a| a.step_a as usize)
                .sum::<usize>(),
            alignment
                .path()
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
                alignment.seq_a()[a..a + step.step_a as usize]
                    .iter()
                    .map(|a| a.aminoacid.pro_forma_definition())
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
                alignment.seq_b()[b..b + step.step_b as usize]
                    .iter()
                    .map(|a| a.aminoacid.pro_forma_definition())
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
                        .fg(imgt.and_then(|imgt| {
                            imgt.get_annotations(a).next().and_then(|a| a.fg_color())
                        }))
                        .or_fg(
                            a_glycan
                                .contains(&a)
                                .then_some(Annotation::NGlycan)
                                .and_then(|a| a.fg_color()),
                        )
                        .maybe_style(
                            (alignment.seq_a()[a..a + step.step_a as usize]
                                .iter()
                                .any(|a| !a.modifications.is_empty()))
                            .then_some(Styles::Underline),
                        ),
                ),
                (
                    b_str[s],
                    Styling::none()
                        .fg(b_glycan
                            .contains(&b)
                            .then_some(Annotation::NGlycan)
                            .and_then(|a| a.fg_color()))
                        .maybe_style(
                            (alignment.seq_b()[b..b + step.step_b as usize]
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
    // End context
    if context {
        let len = (alignment.seq_a().len() - a).max(alignment.seq_b().len() - b);

        for index in 0..len {
            let a_index = (a + index < alignment.seq_a().len()).then_some(a + index);
            let b_index = (b + index < alignment.seq_b().len()).then_some(b + index);
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
                &a_name,
                None,
                None,
                None,
                (
                    number_tail.pop().unwrap_or(' '),
                    Styling::with_style(Styles::Dimmed),
                ),
                (
                    a_index.map_or(' ', |a| {
                        alignment.seq_a()[a]
                            .aminoacid
                            .one_letter_code()
                            .unwrap_or('X')
                    }),
                    Styling::with_style(Styles::Dimmed).maybe_style(a_index.and_then(|a| {
                        alignment.seq_a()[a..a + 1]
                            .iter()
                            .any(|a| !a.modifications.is_empty())
                            .then_some(Styles::Underline)
                    })),
                ),
                (
                    b_index.map_or(' ', |b| {
                        alignment.seq_b()[b]
                            .aminoacid
                            .one_letter_code()
                            .unwrap_or('X')
                    }),
                    Styling::with_style(Styles::Dimmed).maybe_style(b_index.and_then(|b| {
                        alignment.seq_b()[b..b + 1]
                            .iter()
                            .any(|a| !a.modifications.is_empty())
                            .then_some(Styles::Underline)
                    })),
                ),
                ' ',
            );
        }
    }
    (number_tail, last_region.cloned())
}

pub fn show_alignment_header<A: AtMax<Linear>, B: AtMax<Linear>>(
    alignment: &Alignment<'_, A, B>,
    tolerance: Tolerance<Mass>,
    names: (impl Display, impl Display),
    additional_b_start: Option<usize>,
    full_number: bool,
) {
    let precision = if full_number {
        None
    } else {
        Some(NUMBER_PRECISION)
    };
    let stats = alignment.stats();
    let score = alignment.score();
    println!(
        "Identity: {} {}, Mass similarity: {} {}, Similarity: {} {}, Gaps: {} {}, Score: {} {}, {}\nStart: {} {} {} {}, Path: {}\n{}\n",
        display_with_precision(stats.identity(), precision).bright_blue(),
        format!("({}/{})", stats.identical, stats.length).dimmed(),
        display_with_precision(stats.mass_similarity(), precision).blue(),
        format!("({}/{})", stats.mass_similar, stats.length).dimmed(),
        display_with_precision(stats.similarity(), precision).blue(),
        format!("({}/{})", stats.similar, stats.length).dimmed(),
        display_with_precision(stats.gaps_fraction(), precision).cyan(),
        format!("({}/{})", stats.gaps, stats.length).dimmed(),
        display_with_precision(score.normalised.0, precision).green(),
        format!("({}/{})", score.absolute, score.max).dimmed(),
        if alignment
        .mass_difference().value==0.0 {
            "Equal mass".yellow().to_string()
        } else {
            let (num, unit) = relative_notation(alignment.ppm().value * 1e6, 3); // ratio to ppm
            format!("Mass difference: {} {} {}",
                display_mass(alignment.mass_difference(), true, precision),
                num.yellow(),
                unit,)
        },
        names.0,
        alignment.start_a().to_string().magenta(),
        names.1,
        (additional_b_start.unwrap_or_default() + alignment.start_b()).to_string().magenta(),
        alignment.short().dimmed(),
        {
            format!("Tolerance: {tolerance}, Alignment: {} ({}), Maximal isobaric step: {}",
            alignment.align_type().description(),
            alignment.align_type().symbol(),
            alignment.max_step()).dimmed()
        },
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
    omit_headers: bool,
    a_names: HashSet<String>,
    b_name: String,
}

impl CombinedLines {
    fn new(
        line_width: usize,
        only_display_a: bool,
        omit_headers: bool,
        b_name: impl Into<String>,
    ) -> Self {
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
            omit_headers,
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
            b.0.apply(&b.1.clone().or_fg(color_fg).bg(background_colour))
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
        if !self.omit_headers {
            println!("{}", self.numbers);
        }
        let padding = if self.lines > 0 {
            " ".repeat(self.line_width - self.chars)
        } else {
            String::new()
        };
        if self.a_content {
            print!("{}", self.a,);
            if self.omit_headers {
                println!();
            } else {
                println!("{} {}", padding, self.a_names.iter().join(" / ").dimmed(),);
            }
        }
        if !self.only_display_a && self.b_content {
            print!("{}", self.b,);
            if self.omit_headers {
                println!();
            } else {
                println!("{} {}", padding, self.b_name.dimmed(),);
            }
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

pub fn table<const N: usize>(
    data: &[[String; N]],
    header: bool,
    styling: &[Styling; N],
    display_csv: bool,
) {
    if display_csv {
        let display_cell = |cell: &str| {
            let content = cell.replace('\"', "\'");
            if content.contains(',') {
                format!("\"{content}\"")
            } else {
                content
            }
        };
        for line in data {
            if !line.is_empty() {
                print!("{}", display_cell(&line[0]));
            }

            for cell in line.iter().skip(1) {
                print!(",{}", display_cell(cell))
            }
            println!();
        }
        println!();
    } else {
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
}

pub fn display_mass(value: Mass, colour: bool, precision: Option<usize>) -> String {
    let (num, suf) = engineering_notation(value.value, precision);
    format!(
        "{} {}Da",
        if colour {
            num.yellow().to_string()
        } else {
            num.to_string()
        },
        suf.map_or(String::new(), |s| s.to_string())
    )
}

fn display_with_precision(n: f64, precision: Option<usize>) -> String {
    if let Some(precision) = precision {
        format!("{n:.precision$}")
    } else {
        format!("{n}",)
    }
}

/// Display the given value in engineering notation eg `1000` -> `10 k`, with the given number of decimal points and returns the suffix separately.
/// /// A value of `0.0` will result in the lowest possible suffix `0.0 q`.
fn engineering_notation(value: f64, precision: Option<usize>) -> (String, Option<char>) {
    const BIG_SUFFIXES: &[char] = &[' ', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y', 'R', 'Q'];
    const SMALL_SUFFIXES: &[char] = &[' ', 'm', 'μ', 'n', 'p', 'f', 'a', 'z', 'y', 'r', 'q'];
    let base = if value == 0.0 {
        0
    } else {
        ((value.abs().log10() / 3.0).floor() as isize).clamp(
            -(SMALL_SUFFIXES.len() as isize - 1),
            BIG_SUFFIXES.len() as isize - 1,
        )
    };
    match base.cmp(&0) {
        Ordering::Less => (
            display_with_precision(
                value * 10_usize.pow(base.unsigned_abs() as u32 * 3) as f64,
                precision,
            ),
            Some(SMALL_SUFFIXES[base.unsigned_abs()]),
        ),
        Ordering::Equal => (display_with_precision(value, precision), None),
        Ordering::Greater => (
            display_with_precision(value / 10_usize.pow(base as u32 * 3) as f64, precision),
            Some(BIG_SUFFIXES[base as usize]),
        ),
    }
}

/// Display the given relative value in nice notation eg `1000 ppm` -> `10 ‰`, with the given number of decimal points and returns the suffix separately.
/// A value of `0.0` will result in the lowest possible suffix `0.0 ppq`.
fn relative_notation(ppm: f64, precision: usize) -> (String, &'static str) {
    if ppm >= 1.0e6 {
        (format!("{:.precision$}", ppm * 1e-6), "⅟")
    } else if ppm >= 1.0e3 {
        (format!("{:.precision$}", ppm * 1e-3), "‰")
    } else if ppm <= 1.0e-6 {
        (format!("{:.precision$}", ppm * 1e9), "ppq")
    } else if ppm <= 1.0e-3 {
        (format!("{:.precision$}", ppm * 1e6), "ppt")
    } else if ppm <= 1.0 {
        (format!("{:.precision$}", ppm * 1e3), "ppb")
    } else {
        (format!("{:.precision$}", ppm), "ppm")
    }
}

fn find_possible_n_glycan_locations<A>(sequence: &Peptidoform<A>) -> Vec<usize> {
    let mut result = Vec::new();
    for (index, aa) in sequence.sequence().windows(3).enumerate() {
        if let (AminoAcid::Asparagine, AminoAcid::Serine | AminoAcid::Threonine) =
            (aa[0].aminoacid.aminoacid(), aa[2].aminoacid.aminoacid())
        {
            if aa[1].aminoacid.aminoacid() != AminoAcid::Proline {
                result.push(index);
            }
        }
    }
    result
}
