use bio::alignment::pairwise::{Aligner, Scoring};
use clap::Parser;
use colored::Colorize;
use rustyms::{
    find_isobaric_sets, ComplexPeptide, CustomError, LinearPeptide, MassTolerance, Modification,
    ReturnModification,
};
use std::{fmt::Display, io::Write, process::exit};

mod render;
mod stats;

use render::*;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about)]
struct Args {
    /// First sequence
    #[arg()]
    x: String,

    /// Second sequence
    #[arg()]
    y: Option<String>,

    /// A fasta database file to open to align the sequence to
    #[arg(short, long)]
    file: Option<String>,

    /// Use global alignment, default
    #[arg(short, long)]
    global: bool,

    /// Use semi-global alignment
    #[arg(short, long)]
    semi_global: bool,

    /// Use local alignment
    #[arg(short, long)]
    local: bool,

    /// Use normal alignment (instead of the default of Mass alignment)
    #[arg(long)]
    normal: bool,

    /// The number of characters to show on a single line in the alignment
    #[arg(short = 'n', long, default_value_t = 50)]
    line_width: usize,

    /// The maximal number of isobaric sets the generate, use `all` to generate all options
    #[arg(short, long, default_value_t = IsobaricNumber::Limited(25), value_parser=options_parse)]
    isobaric: IsobaricNumber,

    /// All possible modifications that will be used in the isobaric sets generation, separated by commas `,`
    #[arg(short, long, default_value_t = Modifications::None, value_parser=modifications_parse)]
    modifications: Modifications,

    /// The tolerance for the isobaric set search and the definition for isobaric sets in the alignment, use `<x>ppm` or `<x>da` to control the unit, eg `10.0ppm` or `2.3da`
    #[arg(short, long, default_value_t = MassTolerance::Ppm(10.0), value_parser=mass_tolerance_parse)]
    tolerance: MassTolerance,
}

#[derive(Debug, Clone)]
enum IsobaricNumber {
    All,
    Limited(usize),
}
impl Display for IsobaricNumber {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::All => write!(f, "all"),
            Self::Limited(limit) => write!(f, "{limit}"),
        }
    }
}
fn mass_tolerance_parse(input: &str) -> Result<MassTolerance, &'static str> {
    input.parse().map_err(|()| "Invalid tolerance parameter")
}
fn options_parse(input: &str) -> Result<IsobaricNumber, &'static str> {
    if input.to_lowercase() == "all" {
        Ok(IsobaricNumber::All)
    } else {
        input
            .parse::<usize>()
            .map(IsobaricNumber::Limited)
            .map_err(|_| "Invalid options parameter")
    }
}
#[derive(Debug, Clone)]
enum Modifications {
    None,
    Some(Vec<Modification>),
}
impl Modifications {
    fn mods(&self) -> &[Modification] {
        match self {
            Self::None => &[],
            Self::Some(m) => m,
        }
    }
}
impl Display for Modifications {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::None => write!(f, ""),
            Self::Some(mods) => {
                let mut start = true;
                for m in mods {
                    write!(f, "{}{}", if start { "" } else { "," }, m).unwrap();
                    start = false;
                }
                Ok(())
            }
        }
    }
}
fn modifications_parse(input: &str) -> Result<Modifications, &'static str> {
    if input.is_empty() {
        Ok(Modifications::None)
    } else {
        input
            .trim_end_matches(',')
            .split(',')
            .map(|m| {
                Modification::try_from(m, 0..m.len(), &mut Vec::new()).map(|m| match m {
                ReturnModification::Defined(d) => d,
                _ => {
                    panic!("Can not define ambiguous modifications for the modifications parameter")
                }
            })
            })
            .collect::<Result<Vec<Modification>, CustomError>>()
            .map(Modifications::Some)
            .map_err(|e| {
                print!("{e}");
                "Failed to parse modification"
            })
    }
}

fn main() {
    let args = Args::parse();
    if args.global as u8 + args.semi_global as u8 + args.local as u8 > 1 {
        panic!("Cannot have multiple alignment types at the same time")
    }
    if args.y.is_some() as u8 + args.file.is_some() as u8 > 1 {
        panic!("Cannot have multiple secondary sequence sources (y + file) at the same time")
    }
    if let Some(y) = &args.y {
        if !args.normal {
            let a = parse_peptide(&args.x);
            let b = parse_peptide(y);
            let ty = if args.local {
                rustyms::align::Type::Local
            } else if args.semi_global {
                rustyms::align::Type::GlobalForB
            } else {
                rustyms::align::Type::Global
            };
            let alignment = rustyms::align::align(
                a.clone().assume_linear(),
                b.clone().assume_linear(),
                rustyms::align::BLOSUM62,
                args.tolerance,
                ty,
            );
            show_mass_alignment(&alignment, args.line_width, args.tolerance);
        } else if args.x.contains(',') {
            for (x, y) in args.x.split(',').zip(y.split(',')) {
                align(&args, x.as_bytes(), y.as_bytes());
            }
        } else {
            align(&args, args.x.as_bytes(), y.as_bytes());
        }
    } else if let Some(path) = args.file {
        if args.normal {
            panic!("Can only do the peptide to database matching based on mass")
        }
        let sequences = rustyms::identifications::FastaData::parse_file(&path).unwrap();
        let search_sequence = parse_peptide(&args.x);
        let ty = if args.local {
            rustyms::align::Type::Local
        } else if args.semi_global {
            rustyms::align::Type::GlobalForB
        } else {
            rustyms::align::Type::Global
        };
        let mut alignments: Vec<_> = sequences
            .into_iter()
            .map(|seq| {
                let a = seq.peptide.clone();
                (
                    seq,
                    rustyms::align::align(
                        a,
                        search_sequence.clone().assume_linear(),
                        rustyms::align::BLOSUM62,
                        args.tolerance,
                        ty,
                    ),
                )
            })
            .collect();
        alignments.sort_unstable_by_key(|a| -a.1.absolute_score);
        let best = alignments[0].1.clone();
        let selected: Vec<_> = alignments.into_iter().take(10).collect();
        let mut data = vec![(
            "Rank".to_string(),
            "Database id".to_string(),
            "Score".to_string(),
            "Normalised score".to_string(),
            "Identity".to_string(),
            "Similarity".to_string(),
            "Gap".to_string(),
        )];
        for (rank, (fasta, alignment)) in selected.into_iter().enumerate() {
            let stats = alignment.stats();
            data.push((
                (rank + 1).to_string(),
                fasta.id,
                alignment.absolute_score.to_string(),
                format!("{:.3}", alignment.normalised_score),
                format!("{:.2}%", stats.0 as f64 / stats.3 as f64 * 100.0),
                format!("{:.2}%", stats.1 as f64 / stats.3 as f64 * 100.0),
                format!("{:.2}%", stats.2 as f64 / stats.3 as f64 * 100.0),
            ));
        }
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
        println!("Alignment for the best match: ");
        show_mass_alignment(&best, args.line_width, args.tolerance);
    } else {
        single_stats(
            &args,
            ComplexPeptide::pro_forma(&args.x)
                .unwrap_or_else(|e| {
                    panic!("Sequence is not a valid Pro Forma sequence\nMessage: {e}\n")
                })
                .assume_linear(),
        )
    }
}

/// Parses the peptide or shows an error and exits
fn parse_peptide(input: &str) -> ComplexPeptide {
    match ComplexPeptide::pro_forma(input) {
        Ok(v) => v,
        Err(e) => {
            println!("{}", e);
            exit(1);
        }
    }
}

fn align(args: &Args, x: &[u8], y: &[u8]) {
    let scoring = get_blosum62(-12, -1);
    let mut aligner = Aligner::with_capacity_and_scoring(x.len(), y.len(), scoring);
    let alignment = if args.local {
        aligner.local(x, y)
    } else if args.semi_global {
        aligner.semiglobal(x, y)
    } else {
        aligner.global(x, y)
    };
    show_alignment(&alignment, x, y, args.semi_global, args.line_width);
}

pub fn get_blosum62(gap_open: i32, gap_extend: i32) -> Scoring<impl Fn(u8, u8) -> i32> {
    const TRANSLATION_TABLE: &[usize] = include!("translation_table.txt");
    const BLOSUM62: &[&[i32]] = include!("blosum62.txt");
    let match_fn = |a: u8, b: u8| {
        BLOSUM62[if a > 64 && a < 91 {
            TRANSLATION_TABLE[(a - 65) as usize]
        } else {
            23 // Align as .
        }][if b > 64 && b < 91 {
            TRANSLATION_TABLE[(b - 65) as usize]
        } else {
            23 // Align as .
        }]
    };
    Scoring::new(gap_open, gap_extend, match_fn)
}

fn single_stats(args: &Args, seq: LinearPeptide) {
    if let Some(complete) = seq.formula().and_then(|f| f.monoisotopic_mass()) {
        let bare = seq.bare_formula().unwrap().monoisotopic_mass().unwrap();
        println!(
            "Full mass: {} Da | {} Da | {} Da {}",
            format!("{:.2}", complete.value).yellow(),
            format!(
                "{:.2}",
                seq.formula().unwrap().average_weight().unwrap().value
            )
            .yellow(),
            format!(
                "{:.2}",
                seq.formula().unwrap().most_abundant_mass().unwrap().value
            )
            .yellow(),
            "(monoisotopic | average | most abundant)".dimmed(),
        );
        println!(
            "Bare mass: {} Da {}",
            format!("{:.2}", bare.value).yellow(),
            "(no N/C terminal taken into account)".dimmed(),
        );
        println!(
            "Composition: {} {}",
            seq.bare_formula().unwrap().hill_notation_fancy().green(),
            "(no N/C terminal taken into account)".dimmed(),
        );
        if !matches!(args.isobaric, IsobaricNumber::Limited(0)) {
            match args.isobaric {
                IsobaricNumber::All => {
                    println!(
                        "Isobaric options {}: ",
                        format!("(all, tolerance {})", args.tolerance).dimmed()
                    );
                    let _ = std::io::stdout().flush();
                    for set in find_isobaric_sets(bare, args.tolerance, args.modifications.mods()) {
                        print!("{}, ", format!("{set}").blue());
                        let _ = std::io::stdout().flush();
                    }
                }
                IsobaricNumber::Limited(limit) => {
                    println!(
                        "Isobaric options: {}",
                        format!("(limited to {}, tolerance {})", limit, args.tolerance).dimmed()
                    );
                    let _ = std::io::stdout().flush();
                    for set in find_isobaric_sets(bare, args.tolerance, args.modifications.mods())
                        .take(limit)
                    {
                        print!("{}, ", format!("{set}").blue());
                        let _ = std::io::stdout().flush();
                    }
                }
            }
        }
    } else {
        println!("The sequence has no defined mass");
    }
}
