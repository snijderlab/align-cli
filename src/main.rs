use bio::alignment::pairwise::{Aligner, Scoring};
use clap::Parser;
use colored::Colorize;
use rustyms::ComplexPeptide;

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

    /// Use mass alignment, interpreting the sequences as ProForma
    #[arg(short, long)]
    mass: bool,

    /// The number of characters to show on a single line in the alignment
    #[arg(short = 'n', long, default_value_t = 50)]
    line_width: usize,
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
        if args.mass {
            let a = ComplexPeptide::pro_forma(&args.x);
            let b = ComplexPeptide::pro_forma(y);
            if let (Ok(a), Ok(b)) = (&a, &b) {
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
                    ty,
                );
                show_mass_alignment(&alignment, args.line_width);
            } else {
                println!("{}", "Error".red());
                if a.is_err() {
                    println!(
                        "Sequence A is not a valid Pro Forma sequence\nMessage: {}\n",
                        a.err().unwrap()
                    )
                }
                if b.is_err() {
                    println!(
                        "Sequence B is not a valid Pro Forma sequence\nMessage: {}\n",
                        b.err().unwrap()
                    )
                }
            }
        } else if args.x.contains(',') {
            for (x, y) in args.x.split(',').zip(y.split(',')) {
                align(&args, x.as_bytes(), y.as_bytes());
            }
        } else {
            align(&args, args.x.as_bytes(), y.as_bytes());
        }
    } else if let Some(path) = args.file {
        if !args.mass {
            panic!("Can only do the peptide to database matching based on mass")
        }
        let sequences = rustyms::identifications::FastaData::parse_reads(&path).unwrap();
        let mut alignments: Vec<_> = sequences
            .into_iter()
            .map(|seq| {
                let a = ComplexPeptide::pro_forma(&args.x);
                let b = seq.sequence.clone();
                if let Ok(a) = &a {
                    let ty = if args.local {
                        rustyms::align::Type::Local
                    } else if args.semi_global {
                        rustyms::align::Type::GlobalForB
                    } else {
                        rustyms::align::Type::Global
                    };
                    (
                        seq,
                        rustyms::align::align(
                            a.clone().assume_linear(),
                            b,
                            rustyms::align::BLOSUM62,
                            ty,
                        ),
                    )
                } else {
                    println!("{}", "Error".red());
                    if a.is_err() {
                        println!(
                            "Sequence A is not a valid Pro Forma sequence\nMessage: {}\n",
                            a.err().unwrap()
                        )
                    }
                    panic!()
                }
            })
            .collect();
        alignments.sort_unstable_by_key(|a| -a.1.score);
        let best = alignments[0].1.clone();
        let selected: Vec<_> = alignments.into_iter().take(10).collect();
        let mut data = vec![(
            "Rank".to_string(),
            "Database id".to_string(),
            "Score".to_string(),
            "Identity".to_string(),
            "Gap".to_string(),
        )];
        for (rank, (fasta, alignment)) in selected.into_iter().enumerate() {
            let stats = alignment.stats();
            data.push((
                (rank + 1).to_string(),
                fasta.id,
                alignment.score.to_string(),
                format!("{:.2}%", stats.0 as f64 / stats.2 as f64 * 100.0),
                format!("{:.2}%", stats.1 as f64 / stats.2 as f64 * 100.0),
            ));
        }
        let sizes = data
            .iter()
            .fold((0, 0, 0, 0, 0), |(aa, ab, ac, ad, ae), (a, b, c, d, e)| {
                (
                    aa.max(a.len()),
                    ab.max(b.len()),
                    ac.max(c.len()),
                    ad.max(d.len()),
                    ae.max(e.len()),
                )
            });
        println!(
            "┌{}┬{}┬{}┬{}┬{}┐",
            "─".repeat(sizes.0),
            "─".repeat(sizes.1),
            "─".repeat(sizes.2),
            "─".repeat(sizes.3),
            "─".repeat(sizes.4)
        );
        for (a, b, c, d, e) in data {
            println!(
                "│{:w0$}│{:w1$}│{:w2$}│{:w3$}│{:w4$}│",
                a,
                b,
                c,
                d,
                e,
                w0 = sizes.0,
                w1 = sizes.1,
                w2 = sizes.2,
                w3 = sizes.3,
                w4 = sizes.4
            );
        }
        println!(
            "└{}┴{}┴{}┴{}┴{}┘",
            "─".repeat(sizes.0),
            "─".repeat(sizes.1),
            "─".repeat(sizes.2),
            "─".repeat(sizes.3),
            "─".repeat(sizes.4)
        );
        println!("Alignment for the best match: ");
        show_mass_alignment(&best, args.line_width);
    } else {
        single_stats(
            &args,
            ComplexPeptide::pro_forma(&args.x).unwrap_or_else(|e| {
                panic!("Sequence is not a valid Pro Forma sequence\nMessage: {e}\n")
            }),
        )
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

fn single_stats(_args: &Args, seq: ComplexPeptide) {
    println!(
        "Mass: {}",
        seq.assume_linear()
            .formula()
            .and_then(|f| f.monoisotopic_mass().map(|m| format!("{:.2} Da", m.value)))
            .unwrap_or("Undefined".to_string())
    );
}
