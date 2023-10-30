use bio::alignment::pairwise::{Aligner, Scoring};
use clap::{Args, Parser};
use colored::Colorize;
use rustyms::{
    error::CustomError,
    find_isobaric_sets,
    modification::{GnoComposition, ReturnModification},
    ontologies::*,
    placement_rule::*,
    AminoAcid, Chemical, ComplexPeptide, LinearPeptide, MassTolerance, Modification,
};
use std::{fmt::Display, io::Write, process::exit};

mod render;
mod stats;

use render::*;

#[derive(Parser, Debug)]
#[command(author, version, about)]
#[command(long_about = "It supports four distinct use cases:

1) Align two sequences `align <X> <Y>`, this shows the best alignment for these two sequences.

2) Align a single peptide to a database `align <X> --file <FILE.fasta>`, this shows the scores for the best matches for this peptide alongside the alignment for the best match.

3) Get information about a single sequence `align <sequence>`, this shows many basic properties (like mass) and generates isobaric sequences to this sequence.

4) Get information about a single modification, `align --modification <MODIFICATION>`, this shows basic properties, and if it is a mass shift, eg `+58.01`, it shows all modifications that have the same mass within the tolerance.")]
struct Cli {
    /// First sequence
    #[arg()]
    x: Option<String>,

    /// Second sequence
    #[arg()]
    y: Option<String>,

    /// A fasta database file to open to align the sequence to, only provide a single sequence for this mode
    #[arg(short, long)]
    file: Option<String>,

    #[command(flatten)]
    alignment_type: AlignmentType,

    /// Use normal alignment (instead of the default of Mass alignment) this uses Smith Waterman or Needleman Wunsch algorithms (based on the alignment mode)
    /// using the BLOSUM62 scoring table.
    #[arg(long)]
    normal: bool,

    /// The number of characters to show on a single line in the alignment
    #[arg(short = 'n', long, default_value_t = 50)]
    line_width: usize,

    /// The maximal number of isobaric sets the generate, use `all` to generate all options
    #[arg(short, long, default_value_t = IsobaricNumber::Limited(25), value_parser=options_parse)]
    isobaric: IsobaricNumber,

    /// All possible fixed modifications that will be used in the isobaric sets generation, separated by commas `,`
    #[arg(short = 'F', long, default_value_t = Modifications::None, value_parser=modifications_parse)]
    fixed: Modifications,

    /// All possible variable modifications that will be used in the isobaric sets generation, separated by commas `,`
    #[arg(short, long, default_value_t = Modifications::None, value_parser=modifications_parse)]
    variable: Modifications,

    /// The tolerance for the isobaric set search and the definition for isobaric sets in the alignment, use `<x>ppm` or `<x>da` to control the unit, eg `10.0ppm` or `2.3da`
    #[arg(short, long, default_value_t = MassTolerance::Ppm(10.0), value_parser=mass_tolerance_parse)]
    tolerance: MassTolerance,

    /// A modification you want details on, if it is a mass shift modification eg `+58.01` it will show all predefined modifications that are within the tolerance of this mass
    #[arg(short, long, value_parser=modification_parse)]
    modification: Option<Modification>,
}

#[test]
fn verify_cli() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}

#[derive(Args, Debug)]
#[group(multiple = false)]
struct AlignmentType {
    /// Use global alignment [default]
    #[arg(short, long)]
    global: bool,

    /// Use semi-global alignment, meaning that the second sequence has to match fully, while the first sequence can be longer then the alignment.
    /// When the `--file` mode is used this flag indicates that the given sequence can align semi globally to the provided database sequences.
    #[arg(short, long)]
    semi_global: bool,

    /// Use local alignment
    #[arg(short, long)]
    local: bool,
}

impl AlignmentType {
    fn ty(&self) -> rustyms::align::Type {
        if self.local {
            rustyms::align::Type::Local
        } else if self.semi_global {
            rustyms::align::Type::GlobalForB
        } else {
            rustyms::align::Type::Global
        }
    }
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

fn modification_parse(input: &str) -> Result<Modification, String> {
    if input.is_empty() {
        Err("Empty".to_string())
    } else {
        Modification::try_from(input, 0..input.len(), &mut Vec::new())
            .map(|m| match m {
                ReturnModification::Defined(d) => d,
                _ => {
                    panic!("Can not define ambiguous modifications for the modifications parameter")
                }
            })
            .map_err(|err| err.to_string())
    }
}

fn main() {
    let args = Cli::parse();
    if args.y.is_some() as u8 + args.file.is_some() as u8 > 1 {
        panic!("Cannot have multiple secondary sequence sources (y + file) at the same time")
    }
    if let (Some(x), Some(y)) = (&args.x, &args.y) {
        if !args.normal {
            let a = parse_peptide(x);
            let b = parse_peptide(y);
            let ty = args.alignment_type.ty();
            let alignment = rustyms::align::align(
                a.clone().assume_linear(),
                b.clone().assume_linear(),
                rustyms::align::BLOSUM62,
                args.tolerance,
                ty,
            );
            show_mass_alignment(&alignment, args.line_width, args.tolerance);
        } else if x.contains(',') {
            for (x, y) in x.split(',').zip(y.split(',')) {
                align(&args, x.as_bytes(), y.as_bytes());
            }
        } else {
            align(&args, x.as_bytes(), y.as_bytes());
        }
    } else if let (Some(x), Some(path)) = (&args.x, &args.file) {
        if args.normal {
            panic!("Can only do the peptide to database matching based on mass")
        }
        let sequences = rustyms::identifications::FastaData::parse_file(path).unwrap();
        let search_sequence = parse_peptide(x);
        let ty = args.alignment_type.ty();
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
    } else if let Some(x) = &args.x {
        single_stats(
            &args,
            ComplexPeptide::pro_forma(x)
                .unwrap_or_else(|e| {
                    panic!("Sequence is not a valid Pro Forma sequence\nMessage: {e}\n")
                })
                .assume_linear(),
        )
    } else if let Some(modification) = &args.modification {
        modification_stats(modification, args.tolerance);
    } else {
        println!("Please provide an argument to work with, use --help to see all options.")
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

fn align(args: &Cli, x: &[u8], y: &[u8]) {
    let scoring = get_blosum62(-12, -1);
    let mut aligner = Aligner::with_capacity_and_scoring(x.len(), y.len(), scoring);
    let alignment = if args.alignment_type.local {
        aligner.local(x, y)
    } else if args.alignment_type.semi_global {
        aligner.semiglobal(x, y)
    } else {
        aligner.global(x, y)
    };
    show_alignment(
        &alignment,
        x,
        y,
        args.alignment_type.semi_global,
        args.line_width,
    );
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

fn single_stats(args: &Cli, seq: LinearPeptide) {
    if let Some(complete) = seq.formula().and_then(|f| f.monoisotopic_mass()) {
        let bare = seq.bare_formula().unwrap().monoisotopic_mass().unwrap();
        println!(
            "Full mass: {} Da | {} Da {}",
            format!("{:.2}", complete.value).yellow(),
            format!(
                "{:.2}",
                seq.formula().unwrap().average_weight().unwrap().value
            )
            .yellow(),
            "(monoisotopic | average)".dimmed(),
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
                    for set in find_isobaric_sets(
                        bare,
                        args.tolerance,
                        args.fixed.mods(),
                        args.variable.mods(),
                    ) {
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
                    for set in find_isobaric_sets(
                        bare,
                        args.tolerance,
                        args.fixed.mods(),
                        args.variable.mods(),
                    )
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

fn modification_stats(modification: &Modification, tolerance: MassTolerance) {
    if let Modification::Mass(mass) = modification {
        println!(
            "All ontology modifications close to the given monoisotopic mass: {}",
            format!("tolerance: {tolerance}").dimmed()
        );
        for (_, _, modification) in rustyms::ontologies::unimod_ontology()
            .iter()
            .chain(rustyms::ontologies::psimod_ontology())
            .chain(rustyms::ontologies::gnome_ontology())
        {
            if tolerance.within(*mass, modification.formula().monoisotopic_mass().unwrap()) {
                println!(
                    "{} {} {}",
                    modification.to_string().purple(),
                    format!(
                        "{:.4}",
                        modification.formula().monoisotopic_mass().unwrap().value
                    )
                    .blue(),
                    modification.formula().hill_notation_fancy().green(),
                );
            }
        }
    } else if let Some(monoisotopic) = modification.formula().monoisotopic_mass() {
        println!(
            "Full mass: {} Da | {} Da  {}",
            format!("{:.2}", monoisotopic.value).yellow(),
            format!(
                "{:.2}",
                modification.formula().average_weight().unwrap().value
            )
            .yellow(),
            "(monoisotopic | average)".dimmed(),
        );
        println!(
            "Composition: {}",
            modification.formula().hill_notation_fancy().green(),
        );
        if let Modification::Predefined(_, rules, ontology, name, index) = modification {
            println!(
                "Ontology: {}, name: {}, index: {}",
                ontology.to_string().purple(),
                name.green(),
                index.to_string().blue()
            );
            print!("Placement rules: ");

            let mut first = true;
            for rule in rules {
                match rule {
                    PlacementRule::AminoAcid(aa, pos) => {
                        print!(
                            "{}{}@{}",
                            if first { "" } else { ", " },
                            aa.iter().map(AminoAcid::char).collect::<String>().yellow(),
                            pos.to_string().green()
                        )
                    }
                    PlacementRule::PsiModification(index, pos) => {
                        print!(
                            "{}{}@{}",
                            if first { "" } else { ", " },
                            psimod_ontology()
                                .find_id(*index)
                                .unwrap()
                                .to_string()
                                .blue(),
                            pos.to_string().green()
                        )
                    }
                    PlacementRule::Terminal(pos) => {
                        print!(
                            "{}{}",
                            if first { "" } else { ", " },
                            pos.to_string().green()
                        )
                    }
                }
                first = false;
            }
        } else if let Modification::Gno(composition, name) = modification {
            println!(
                "Ontology: {}, name: {}",
                "GNOme".to_string().purple(),
                name.to_uppercase().green(),
            );
            match composition {
                GnoComposition::Mass(_) => {
                    println!("Only mass known")
                }
                GnoComposition::Structure(structure) => {
                    println!("Structure: {}", structure.to_string().green())
                }
            }
        }
    } else {
        println!("{}", "No defined mass".red())
    }
}
