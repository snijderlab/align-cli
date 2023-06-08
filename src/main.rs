use bio::alignment::pairwise::{Aligner, Scoring};
use clap::Parser;
use rustyms::{MonoIsotopic, Peptide};

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
    y: String,

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
    if args.mass {
        let a = Peptide::pro_forma(&args.x).unwrap();
        let b = Peptide::pro_forma(&args.y).unwrap();
        let ty = if args.local {
            rustyms::align::Type::Local
        } else if args.semi_global {
            rustyms::align::Type::GlobalForB
        } else {
            rustyms::align::Type::Global
        };
        let alignment = rustyms::align::align::<MonoIsotopic>(a, b, rustyms::align::BLOSUM62, ty);
        show_mass_alignment(&alignment, args.line_width);
    } else if args.x.contains(',') {
        for (x, y) in args.x.split(',').zip(args.y.split(',')) {
            align(&args, x.as_bytes(), y.as_bytes());
        }
    } else {
        align(&args, args.x.as_bytes(), args.y.as_bytes());
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
