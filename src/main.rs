use bio::alignment::{
    pairwise::{Aligner, Scoring},
    Alignment, AlignmentOperation,
};
use clap::Parser;

/// Strip the sequences from a PDB or mmCIF file.
#[derive(Parser, Debug)]
struct Args {
    /// first alignment
    #[arg()]
    x: String,

    /// Second sequence
    #[arg()]
    y: String,
}

fn main() {
    let args = Args::parse();
    if args.x.contains(',') {
        for (x, y) in args.x.split(',').zip(args.y.split(',')) {
            align(x.as_bytes(), y.as_bytes());
        }
    } else {
        align(args.x.as_bytes(), args.y.as_bytes());
    }
}

fn align(x: &[u8], y: &[u8]) {
    let scoring = get_blosum62(-12, -1);
    let mut aligner = Aligner::with_capacity_and_scoring(x.len(), y.len(), scoring);
    let alignment = aligner.global(x, y);
    let stats = score_stats(&alignment, x, y);
    println!("Identity: {}, Gaps: {}", stats.0, stats.1);
    //println!("{:#?}", alignment);
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

pub fn score_stats(alignment: &Alignment, sequence_x: &[u8], sequence_y: &[u8]) -> (f64, f64) {
    let x_len = sequence_x.len();
    let y_len = sequence_y.len();
    let mut x = alignment.xstart;
    let mut y = alignment.ystart;
    let mut identical = 0;
    let mut gaps = 0;
    for step in &alignment.operations {
        match step {
            AlignmentOperation::Del => {
                y += 1;
                gaps += 1;
            }
            AlignmentOperation::Ins => {
                x += 1;
                gaps += 1;
            }
            AlignmentOperation::Subst => {
                x += 1;
                y += 1;
            }
            AlignmentOperation::Match => {
                x += 1;
                y += 1;
                identical += 1;
            }
            AlignmentOperation::Xclip(_) => todo!(),
            AlignmentOperation::Yclip(_) => todo!(),
        }
    }
    debug_assert!(x == alignment.xend);
    debug_assert!(y == alignment.yend);
    (
        identical as f64 / (x_len).min(y_len) as f64,
        (gaps - (x_len as i32 - y_len as i32).abs()) as f64 / x_len.min(y_len) as f64,
    )
}
