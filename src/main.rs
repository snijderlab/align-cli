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

    /// Use global alignment, default
    #[arg(short, long)]
    global: bool,

    /// Use semi-global alignment
    #[arg(short, long)]
    semi_global: bool,

    /// Use local alignment
    #[arg(short, long)]
    local: bool,
}

fn main() {
    let args = Args::parse();
    if args.global & args.semi_global || args.global && args.local || args.semi_global && args.local
    {
        panic!("Cannot have multiple alignment types at the same time")
    }
    if args.x.contains(',') {
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
    show_alignment(&alignment, x, y);
}

fn show_alignment(alignment: &Alignment, sequence_x: &[u8], sequence_y: &[u8]) {
    let stats = score_stats(alignment, sequence_x, sequence_y);

    println!(
        "Identity: {}, Gaps: {}, Score: {}, CIGAR: {}",
        stats.0,
        stats.1,
        alignment.score,
        alignment.cigar(false)
    );

    let mut top = Vec::new();
    let mut bottom = Vec::new();
    let mut x = alignment.xstart;
    let mut y = alignment.ystart;
    for step in &alignment.operations {
        match step {
            AlignmentOperation::Del => {
                top.push(b'-');
                bottom.push(sequence_y[y]);
                y += 1;
            }
            AlignmentOperation::Ins => {
                top.push(sequence_x[x]);
                bottom.push(b'-');
                x += 1;
            }
            AlignmentOperation::Subst => {
                top.push(sequence_x[x]);
                bottom.push(sequence_y[y]);
                x += 1;
                y += 1;
            }
            AlignmentOperation::Match => {
                top.push(sequence_x[x]);
                bottom.push(sequence_y[y]);
                x += 1;
                y += 1;
            }
            AlignmentOperation::Xclip(_) => todo!(),
            AlignmentOperation::Yclip(_) => todo!(),
        }
    }
    const CHUNK_SIZE: usize = 50;
    for section in top.chunks(CHUNK_SIZE).zip(bottom.chunks(CHUNK_SIZE)) {
        println!("{}", String::from_utf8_lossy(section.0));
        println!("{}\n", String::from_utf8_lossy(section.1));
    }
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
        identical as f64 / (x_len).max(y_len) as f64,
        gaps as f64 / x_len.max(y_len) as f64,
    )
}
