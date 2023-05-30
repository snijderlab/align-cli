use bio::alignment::{
    pairwise::{Aligner, Scoring},
    Alignment, AlignmentOperation,
};
use clap::Parser;
use colored::Colorize;
use std::fmt::Write;

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

    /// The number of characters to show on a single line in the alignment
    #[arg(short = 'n', long, default_value_t = 50)]
    line_width: usize,
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
    show_alignment(&alignment, x, y, args.semi_global, args.line_width);
}

fn show_alignment(
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

pub fn score_stats(
    alignment: &Alignment,
    sequence_x: &[u8],
    sequence_y: &[u8],
) -> (usize, usize, usize, usize) {
    let x_len = sequence_x.len();
    let y_len = sequence_y.len();
    let mut x = alignment.xstart;
    let mut y = alignment.ystart;
    let mut identical = 0;
    let mut similar = 0;
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
                if SIMILAR.contains(&(sequence_x[x], sequence_y[y])) {
                    similar += 1;
                }
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
    (identical, similar + identical, gaps, (x_len).max(y_len))
}

fn number_length(i: usize) -> usize {
    if i == 0 {
        1
    } else {
        i.ilog10() as usize + 1
    }
}

const SIMILAR: &[(u8, u8)] = &[(b'I', b'L'), (b'L', b'I'), (b'D', b'N'), (b'N', b'D')];

#[test]
fn number_length_test() {
    assert_eq!(number_length(0), 1);
    assert_eq!(number_length(1), 1);
    assert_eq!(number_length(9), 1);
    assert_eq!(number_length(10), 2);
    assert_eq!(number_length(11), 2);
    assert_eq!(number_length(99), 2);
    assert_eq!(number_length(100), 3);
    assert_eq!(number_length(1000), 4);
    assert_eq!(number_length(10000), 5);
}
