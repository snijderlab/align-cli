use bio::alignment::{Alignment, AlignmentOperation};

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

pub fn number_length(i: usize) -> usize {
    if i == 0 {
        1
    } else {
        i.ilog10() as usize + 1
    }
}

pub const SIMILAR: &[(u8, u8)] = &[(b'I', b'L'), (b'L', b'I'), (b'D', b'N'), (b'N', b'D')];

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
