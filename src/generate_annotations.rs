use imgt::Allele;
use mzalign::{Alignment, MatchType, Piece};
use mzcore::sequence::{Annotation, Region};

pub fn generate_annotations<A, B>(
    alignments: &[(Allele, Alignment<A, B>)],
) -> (Vec<(Region, usize)>, Vec<(Annotation, usize)>) {
    // Show annotation and regions for fasta
    let mut annotations = Vec::new();
    let mut regions = Vec::new();
    let mut a_regions: Vec<_> = alignments
        .iter()
        .map(|(a, al)| (a, (al.start_b() != 0).then_some((None, al.start_a()))))
        .flat_map(|(a, start)| {
            start
                .into_iter()
                .chain(a.regions.iter().map(|(r, l)| (Some(r.clone()), *l)))
        })
        .collect();
    a_regions.reverse();

    let mut index_a = 0;
    let mut index_b = 0;
    let mut len_a = 0;
    let mut len_b = 0;
    let mut offset_a = 0;
    let mut last_region = None;
    let mut last_index = None;
    for (index, path) in alignments
        .iter()
        .enumerate()
        .map(|(i, (_, al))| {
            (
                i,
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
        .flat_map(|(i, al, a)| {
            a.into_iter()
                .chain(al.path().iter().cloned())
                .map(move |p| (i, p))
        })
    {
        if last_index != Some(index) {
            last_index = Some(index);
            offset_a = index_a;
        }
        index_a += path.step_a as usize;
        index_b += path.step_b as usize;
        len_a += path.step_a as usize;
        len_b += path.step_b as usize;
        if let Some((r, l)) = a_regions.last().cloned()
            && l <= len_a
        {
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

        // Only transfer annotations if the sequence is the same
        // Annotations are either conserved or N-glycan so the amino acid has to stay the same for them to make sense
        if matches!(
            path.match_type,
            MatchType::FullIdentity | MatchType::IdentityMassMismatch
        ) {
            for (a, annotation_index) in alignments[index].0.annotations {
                if index_a - offset_a + path.step_a as usize == *annotation_index {
                    annotations.push((a.clone(), index_b + path.step_b as usize));
                }
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

    (regions, annotations)
}

#[cfg(test)]
mod tests {
    use imgt::Gene;
    use itertools::Itertools;
    use mzalign::{AlignScoring, AlignType, Alignment, Side};
    use mzcore::{ontology::STATIC_ONTOLOGIES, prelude::Peptidoform};

    use crate::generate_annotations::generate_annotations;

    #[test]
    fn domain() {
        let ab = Peptidoform::pro_forma("QVTLRESGPVRVKPTLTETLTCAGSGFPLSDTGVRAGSGFSLGDPGVGVSWIRQPPGKALEWLAHIFSDDEKFYNASLKTRLTVSKDTSKGQVVLRLTNMDPVDTATYFCARVGRGYDSESGFHDKAMVWFDSWGKGTQVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDK", &STATIC_ONTOLOGIES).unwrap().0.into_simple_linear().unwrap();
        let ab1 = Peptidoform::pro_forma("GRGYDSESGFHDKAMVWFDSWGKGTQVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDK", &STATIC_ONTOLOGIES).unwrap().0.into_simple_linear().unwrap();
        let ab2 = Peptidoform::pro_forma("ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDK", &STATIC_ONTOLOGIES).unwrap().0.into_simple_linear().unwrap();
        let v = &imgt::STATIC_IMGT
            .data()
            .get(&imgt::Species::HomoSapiens)
            .unwrap()
            .find_allele(Gene::from_imgt_name("IGHV2-26*01").unwrap(), Some(1))
            .unwrap();
        let j = &imgt::STATIC_IMGT
            .data()
            .get(&imgt::Species::HomoSapiens)
            .unwrap()
            .find_allele(Gene::from_imgt_name("IGHJ5*01").unwrap(), Some(1))
            .unwrap();
        let c = &imgt::STATIC_IMGT
            .data()
            .get(&imgt::Species::HomoSapiens)
            .unwrap()
            .find_allele(Gene::from_imgt_name("IGHG1*01").unwrap(), Some(1))
            .unwrap();
        let alignments = vec![
            (
                v.clone(),
                Alignment::create_from_path(
                    v,
                    &ab,
                    0,
                    0,
                    "4=1X5=1X4=3r4=9I1=1I1=3I5=5X21=1X3=1X1=2X3=1X3=1X6=1X4=2X11=1X3=1X",
                    AlignScoring::default(),
                    AlignType {
                        left: Side::Specified { a: true, b: true },
                        right: Side::EitherGlobal,
                    },
                    4,
                )
                .unwrap(),
            ),
            (
                j.clone(),
                Alignment::create_from_path(
                    j,
                    &ab1,
                    0,
                    11,
                    "1=1X3I6=1X2=1X5=",
                    AlignScoring::default(),
                    AlignType {
                        left: Side::Specified { a: true, b: false },
                        right: Side::EitherGlobal,
                    },
                    4,
                )
                .unwrap(),
            ),
            (
                c.clone(),
                Alignment::create_from_path(
                    c,
                    &ab2,
                    0,
                    0,
                    "105=",
                    AlignScoring::default(),
                    AlignType {
                        left: Side::Specified { a: true, b: true },
                        right: Side::EitherGlobal,
                    },
                    4,
                )
                .unwrap(),
            ),
        ];
        let (regions, annotations) = generate_annotations(&alignments);
        assert_eq!(
            regions.iter().map(|(r, l)| format!("{r}:{l}")).join(";"),
            "FR1:38;CDR1:10;FR2:17;CDR2:7;FR3:38;CDR3:23;FR4:11;CH1:98;H:7"
        );
        assert_eq!(
            annotations
                .iter()
                .map(|(r, l)| format!("{r}:{l}"))
                .join(";"),
            "Conserved:21;Conserved:50;Conserved:109;Conserved:133;Conserved:134;Conserved:136"
        );
    }
}
