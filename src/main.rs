use clap::Parser;
use colored::{Colorize, Styles};
use imgt_germlines::{Allele, GeneType, Selection};
use rayon::prelude::*;
use rustyms::{
    align::{Alignment, MatchType, Piece, Type},
    find_isobaric_sets,
    modification::GnoComposition,
    ontologies::*,
    placement_rule::*,
    AminoAcid, Chemical, ComplexPeptide, LinearPeptide, MassTolerance, Modification,
};
use std::{collections::HashSet, io::Write, process::exit};

mod cli;
mod legend;
mod render;
mod styling;

use cli::*;
use render::*;
use styling::*;

fn main() {
    let args = Cli::parse();
    if let (Some(x), Some(y)) = (&args.x, &args.second.y) {
        let a = parse_peptide(x);
        let b = parse_peptide(y);
        let alignment = align(
            a.clone(),
            b.clone(),
            args.tolerance,
            args.alignment_type.ty(),
            args.normal,
        );
        show_annotated_mass_alignment(
            &alignment,
            args.tolerance,
            None,
            args.line_width,
            args.context,
            false,
        );
    } else if let (Some(x), Some(path)) = (&args.x, &args.second.file) {
        let sequences = rustyms::identifications::FastaData::parse_file(path).unwrap();
        let search_sequence = parse_peptide(x);
        let mut alignments: Vec<_> = sequences
            .into_par_iter()
            .map(|seq| {
                let a = seq.peptide.clone();
                (
                    seq,
                    align(
                        a,
                        search_sequence.clone(),
                        args.tolerance,
                        args.alignment_type.ty(),
                        args.normal,
                    ),
                )
            })
            .collect();
        alignments.sort_unstable_by(|a, b| b.1.normalised_score.total_cmp(&a.1.normalised_score));
        let selected: Vec<_> = alignments.into_iter().take(args.number_of_hits).collect();
        let mut data = vec![[
            "Rank".to_string(),
            "Database id".to_string(),
            "Score".to_string(),
            "Normalised score".to_string(),
            "Identity".to_string(),
            "Similarity".to_string(),
            "Gap".to_string(),
        ]];
        for (rank, (fasta, alignment)) in selected.iter().enumerate() {
            let stats = alignment.stats();
            data.push([
                (rank + 1).to_string(),
                fasta.id.clone(),
                alignment.absolute_score.to_string(),
                format!("{:.3}", alignment.normalised_score),
                format!("{:.2}%", stats.0 as f64 / stats.3 as f64 * 100.0),
                format!("{:.2}%", stats.1 as f64 / stats.3 as f64 * 100.0),
                format!("{:.2}%", stats.2 as f64 / stats.3 as f64 * 100.0),
            ]);
        }
        table(
            &data,
            true,
            &[
                Styling::with_style(Styles::Dimmed),
                Styling::none(),
                Styling::none(),
                Styling::none(),
                Styling::none(),
                Styling::none(),
                Styling::none(),
            ],
        );
        println!(
            "{} ({})",
            "Alignment for the best match".underline().italic(),
            selected[0].0.id.dimmed()
        );
        show_annotated_mass_alignment(
            &selected[0].1,
            args.tolerance,
            None,
            args.line_width,
            args.context,
            false,
        );
    } else if let (Some(x), Some(IMGTSelection::Search(selection))) = (&args.x, &args.second.imgt) {
        let seq_b = parse_peptide(x);
        let mut alignments: Vec<_> = selection
            .par_germlines()
            .map(|seq| {
                let a = seq.sequence.clone();
                (
                    seq,
                    align(
                        a,
                        seq_b.clone(),
                        args.tolerance,
                        args.alignment_type.ty(),
                        args.normal,
                    ),
                )
            })
            .collect();
        alignments.sort_unstable_by(|a, b| b.1.normalised_score.total_cmp(&a.1.normalised_score));
        let selected: Vec<_> = alignments.into_iter().take(args.number_of_hits).collect();
        let mut data = vec![[
            "Rank".to_string(),
            "Species".to_string(),
            "IMGT Name".to_string(),
            "Score".to_string(),
            "Normalised score".to_string(),
            "Identity".to_string(),
            "Similarity".to_string(),
            "Gap".to_string(),
        ]];
        for (rank, (imgt, alignment)) in selected.iter().enumerate() {
            let stats = alignment.stats();
            data.push([
                (rank + 1).to_string(),
                imgt.species.scientific_name().to_string(),
                format!("{} / {}", imgt.name(), imgt.fancy_name()),
                alignment.absolute_score.to_string(),
                format!("{:.3}", alignment.normalised_score),
                format!("{:.3}", stats.0 as f64 / stats.3 as f64),
                format!("{:.3}", stats.1 as f64 / stats.3 as f64),
                format!("{:.3}", stats.2 as f64 / stats.3 as f64),
            ]);
        }
        table(
            &data,
            true,
            &[
                Styling::with_style(Styles::Dimmed),
                Styling::none(),
                Styling::none(),
                Styling::none(),
                Styling::none(),
                Styling::none(),
                Styling::none(),
                Styling::none(),
            ],
        );
        println!(
            "{} ({} {})",
            "Alignment for the best match".underline().italic(),
            selected[0].0.species.scientific_name().dimmed(),
            format!("{} / {}", selected[0].0.name(), selected[0].0.fancy_name()).dimmed(),
        );
        show_annotated_mass_alignment(
            &selected[0].1,
            args.tolerance,
            Some(&selected[0].0),
            args.line_width,
            args.context,
            false,
        );
    } else if let (Some(x), Some(IMGTSelection::Domain(species, chains, allele))) =
        (&args.x, &args.second.imgt)
    {
        let seq_b = parse_peptide(x);
        let v_selection = Selection {
            species: species.clone(),
            chains: chains.clone(),
            genes: Some(HashSet::from([GeneType::V])),
            allele: allele.clone(),
        };

        let mut v_alignments: Vec<_> = v_selection
            .par_germlines()
            .map(|seq| {
                let a = seq.sequence.clone();
                (
                    seq,
                    align(
                        a,
                        seq_b.clone(),
                        args.tolerance,
                        Type::GlobalForA,
                        args.normal,
                    ),
                )
            })
            .collect();
        v_alignments.sort_unstable_by(|a, b| b.1.normalised_score.total_cmp(&a.1.normalised_score));
        let selected: Vec<_> = v_alignments.into_iter().take(args.number_of_hits).collect();
        let mut data = vec![[
            "Rank".to_string(),
            "Species".to_string(),
            "IMGT Name".to_string(),
            "Score".to_string(),
            "Normalised score".to_string(),
            "Identity".to_string(),
            "Similarity".to_string(),
            "Gap".to_string(),
        ]];
        for (rank, (imgt, alignment)) in selected.iter().enumerate() {
            let stats = alignment.stats();
            data.push([
                (rank + 1).to_string(),
                imgt.species.scientific_name().to_string(),
                format!("{} / {}", imgt.name(), imgt.fancy_name()),
                alignment.absolute_score.to_string(),
                format!("{:.3}", alignment.normalised_score),
                format!("{:.3}", stats.0 as f64 / stats.3 as f64),
                format!("{:.3}", stats.1 as f64 / stats.3 as f64),
                format!("{:.3}", stats.2 as f64 / stats.3 as f64),
            ]);
        }
        table(
            &data,
            true,
            &[
                Styling::with_style(Styles::Dimmed),
                Styling::none(),
                Styling::none(),
                Styling::none(),
                Styling::none(),
                Styling::none(),
                Styling::none(),
                Styling::none(),
            ],
        );
        let j_selection = Selection {
            species: Some(HashSet::from([selected[0].0.species])),
            chains: chains.clone(),
            genes: Some(HashSet::from([GeneType::J])),
            allele: allele.clone(),
        };
        let v_length = selected[0].1.start_a + selected[0].1.len_a();
        let mut j_alignments: Vec<_> = j_selection
            .par_germlines()
            .map(|seq| {
                let a = seq.sequence.clone();
                (
                    seq,
                    align(
                        a,
                        seq_b.clone().sequence.into_iter().skip(v_length).into(),
                        args.tolerance,
                        Type::GlobalForA,
                        args.normal,
                    ),
                )
            })
            .collect();
        j_alignments.sort_unstable_by(|a, b| b.1.normalised_score.total_cmp(&a.1.normalised_score));
        let j_alignment = j_alignments.first().unwrap();
        println!(
            "Best scoring species: {}",
            selected[0].0.species.scientific_name().green(),
        );
        let v_alignment = Alignment {
            seq_b: selected[0]
                .1
                .seq_b
                .sequence
                .iter()
                .take(v_length + 1)
                .cloned()
                .into(),
            ..selected[0].1.clone()
        };
        show_chained_annotated_mass_alignment(
            &(&selected[0].0, v_alignment),
            j_alignment,
            args.tolerance,
            args.line_width,
            args.context,
        );
    } else if let (Some(x), Some(IMGTSelection::Gene(species, gene, allele))) =
        (&args.x, &args.second.imgt)
    {
        if let Some(allele) = imgt_germlines::get_germline(*species, gene.clone(), *allele) {
            let alignment = align(
                allele.sequence.clone(),
                parse_peptide(x),
                args.tolerance,
                args.alignment_type.ty(),
                args.normal,
            );
            println!(
                "Selected: {} {}",
                allele.species.scientific_name().to_string().purple(),
                format!("{} / {}", allele.name(), allele.fancy_name()).purple(),
            );
            show_annotated_mass_alignment(
                &alignment,
                args.tolerance,
                Some(&allele),
                args.line_width,
                args.context,
                false,
            );
        } else {
            println!("Could not find specified germline")
        }
    } else if let Some(x) = &args.x {
        single_stats(&args, parse_peptide(x))
    } else if let Some(modification) = &args.modification {
        modification_stats(modification, args.tolerance);
    } else if let Some(IMGTSelection::Gene(species, gene, allele)) = &args.second.imgt {
        if let Some(allele) = imgt_germlines::get_germline(*species, gene.clone(), *allele) {
            display_germline(allele, &args);
        } else {
            println!("Could not find specified germline")
        }
    } else if let Some(IMGTSelection::Search(selection)) = &args.second.imgt {
        let mut first = true;
        for allele in selection.germlines() {
            if !first {
                println!();
            } else {
                first = false;
            }
            display_germline(allele, &args);
        }
    } else {
        println!("Please provide an argument to work with, use --help to see all options.")
    }
}

/// Parses the peptide and assumes it to be linear or shows an error and exit
fn parse_peptide(input: &str) -> LinearPeptide {
    match ComplexPeptide::pro_forma(input) {
        Ok(v) => v.assume_linear(),
        Err(e) => {
            println!("Peptide is not valid ProForma: {}", e);
            exit(1);
        }
    }
}

fn single_stats(args: &Cli, seq: LinearPeptide) {
    if let Some(complete) = seq.formula().and_then(|f| f.monoisotopic_mass()) {
        let bare = seq.bare_formula().unwrap().monoisotopic_mass().unwrap();
        println!(
            "Full mass: {} Da | {} Da {}",
            format!("{:.4}", complete.value).yellow(),
            format!(
                "{:.4}",
                seq.formula().unwrap().average_weight().unwrap().value
            )
            .yellow(),
            "(monoisotopic | average)".dimmed(),
        );
        println!(
            "Bare mass: {} Da {}",
            format!("{:.4}", bare.value).yellow(),
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
                        args.amino_acids
                            .as_deref()
                            .unwrap_or(AminoAcid::UNIQUE_MASS_AMINO_ACIDS),
                        args.fixed.mods(),
                        args.variable.mods(),
                        args.include.as_ref(),
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
                        args.amino_acids
                            .as_deref()
                            .unwrap_or(AminoAcid::UNIQUE_MASS_AMINO_ACIDS),
                        args.fixed.mods(),
                        args.variable.mods(),
                        args.include.as_ref(),
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
            format!("{:.4}", monoisotopic.value).yellow(),
            format!(
                "{:.4}",
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

fn display_germline(allele: Allele, args: &Cli) {
    let alignment = Alignment {
        absolute_score: 0,
        normalised_score: 0.0,
        path: vec![Piece::new(0, 0, MatchType::FullIdentity, 1, 0); allele.sequence.len()],
        start_a: 0,
        start_b: 0,
        seq_a: allele.sequence.clone(),
        seq_b: LinearPeptide::default(),
        ty: rustyms::align::Type::Global,
    };
    println!(
        "{} ({}) {}",
        allele.species.scientific_name().to_string().purple(),
        allele.species.common_name(),
        format!("{} / {}", allele.name(), allele.fancy_name()).purple(),
    );
    show_annotated_mass_alignment(
        &alignment,
        args.tolerance,
        Some(&allele),
        args.line_width,
        args.context,
        true,
    );
}

fn align(
    seq_a: LinearPeptide,
    seq_b: LinearPeptide,
    tolerance: MassTolerance,
    ty: Type,
    normal: bool,
) -> Alignment {
    if normal {
        rustyms::align::align::<1>(seq_a, seq_b, rustyms::align::BLOSUM62, tolerance, ty)
    } else {
        rustyms::align::align::<4>(seq_a, seq_b, rustyms::align::BLOSUM62, tolerance, ty)
    }
}
