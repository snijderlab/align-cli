use clap::Parser;
use colored::{Color, Colorize, Styles};
use imgt::{Allele, GeneType};
use itertools::Itertools;
use rayon::prelude::*;
use rustyms::align::par_consecutive_align;
use rustyms::{
    align::*,
    find_isobaric_sets, imgt,
    modification::{GnoComposition, Ontology},
    placement_rule::*,
    AminoAcid, Chemical, ComplexPeptide, LinearPeptide, MassComparable, Modification,
    MolecularFormula, Multi, MultiChemical, Tolerance,
};
use std::{
    collections::HashSet,
    io::{BufWriter, Write},
    path::Path,
    process::exit,
};

mod cli;
mod legend;
mod render;
mod styling;

use cli::*;
use render::*;
use styling::*;

fn main() {
    let args = Cli::parse();
    if let (Some(a), Some(b)) = (&args.a, &args.second.b) {
        let a = LinearPeptide::pro_forma(a).unwrap();
        let b = LinearPeptide::pro_forma(b).unwrap();
        let alignment = align(
            &a,
            &b,
            args.tolerance,
            args.alignment_type.ty(),
            args.scoring_matrix.matrix(),
            args.alignment_kind,
        );
        show_annotated_mass_alignment(
            &alignment,
            args.tolerance,
            None,
            args.line_width,
            args.context,
            false,
            ("A", "B"),
        );
    } else if let (Some(b), Some(path)) = (&args.a, &args.second.file) {
        let sequences = rustyms::identification::FastaData::parse_file(path).unwrap();
        let search_sequence = LinearPeptide::pro_forma(b).unwrap();
        let mut alignments: Vec<_> = sequences
            .into_par_iter()
            .map(|seq| {
                let sequence = seq.peptide.sequence.iter().cloned().collect();
                let alignment = align(
                    &sequence,
                    &search_sequence,
                    args.tolerance,
                    args.alignment_type.ty(),
                    args.scoring_matrix.matrix(),
                    args.alignment_kind,
                );
                (seq, alignment.to_owned())
            })
            .collect();
        alignments.sort_unstable_by(|a, b| b.1.cmp(&a.1));
        let selected: Vec<_> = alignments.into_iter().take(args.number_of_hits).collect();
        let mut data = vec![[
            String::new(),
            "Id".to_string(),
            "Score".to_string(),
            "Normalised score".to_string(),
            "Identity".to_string(),
            "Mass similarity".to_string(),
            "Gap".to_string(),
        ]];
        for (rank, (fasta, alignment)) in selected.iter().enumerate() {
            let stats = alignment.stats();
            data.push([
                (rank + 1).to_string(),
                fasta.id.clone(),
                alignment.score().absolute.to_string(),
                format!("{:.3}", alignment.normalised_score()),
                format!("{:.2}%", stats.identity() * 100.0),
                format!("{:.2}%", stats.mass_similarity() * 100.0),
                format!("{:.2}%", stats.gaps_fraction() * 100.0),
            ]);
        }
        table(
            &data,
            true,
            &[
                Styling::with_style(Styles::Dimmed),
                Styling::none(),
                Styling::none(),
                Styling::with_fg(Some(Color::Green)),
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
            (&selected[0].0.id, "Query"),
        );
    } else if let (Some(x), Some(IMGTSelection::Search(selection))) = (&args.a, &args.second.imgt) {
        let seq_b = parse_peptide(x);
        let mut alignments: Vec<_> = selection
            .clone()
            .par_germlines()
            .map(|seq| {
                let alignment = align(
                    seq.sequence,
                    &seq_b,
                    args.tolerance,
                    args.alignment_type.ty(),
                    args.scoring_matrix.matrix(),
                    args.alignment_kind,
                );
                (seq, alignment)
            })
            .collect();
        alignments.sort_unstable_by(|a, b| b.1.cmp(&a.1));
        let selected: Vec<_> = alignments.into_iter().take(args.number_of_hits).collect();
        let mut data = vec![[
            String::new(),
            "Species".to_string(),
            "IMGT name".to_string(),
            "Alternative name".to_string(),
            "Score".to_string(),
            "Normalised score".to_string(),
            "Identity".to_string(),
            "Mass similarity".to_string(),
            "Gap".to_string(),
        ]];
        for (rank, (imgt, alignment)) in selected.iter().enumerate() {
            let stats = alignment.stats();
            data.push([
                (rank + 1).to_string(),
                imgt.species.scientific_name().to_string(),
                imgt.name(),
                imgt.fancy_name(),
                alignment.score().absolute.to_string(),
                format!("{:.3}", alignment.normalised_score()),
                format!("{:.2}%", stats.identity() * 100.0),
                format!("{:.2}%", stats.mass_similarity() * 100.0),
                format!("{:.2}%", stats.gaps_fraction() * 100.0),
            ]);
        }
        table(
            &data,
            true,
            &[
                Styling::with_style(Styles::Dimmed),
                Styling::none(),
                Styling::none(),
                Styling::with_style(Styles::Dimmed),
                Styling::none(),
                Styling::with_fg(Some(Color::Green)),
                Styling::none(),
                Styling::none(),
                Styling::none(),
            ],
        );
        println!(
            "{} ({} {} {})",
            "Alignment for the best match".underline().italic(),
            selected[0].0.species.scientific_name().dimmed(),
            selected[0].0.species.common_name().dimmed(),
            format!("{} / {}", selected[0].0.name(), selected[0].0.fancy_name()).dimmed(),
        );
        show_annotated_mass_alignment(
            &selected[0].1,
            args.tolerance,
            Some(&selected[0].0),
            args.line_width,
            args.context,
            false,
            (selected[0].0.name(), "Query"),
        );
    } else if let (Some(x), Some(IMGTSelection::Domain(species, chains, allele))) =
        (&args.a, &args.second.imgt)
    {
        let scores = consecutive_align(
            &LinearPeptide::pro_forma(x).unwrap(),
            species.clone(),
            chains.clone(),
            allele.clone(),
            args.tolerance,
            args.scoring_matrix.matrix(),
            args.number_of_hits,
            args.alignment_kind,
        );

        for gene in &scores {
            let mut data = vec![[
                String::new(),
                "Species".to_string(),
                "IMGT name".to_string(),
                "Alternative name".to_string(),
                "Score".to_string(),
                "Normalised score".to_string(),
                "Identity".to_string(),
                "Mass similarity".to_string(),
                "Gap".to_string(),
            ]];
            for (rank, (imgt, alignment)) in gene.iter().enumerate() {
                let stats = alignment.stats();
                data.push([
                    (rank + 1).to_string(),
                    imgt.species.scientific_name().to_string(),
                    imgt.name(),
                    imgt.fancy_name(),
                    alignment.score().absolute.to_string(),
                    format!("{:.3}", alignment.normalised_score()),
                    format!("{:.2}%", stats.identity() * 100.0),
                    format!("{:.2}%", stats.mass_similarity() * 100.0),
                    format!("{:.2}%", stats.gaps_fraction() * 100.0),
                ]);
            }
            table(
                &data,
                true,
                &[
                    Styling::with_style(Styles::Dimmed),
                    Styling::none(),
                    Styling::none(),
                    Styling::with_style(Styles::Dimmed),
                    Styling::none(),
                    Styling::with_fg(Some(Color::Green)),
                    Styling::none(),
                    Styling::none(),
                    Styling::none(),
                ],
            );
        }

        let tops = scores
            .into_iter()
            .map(|options| options[0].clone())
            .collect_vec();
        show_chained_annotated_mass_alignment(&tops, args.tolerance, args.line_width, args.context);
    } else if let (Some(x), Some(IMGTSelection::Gene(species, gene, allele))) =
        (&args.a, &args.second.imgt)
    {
        if let Some(allele) = imgt::get_germline(*species, gene.clone(), *allele) {
            let b = LinearPeptide::pro_forma(x).unwrap();
            let alignment = align(
                allele.sequence,
                &b,
                args.tolerance,
                args.alignment_type.ty(),
                args.scoring_matrix.matrix(),
                args.alignment_kind,
            );
            println!(
                "Selected: {} {} {}",
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
                false,
                (allele.name(), "Query"),
            );
        } else {
            println!("Could not find specified germline")
        }
    } else if let Some(x) = &args.a {
        single_stats(&args, parse_peptide(x))
    } else if let Some(modification) = &args.modification {
        modification_stats(modification, args.tolerance);
    } else if let Some(file) = &args.second.csv {
        let csv = rustyms::csv::parse_csv(file, b',', None).unwrap();
        let output = std::fs::File::create(
            Path::new(file).with_file_name(
                Path::new(file)
                    .file_name()
                    .unwrap_or_default()
                    .to_str()
                    .unwrap_or_default()
                    .to_owned()
                    + "_output.csv",
            ),
        )
        .unwrap();
        let mut writer = BufWriter::new(output);
        let mut first = true;
        for line in csv {
            let line = line.unwrap();
            if first {
                writeln!(
                    writer,
                    "{},path,score,absolute score,maximal score,identical,mass similar,gaps,length",
                    line.headers().join(",")
                )
                .unwrap();
                first = false;
            }
            let a = LinearPeptide::pro_forma(line.index_column("a").unwrap().0).unwrap();
            let b = LinearPeptide::pro_forma(line.index_column("b").unwrap().0).unwrap();
            let alignment = align(
                &a,
                &b,
                args.tolerance,
                args.alignment_type.ty(),
                args.scoring_matrix.matrix(),
                args.alignment_kind,
            );
            let stats = alignment.stats();
            let score = alignment.score();
            writeln!(
                writer,
                "{},{},{},{},{},{},{},{},{}",
                line.line(),
                alignment.short(),
                score.normalised,
                score.absolute,
                score.max,
                stats.identical,
                stats.mass_similar,
                stats.gaps,
                stats.length
            )
            .unwrap();
        }
    } else if let Some(IMGTSelection::Gene(species, gene, allele)) = &args.second.imgt {
        if let Some(allele) = imgt::get_germline(*species, gene.clone(), *allele) {
            display_germline(allele, &args);
        } else {
            println!("Could not find specified germline")
        }
    } else if let Some(IMGTSelection::Search(selection)) = &args.second.imgt {
        let mut first = true;
        for allele in selection.clone().germlines() {
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
        Ok(v) => v
            .singular()
            .expect("Expected a singular peptide but a chimeric peptide was supplied"),
        Err(e) => {
            println!("Peptide is not valid ProForma: {}", e);
            exit(1);
        }
    }
}

fn single_stats(args: &Cli, seq: LinearPeptide) {
    let full_formulas = seq.formulas().unique();
    let bare_formulas = seq.bare_formulas().unique();
    print_multi_formula(&full_formulas, "Full", "");
    print_multi_formula(&bare_formulas, "Bare", "no N/C terminal taken into account");
    let multiple = full_formulas.len() > 1;

    let bare = seq
        .bare_formulas()
        .mass_bounds()
        .into_option()
        .expect("No masses for peptide")
        .0
        .monoisotopic_mass();
    println!();
    if multiple {
        println!("{}", "Multiple precursor masses found, it will generate isobaric options based on the lowest bare mass".dimmed().italic());
    }
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
}

fn print_multi_formula(formulas: &Multi<MolecularFormula>, prefix: &str, suffix: &str) {
    let multiple = formulas.len() > 1;
    print!("{}: ", prefix);
    if multiple {
        println!(
            "{}",
            if suffix.is_empty() {
                String::new()
            } else {
                format!("({suffix})").dimmed().to_string()
            }
        )
    }
    let mut lengths = (0, 0, 0);
    let mut rows = Vec::with_capacity(formulas.len());
    for formula in formulas.iter() {
        let row = (
            formula.hill_notation_fancy().green(),
            display_mass(formula.monoisotopic_mass(), true),
            display_mass(formula.average_weight(), true),
        );
        lengths = (
            lengths.0.max(row.0.chars().count()),
            lengths.1.max(row.1.chars().count()),
            lengths.2.max(row.2.chars().count()),
        );
        rows.push(row);
    }
    for formula in formulas.iter() {
        if multiple {
            print!("  ");
        }
        print!(
            "{:3$} {:4$} {:5$}",
            formula.hill_notation_fancy().green(),
            display_mass(formula.monoisotopic_mass(), true),
            display_mass(formula.average_weight(), true),
            lengths.0,
            lengths.1,
            lengths.2,
        );
        if multiple {
            println!();
        }
    }
    if multiple {
        print!("  ");
    } else {
        print!(" ");
    }
    println!(
        "{}{}",
        "(formula | monoisotopic mass | average weight)".dimmed(),
        if suffix.is_empty() || multiple {
            String::new()
        } else {
            format!(" ({suffix})").dimmed().to_string()
        }
    )
}

fn modification_stats(modification: &Modification, tolerance: Tolerance) {
    if let Modification::Mass(mass) = modification {
        println!(
            "All ontology modifications close to the given monoisotopic mass: {}",
            format!("tolerance: {tolerance}").dimmed()
        );
        let mut data = vec![[
            "Name".to_string(),
            "Id".to_string(),
            "Monoisotopic mass".to_string(),
            "Formula".to_string(),
        ]];
        for ontology in &[Ontology::Unimod, Ontology::Psimod, Ontology::Gnome] {
            for (id, _, modification) in ontology.lookup() {
                if tolerance.within(
                    &mass.into_inner(),
                    &modification.formula().monoisotopic_mass(),
                ) {
                    data.push([
                        modification.to_string(),
                        format!("{}:{}", ontology.name(), id),
                        display_mass(modification.formula().monoisotopic_mass(), false),
                        modification.formula().hill_notation_fancy(),
                    ])
                }
            }
        }
        if data.len() > 1 {
            table(
                &data,
                true,
                &[
                    Styling::with_fg(Some(Color::Magenta)),
                    Styling::with_style(Styles::Dimmed),
                    Styling::with_fg(Some(Color::Yellow)),
                    Styling::with_fg(Some(Color::Green)),
                ],
            );
        } else {
            println!("{}", "No modifications found".red())
        }
    } else if let Modification::Formula(formula) = modification {
        println!(
            "Full mass: {} {} {}\n",
            display_mass(formula.monoisotopic_mass(), true),
            display_mass(formula.average_weight(), true),
            "(monoisotopic | average)".dimmed(),
        );

        println!("All ontology modifications with the same formula:");
        let mut data = vec![["Name".to_string(), "Id".to_string()]];
        for ontology in &[Ontology::Unimod, Ontology::Psimod, Ontology::Gnome] {
            for (id, _, modification) in ontology.lookup() {
                if modification.formula() == *formula {
                    data.push([
                        modification.to_string(),
                        format!("{}:{}", ontology.name(), id),
                    ])
                }
            }
        }
        if data.len() > 1 {
            table(
                &data,
                true,
                &[
                    Styling::with_fg(Some(Color::Magenta)),
                    Styling::with_style(Styles::Dimmed),
                ],
            );
        } else {
            println!("{}", "No modifications found".red())
        }
    } else {
        let monoisotopic = modification.formula().monoisotopic_mass();
        println!(
            "Full mass: {} {} {}",
            display_mass(monoisotopic, true),
            display_mass(modification.formula().average_weight(), true),
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
                            aa.iter().map(|a| a.char()).collect::<String>().yellow(),
                            pos.to_string().green()
                        )
                    }
                    PlacementRule::PsiModification(index, pos) => {
                        print!(
                            "{}{}@{}",
                            if first { "" } else { ", " },
                            Ontology::Psimod.find_id(*index).unwrap().to_string().blue(),
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
    }
}

fn display_germline(allele: Allele, args: &Cli) {
    let alignment = rustyms::align::align::<1>(
        allele.sequence,
        allele.sequence,
        rustyms::align::BLOSUM90,
        Tolerance::new_ppm(10.0),
        rustyms::align::AlignType::GLOBAL,
    );
    println!(
        "{} {} {}",
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
        ("", ""),
    );
}

fn align<'a>(
    seq_a: &'a LinearPeptide,
    seq_b: &'a LinearPeptide,
    tolerance: Tolerance,
    ty: AlignType,
    matrix: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
    kind: AlignmentKind,
) -> RefAlignment<'a> {
    if kind.normal {
        rustyms::align::align::<1>(seq_a, seq_b, matrix, tolerance, ty)
    } else if kind.mass_based_huge {
        rustyms::align::align::<{ u16::MAX }>(seq_a, seq_b, matrix, tolerance, ty)
    } else if kind.mass_based_long {
        rustyms::align::align::<8>(seq_a, seq_b, matrix, tolerance, ty)
    } else {
        rustyms::align::align::<4>(seq_a, seq_b, matrix, tolerance, ty)
    }
}

fn consecutive_align(
    seq: &LinearPeptide,
    species: Option<HashSet<imgt::Species>>,
    chains: Option<HashSet<imgt::ChainType>>,
    allele: imgt::AlleleSelection,
    tolerance: Tolerance,
    matrix: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
    return_number: usize,
    kind: AlignmentKind,
) -> Vec<Vec<(Allele<'static>, OwnedAlignment)>> {
    if kind.normal {
        par_consecutive_align::<1>(
            seq,
            &[
                (GeneType::V, AlignType::new(Some((true, true)), None)),
                (GeneType::J, AlignType::new(Some((true, false)), None)),
                (GeneType::C(None), AlignType::new(Some((true, true)), None)),
            ],
            species.clone(),
            chains.clone(),
            allele.clone(),
            tolerance,
            matrix,
            return_number,
        )
    } else if kind.mass_based_huge {
        par_consecutive_align::<{ u16::MAX }>(
            seq,
            &[
                (GeneType::V, AlignType::new(Some((true, true)), None)),
                (GeneType::J, AlignType::new(Some((true, false)), None)),
                (GeneType::C(None), AlignType::new(Some((true, true)), None)),
            ],
            species.clone(),
            chains.clone(),
            allele.clone(),
            tolerance,
            matrix,
            return_number,
        )
    } else if kind.mass_based_long {
        par_consecutive_align::<8>(
            seq,
            &[
                (GeneType::V, AlignType::new(Some((true, true)), None)),
                (GeneType::J, AlignType::new(Some((true, false)), None)),
                (GeneType::C(None), AlignType::new(Some((true, true)), None)),
            ],
            species.clone(),
            chains.clone(),
            allele.clone(),
            tolerance,
            matrix,
            return_number,
        )
    } else {
        par_consecutive_align::<4>(
            seq,
            &[
                (GeneType::V, AlignType::new(Some((true, true)), None)),
                (GeneType::J, AlignType::new(Some((true, false)), None)),
                (GeneType::C(None), AlignType::new(Some((true, true)), None)),
            ],
            species.clone(),
            chains.clone(),
            allele.clone(),
            tolerance,
            matrix,
            return_number,
        )
    }
}
