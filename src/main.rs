use clap::Parser;
use colored::{Color, Colorize, Styles};
use imgt::{Allele, GeneType};
use itertools::Itertools;
use rayon::prelude::*;
use rustyms::align::par_consecutive_align;
use rustyms::imgt::Selection;
use rustyms::modification::{
    LinkerSpecificity, ModificationId, ModificationSearchResult, SimpleModification,
};
use rustyms::system::Mass;
use rustyms::{
    align::*,
    find_isobaric_sets, imgt,
    modification::{GnoComposition, Ontology},
    placement_rule::*,
    AminoAcid, Chemical, LinearPeptide, MolecularFormula, Multi, Tolerance,
};
use rustyms::{ExtremelySimple, Linear, Simple};
use std::{
    collections::HashSet,
    io::{BufWriter, Write},
    path::Path,
};

/// Define the default precision (in number of digits shown) for number output
const NUMBER_PRECISION: usize = 3;

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
        let a = LinearPeptide::pro_forma(a, None).unwrap().simple().unwrap();
        let b = LinearPeptide::pro_forma(b, None).unwrap().simple().unwrap();
        let alignment = align(
            &a,
            &b,
            args.tolerance,
            args.alignment_type.ty(),
            args.scoring_matrix.matrix(),
            args.alignment_kind,
        );
        show_annotated_mass_alignment(&alignment, None, false, ("A", "B"), &args);
    } else if let (Some(b), Some(path)) = (&args.a, &args.second.file) {
        let sequences = rustyms::identification::FastaData::parse_file(path).unwrap();
        let search_sequence = LinearPeptide::pro_forma(b, None).unwrap().simple().unwrap();
        let mut alignments: Vec<_> = sequences
            .into_par_iter()
            .map(|seq| {
                let sequence: LinearPeptide<Simple> =
                    seq.peptide.sequence.iter().cloned().collect();
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
            None,
            false,
            (&selected[0].0.id, "Query"),
            &args,
        );
    } else if let (Some(x), true) = (&args.a, &args.second.imgt) {
        let seq_b = LinearPeptide::pro_forma(x, None).unwrap().simple().unwrap();
        let mut alignments: Vec<_> = Selection {
            species: args.species.map(|s| HashSet::from([s])),
            chains: args.chains.clone(),
            genes: args.genes.clone(),
            allele: args.allele.clone(),
        }
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
        alignments
            .sort_unstable_by(|a, b| b.1.score().normalised.total_cmp(&a.1.score().normalised));
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
            Some(&selected[0].0),
            false,
            (selected[0].0.name(), "Query"),
            &args,
        );
    } else if let (Some(x), true) = (&args.a, &args.second.domain) {
        let scores = consecutive_align(
            &LinearPeptide::pro_forma(x, None).unwrap().simple().unwrap(),
            args.species.map(|s| HashSet::from([s])),
            args.chains.clone(),
            args.allele.clone(),
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
        show_chained_annotated_mass_alignment(
            &tops,
            args.tolerance,
            args.line_width,
            args.context,
            args.full_number,
        );
    } else if let (Some(x), Some((gene, allele)), Some(species)) =
        (&args.a, &args.second.specific_gene, &args.species)
    {
        if let Some(allele) = imgt::get_germline(*species, gene.clone(), *allele) {
            let b = LinearPeptide::pro_forma(x, None).unwrap().simple().unwrap();
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
                Some(&allele),
                false,
                (allele.name(), "Query"),
                &args,
            );
        } else {
            println!("Could not find specified germline")
        }
    } else if let Some(x) = &args.a {
        single_stats(
            &args,
            LinearPeptide::pro_forma(x, None).unwrap().simple().unwrap(),
        )
    } else if let Some(modification) = &args.modification {
        modification_stats(modification, args.tolerance, args.full_number);
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
            let a = LinearPeptide::pro_forma(line.index_column("a").unwrap().0, None)
                .unwrap()
                .simple()
                .unwrap();
            let b = LinearPeptide::pro_forma(line.index_column("b").unwrap().0, None)
                .unwrap()
                .simple()
                .unwrap();
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
    } else if let (Some((gene, allele)), Some(species)) =
        (&args.second.specific_gene, &args.species)
    {
        if let Some(allele) = imgt::get_germline(*species, gene.clone(), *allele) {
            display_germline(allele, &args);
        } else {
            println!("Could not find specified germline")
        }
    } else if args.second.imgt {
        let mut first = true;
        let selection = Selection {
            species: args.species.map(|s| HashSet::from([s])),
            chains: args.chains.clone(),
            genes: args.genes.clone(),
            allele: args.allele.clone(),
        };
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

fn single_stats(args: &Cli, seq: LinearPeptide<Simple>) {
    let full_formulas = seq.formulas().unique();
    let bare_formulas = seq.bare_formulas().unique();
    print_multi_formula(&full_formulas, "Full", "", args.full_number);
    print_multi_formula(
        &bare_formulas,
        "Bare",
        "no N/C terminal taken into account",
        args.full_number,
    );
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

fn print_multi_formula(
    formulas: &Multi<MolecularFormula>,
    prefix: &str,
    suffix: &str,
    full_number: bool,
) {
    let precision = if full_number {
        None
    } else {
        Some(NUMBER_PRECISION)
    };
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
    let mut lengths = (0, 0, 0, 0);
    let mut rows = Vec::with_capacity(formulas.len());
    for formula in formulas.iter() {
        let row = (
            formula.hill_notation_fancy().green(),
            display_mass(formula.monoisotopic_mass(), true, precision),
            display_mass(formula.average_weight(), true, precision),
            display_mass(formula.most_abundant_mass(), true, precision),
        );
        lengths = (
            lengths.0.max(row.0.chars().count()),
            lengths.1.max(row.1.chars().count()),
            lengths.2.max(row.2.chars().count()),
            lengths.3.max(row.3.chars().count()),
        );
        rows.push(row);
    }
    for formula in formulas.iter() {
        if multiple {
            print!("  ");
        }
        print!(
            "{:4$} {:5$} {:6$} {:7$}",
            formula.hill_notation_fancy().green(),
            display_mass(formula.monoisotopic_mass(), true, precision),
            display_mass(formula.average_weight(), true, precision),
            display_mass(formula.most_abundant_mass(), true, precision),
            lengths.0,
            lengths.1,
            lengths.2,
            lengths.3,
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
        "(formula | monoisotopic mass | average weight | most abundant mass)".dimmed(),
        if suffix.is_empty() || multiple {
            String::new()
        } else {
            format!(" ({suffix})").dimmed().to_string()
        }
    )
}

fn modification_stats(
    modification: &SimpleModification,
    tolerance: Tolerance<Mass>,
    full_number: bool,
) {
    let precision = if full_number {
        None
    } else {
        Some(NUMBER_PRECISION)
    };
    match SimpleModification::search(modification, tolerance, None) {
        ModificationSearchResult::Mass(_mass, tolerance, modifications) => {
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
            for (ontology, id, _name, modification) in modifications {
                data.push([
                    modification.to_string(),
                    format!("{}:{}", ontology.name(), id),
                    display_mass(modification.formula().monoisotopic_mass(), false, precision),
                    modification.formula().hill_notation_fancy(),
                ])
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
        }
        ModificationSearchResult::Formula(formula, modifications) => {
            println!(
                "Full mass: {} {} {}\n",
                display_mass(formula.monoisotopic_mass(), true, precision),
                display_mass(formula.average_weight(), true, precision),
                "(monoisotopic | average)".dimmed(),
            );

            println!("All ontology modifications with the same formula:");
            let mut data = vec![["Name".to_string(), "Id".to_string()]];
            for (ontology, id, _name, modification) in modifications {
                data.push([
                    modification.to_string(),
                    format!("{}:{}", ontology.name(), id),
                ])
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
        }
        ModificationSearchResult::Glycan(_composition, modifications) => {
            println!("All GNOme modifications with the same monosaccharide composition:");
            let mut data = vec![["Name".to_string(), "Structure".to_string()]];
            for (_ontology, _id, _name, modification) in modifications {
                if let SimpleModification::Gno(GnoComposition::Structure(structure), _) =
                    &modification
                {
                    data.push([modification.to_string(), structure.to_string()])
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
        }
        ModificationSearchResult::Single(modification) => {
            let monoisotopic = modification.formula().monoisotopic_mass();
            println!(
                "Full mass: {} {} {} {}",
                display_mass(monoisotopic, true, precision),
                display_mass(modification.formula().average_weight(), true, precision),
                display_mass(modification.formula().most_abundant_mass(), true, precision),
                "(monoisotopic | average | most abundant)".dimmed(),
            );
            if !modification.formula().is_empty() {
                println!(
                    "Composition: {}",
                    modification.formula().hill_notation_fancy().green(),
                );
            }
            match modification {
                SimpleModification::Database {
                    specificities, id, ..
                } => {
                    display_id(&id);
                    println!("Placement rules: ");

                    for rule in specificities {
                        print!("  Locations: ");
                        // Print locations
                        display_placement_rules(&rule.0);
                        // Print neutral losses
                        if !rule.1.is_empty() {
                            print!(
                                ", Neutral losses: {}",
                                rule.1
                                    .iter()
                                    .map(|n| n.hill_notation_fancy().yellow())
                                    .join(", ")
                            );
                        }
                        // Print diagnostic ions
                        if !rule.2.is_empty() {
                            print!(
                                ", Diagnostic ions: {}",
                                rule.2
                                    .iter()
                                    .map(|d| d.0.hill_notation_fancy().green())
                                    .join(", ")
                            );
                        }
                        println!();
                    }
                }
                SimpleModification::Linker {
                    specificities,
                    id,
                    length,
                    ..
                } => {
                    display_id(&id);
                    if let Some(length) = length {
                        println!("Length: {}", length);
                    }
                    println!("Placement rules: ");
                    for specificity in specificities {
                        match specificity {
                            LinkerSpecificity::Symmetric(locations, stubs, diagnostic) => {
                                print!("  Locations: ");
                                display_placement_rules(&locations);
                                if !stubs.is_empty() {
                                    print!(
                                        ", Cleave points: {}",
                                        stubs
                                            .iter()
                                            .map(|(a, b)| format!(
                                                "{} + {}",
                                                a.hill_notation_fancy().yellow(),
                                                b.hill_notation_fancy().yellow()
                                            ))
                                            .join(", ")
                                    );
                                }
                                if !diagnostic.is_empty() {
                                    print!(
                                        ", Diagnostic ions: {}",
                                        diagnostic
                                            .iter()
                                            .map(|d| d.0.hill_notation_fancy().green())
                                            .join(", ")
                                    );
                                }
                            }
                            LinkerSpecificity::Asymmetric(locations, stubs, diagnostic) => {
                                print!("  Left: ");
                                display_placement_rules(&locations.0);
                                print!(", Right: ");
                                display_placement_rules(&locations.1);

                                if !stubs.is_empty() {
                                    print!(
                                        ", Cleave points: {}",
                                        stubs
                                            .iter()
                                            .map(|(a, b)| format!(
                                                "{} + {}",
                                                a.hill_notation_fancy().yellow(),
                                                b.hill_notation_fancy().yellow()
                                            ))
                                            .join(", ")
                                    );
                                }
                                if !diagnostic.is_empty() {
                                    print!(
                                        ", Diagnostic ions: {}",
                                        diagnostic
                                            .iter()
                                            .map(|d| d.0.hill_notation_fancy().green())
                                            .join(", ")
                                    );
                                }
                            }
                        }
                    }
                }
                SimpleModification::Gno(composition, name) => {
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
                _ => (),
            }
        }
    }
}

fn display_placement_rules(rules: &[PlacementRule]) {
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
                    Ontology::Psimod
                        .find_id(*index, None)
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
            PlacementRule::Anywhere => print!("{}", "Anywhere".green()),
        }
        first = false;
    }
}

fn display_id(id: &ModificationId) {
    println!(
        "Ontology: {}, name: {}, index: {}",
        id.ontology.to_string().purple(),
        id.name.green(),
        id.id.to_string().blue()
    );
    if !id.description.is_empty() {
        println!("{}", id.description);
    }
    if !id.cross_ids.is_empty() {
        println!(
            "IDs: {}",
            id.cross_ids
                .iter()
                .map(|(r, i)| format!("{}{}{i}", r.dimmed(), ":".dimmed()))
                .join(", ")
        );
    }
    if !id.synonyms.is_empty() {
        println!("Synonyms: {}", id.synonyms.join(", "));
    }
}

fn display_germline(allele: Allele, args: &Cli) {
    let alignment = rustyms::align::align::<1, ExtremelySimple, ExtremelySimple>(
        allele.sequence,
        allele.sequence,
        rustyms::align::matrix::BLOSUM90,
        Tolerance::new_ppm(10.0),
        rustyms::align::AlignType::GLOBAL,
    );
    println!(
        "{} {} {}",
        allele.species.scientific_name().to_string().purple(),
        allele.species.common_name(),
        format!("{} / {}", allele.name(), allele.fancy_name()).purple(),
    );
    show_annotated_mass_alignment(&alignment, Some(&allele), true, ("", ""), &args);
}

fn align<'a, A: Into<Simple> + Into<Linear>, B: Into<Simple> + Into<Linear>>(
    seq_a: &'a LinearPeptide<A>,
    seq_b: &'a LinearPeptide<B>,
    tolerance: Tolerance<Mass>,
    ty: AlignType,
    matrix: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
    kind: AlignmentKind,
) -> RefAlignment<'a, A, B> {
    if kind.normal {
        rustyms::align::align::<1, A, B>(seq_a, seq_b, matrix, tolerance, ty)
    } else if kind.mass_based_huge {
        rustyms::align::align::<{ u16::MAX }, A, B>(seq_a, seq_b, matrix, tolerance, ty)
    } else if kind.mass_based_long {
        rustyms::align::align::<8, A, B>(seq_a, seq_b, matrix, tolerance, ty)
    } else {
        rustyms::align::align::<4, A, B>(seq_a, seq_b, matrix, tolerance, ty)
    }
}

fn consecutive_align(
    seq: &LinearPeptide<Simple>,
    species: Option<HashSet<imgt::Species>>,
    chains: Option<HashSet<imgt::ChainType>>,
    allele: imgt::AlleleSelection,
    tolerance: Tolerance<Mass>,
    matrix: &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER],
    return_number: usize,
    kind: AlignmentKind,
) -> Vec<Vec<(Allele<'static>, OwnedAlignment)>> {
    if kind.normal {
        par_consecutive_align::<1>(
            seq,
            &[
                (
                    GeneType::V,
                    AlignType {
                        left: Side::Specified { a: true, b: true },
                        right: Side::EitherGlobal,
                    },
                ),
                (
                    GeneType::J,
                    AlignType {
                        left: Side::Specified { a: true, b: false },
                        right: Side::EitherGlobal,
                    },
                ),
                (
                    GeneType::C(None),
                    AlignType {
                        left: Side::Specified { a: true, b: true },
                        right: Side::EitherGlobal,
                    },
                ),
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
                (
                    GeneType::V,
                    AlignType {
                        left: Side::Specified { a: true, b: true },
                        right: Side::EitherGlobal,
                    },
                ),
                (
                    GeneType::J,
                    AlignType {
                        left: Side::Specified { a: true, b: false },
                        right: Side::EitherGlobal,
                    },
                ),
                (
                    GeneType::C(None),
                    AlignType {
                        left: Side::Specified { a: true, b: true },
                        right: Side::EitherGlobal,
                    },
                ),
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
                (
                    GeneType::V,
                    AlignType {
                        left: Side::Specified { a: true, b: true },
                        right: Side::EitherGlobal,
                    },
                ),
                (
                    GeneType::J,
                    AlignType {
                        left: Side::Specified { a: true, b: false },
                        right: Side::EitherGlobal,
                    },
                ),
                (
                    GeneType::C(None),
                    AlignType {
                        left: Side::Specified { a: true, b: true },
                        right: Side::EitherGlobal,
                    },
                ),
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
                (
                    GeneType::V,
                    AlignType {
                        left: Side::Specified { a: true, b: true },
                        right: Side::EitherGlobal,
                    },
                ),
                (
                    GeneType::J,
                    AlignType {
                        left: Side::Specified { a: true, b: false },
                        right: Side::EitherGlobal,
                    },
                ),
                (
                    GeneType::C(None),
                    AlignType {
                        left: Side::Specified { a: true, b: true },
                        right: Side::EitherGlobal,
                    },
                ),
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
