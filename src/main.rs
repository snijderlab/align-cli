use clap::Parser;
use colored::{Color, Colorize, Styles};
use imgt::{Allele, GeneType};
use itertools::Itertools;
use rayon::prelude::*;
use rustyms::align::par_consecutive_align;
use rustyms::imgt::Selection;
use rustyms::system::{dalton, Mass};
use rustyms::{
    align::*,
    find_isobaric_sets, imgt,
    modification::{
        GnoComposition, LinkerSpecificity, ModificationId, Ontology, SimpleModification,
        SimpleModificationInner,
    },
    modification_search_formula, modification_search_glycan, modification_search_mass,
    placement_rule::*,
    AminoAcid, AtMax, Chemical, MassMode, MolecularFormula, Multi, Peptidoform, SimpleLinear,
    Tolerance, UnAmbiguous,
};
use rustyms::{find_formulas, Element};
use std::num::NonZeroU16;
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
        let a = Peptidoform::pro_forma(a, None)
            .unwrap()
            .into_simple_linear()
            .unwrap();
        let b = Peptidoform::pro_forma(b, None)
            .unwrap()
            .into_simple_linear()
            .unwrap();
        let alignment = align(
            &a,
            &b,
            args.scoring(),
            args.alignment_type.ty(),
            args.alignment_kind,
        );
        show_annotated_mass_alignment::<_, _, Allele>(
            &alignment,
            None,
            false,
            false,
            ("A", "B"),
            &args,
        );
    } else if let (Some(b), Some(path)) = (&args.a, &args.second.file) {
        let sequences = rustyms::identification::FastaData::parse_file(path).unwrap();
        let search_sequence = Peptidoform::pro_forma(b, None)
            .unwrap()
            .into_simple_linear()
            .unwrap();
        let mut alignments: Vec<_> = sequences
            .into_par_iter()
            .map(|seq| {
                let sequence = seq.peptide().clone();
                let alignment = align(
                    &sequence,
                    &search_sequence,
                    args.scoring(),
                    args.alignment_type.ty(),
                    args.alignment_kind,
                );
                (seq, alignment.to_owned())
            })
            .filter(|s| !s.1.normalised_score().is_nan())
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
                fasta.identifier().to_string(),
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
            selected[0].0.identifier().to_string().dimmed()
        );
        show_annotated_mass_alignment(
            &selected[0].1,
            Some(&selected[0].0),
            false,
            false,
            (&selected[0].0.identifier().to_string(), "Query"),
            &args,
        );
    } else if let (Some(x), true) = (&args.a, &args.second.imgt) {
        let seq_b = Peptidoform::pro_forma(x, None)
            .unwrap()
            .into_simple_linear()
            .unwrap();
        let mut alignments: Vec<_> = Selection {
            species: args.species.map(|s| HashSet::from([s])),
            chains: args.chains.clone(),
            genes: args.genes.clone(),
            allele: args.allele,
        }
        .par_germlines()
        .map(|seq| {
            let alignment = align(
                seq.sequence,
                &seq_b,
                args.scoring(),
                args.alignment_type.ty(),
                args.alignment_kind,
            );
            (seq, alignment)
        })
        .filter(|s| !s.1.normalised_score().is_nan())
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
            false,
            (selected[0].0.name(), "Query"),
            &args,
        );
    } else if let (Some(x), true) = (&args.a, &args.second.domain) {
        let scores = consecutive_align(
            &Peptidoform::pro_forma(x, None)
                .unwrap()
                .into_simple_linear()
                .unwrap(),
            args.species.map(|s| HashSet::from([s])),
            args.chains.clone(),
            args.allele,
            args.scoring(),
            args.number_of_hits,
            args.alignment_kind,
        );

        for gene in &scores.alignments {
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
            .alignments
            .into_iter()
            .map(|options| options[0].clone())
            .collect_vec();
        show_chained_annotated_mass_alignment(
            &tops,
            args.tolerance,
            args.line_width,
            args.context,
            args.full_number,
            args.generate_annotation,
        );
    } else if let (Some(x), Some((gene, allele)), Some(species)) =
        (&args.a, &args.second.specific_gene, &args.species)
    {
        if let Some(allele) = imgt::get_germline(*species, gene.clone(), *allele) {
            let b = Peptidoform::pro_forma(x, None)
                .unwrap()
                .into_simple_linear()
                .unwrap();
            let alignment = align(
                allele.sequence,
                &b,
                args.scoring(),
                args.alignment_type.ty(),
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
            Peptidoform::pro_forma(x, None)
                .unwrap()
                .into_simple_linear()
                .unwrap(),
        )
    } else if let Some(modification) = &args.modification {
        modification_stats(
            modification,
            args.tolerance,
            args.full_number,
            args.mass_mode,
            args.positions.as_deref(),
        );
    } else if let Some(file) = &args.second.csv {
        let csv = rustyms::csv::parse_csv(file, b',', None).expect("Failed to parse CSV file");
        let output = std::fs::File::create(
            Path::new(file).with_file_name(
                Path::new(file)
                    .file_name()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .to_string()
                    + "_output.csv",
            ),
        )
        .unwrap();
        let mut writer = BufWriter::new(output);
        let mut first = true;
        for line in csv {
            let line = line.expect("Failed to read CSV line");
            if first {
                writeln!(
                    writer,
                    "{},path,score,absolute score,maximal score,identical,mass similar,gaps,length",
                    line.headers().join(",")
                )
                .unwrap();
                first = false;
            }
            let a = Peptidoform::pro_forma(line.index_column("a").unwrap().0, None)
                .unwrap()
                .into_simple_linear()
                .unwrap();
            let b = Peptidoform::pro_forma(line.index_column("b").unwrap().0, None)
                .unwrap()
                .into_simple_linear()
                .unwrap();
            let alignment = align(
                &a,
                &b,
                args.scoring(),
                args.alignment_type.ty(),
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
            allele: args.allele,
        };
        for allele in selection.germlines() {
            if !first {
                println!();
            } else {
                first = false;
            }
            display_germline(allele, &args);
        }
    } else if let Some(target) = args.formula_target {
        const DEFAULT_ELEMENTS: &[(Element, Option<NonZeroU16>)] = &[
            (Element::H, None),
            (Element::C, None),
            (Element::O, None),
            (Element::N, None),
            (Element::S, None),
        ];
        let mut data = vec![["Formula".to_string(), "Mass".to_string()]];
        for formula in find_formulas(
            target.0,
            Tolerance::Absolute(Mass::new::<dalton>(10.0_f64.powf(-(target.1 as f64)) / 2.0)),
            DEFAULT_ELEMENTS,
        )
        .iter()
        {
            data.push([
                formula.hill_notation_fancy(),
                display_mass(
                    formula.monoisotopic_mass(),
                    false,
                    Some((target.1 + 2).max(3)),
                ),
            ]);
        }
        table(&data, true, &[Styling::none(), Styling::none()]);
    } else {
        println!("Please provide an argument to work with, use --help to see all options.")
    }
}

fn single_stats(args: &Cli, seq: Peptidoform<SimpleLinear>) {
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
        .mass(args.mass_mode);
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
    mass_mode: MassMode,
    positions: Option<&[(Vec<AminoAcid>, Position)]>,
) {
    let precision = if full_number {
        None
    } else {
        Some(NUMBER_PRECISION)
    };
    match &**modification {
        SimpleModificationInner::Mass(m)
        | SimpleModificationInner::Gno {
            composition: GnoComposition::Weight(m),
            ..
        } => {
            println!(
                "All ontology modifications close to the given {mass_mode}: {}",
                format!("tolerance: {tolerance}").dimmed()
            );
            let mut data = vec![[
                "Name".to_string(),
                "Id".to_string(),
                mass_mode.to_string(),
                "Formula".to_string(),
            ]];
            for (ontology, id, _name, modification) in
                modification_search_mass(m.into_inner(), tolerance, positions, mass_mode, None)
            {
                data.push([
                    modification.to_string(),
                    format!(
                        "{}{}",
                        ontology.name(),
                        id.map_or(String::new(), |id| format!(":{id}")),
                    ),
                    display_mass(modification.formula().mass(mass_mode), false, precision),
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
        SimpleModificationInner::Formula(f) => {
            display_single_mod(modification, precision);

            println!("\nAll ontology modifications with the same formula:");
            let mut data = vec![["Name".to_string(), "Id".to_string()]];
            for (ontology, id, _name, modification) in modification_search_formula(f, None) {
                data.push([
                    modification.to_string(),
                    format!(
                        "{}{}",
                        id.map_or(String::new(), |id| format!("{id}:")),
                        ontology.name()
                    ),
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
        SimpleModificationInner::Glycan(ref g)
        | SimpleModificationInner::Gno {
            composition: GnoComposition::Composition(ref g),
            ..
        } => {
            display_single_mod(modification, precision);

            println!("\nAll GNOme modifications with the same monosaccharide composition:");
            let mut data = vec![["Name".to_string(), "Definition".to_string()]];
            for (_ontology, _id, _name, modification) in modification_search_glycan(g, true) {
                if let SimpleModificationInner::Gno {
                    composition: GnoComposition::Topology(structure),
                    ..
                } = &*modification
                {
                    data.push([modification.to_string(), structure.to_string()])
                } else if let SimpleModificationInner::Gno {
                    composition: GnoComposition::Composition(composition),
                    ..
                } = &*modification
                {
                    data.push([
                        modification.to_string(),
                        composition
                            .iter()
                            .map(|(sug, amount)| format!("{sug}{amount}"))
                            .join(""),
                    ])
                }
            }
            if data.len() > 1 {
                table(
                    &data,
                    true,
                    &[Styling::with_fg(Some(Color::Magenta)), Styling::none()],
                );
            } else {
                println!("{}", "No modifications found".red())
            }
        }
        modification => display_single_mod(modification, precision),
    }
}

fn display_single_mod(modification: &SimpleModificationInner, precision: Option<usize>) {
    println!(
        "Full mass: {} {} {} {}",
        display_mass(modification.formula().monoisotopic_mass(), true, precision),
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
        SimpleModificationInner::Database {
            specificities, id, ..
        } => {
            display_id(id);
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
        SimpleModificationInner::Linker {
            specificities,
            id,
            length,
            ..
        } => {
            display_id(id);
            if let Some(length) = length {
                println!("Length: {}", length);
            }
            println!("Placement rules: ");
            for specificity in specificities {
                match specificity {
                    LinkerSpecificity::Symmetric(locations, stubs, diagnostic) => {
                        print!("  Locations: ");
                        display_placement_rules(locations);
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
        SimpleModificationInner::Gno {
            composition,
            id,
            structure_score,
            subsumption_level,
            motif,
            taxonomy,
            glycomeatlas,
        } => {
            display_id(id);
            if let Some(score) = structure_score {
                println!("Structure score: {}", score.to_string().blue());
            }
            println!("Subsumption: {}", subsumption_level.to_string().green());
            println!(
                "Motif: {}",
                motif
                    .iter()
                    .map(|(name, id)| format!("{name}:{id}"))
                    .join(", ")
            );
            println!(
                "Taxonomy: {}",
                taxonomy
                    .iter()
                    .map(|(name, id)| format!("{name}:{id}"))
                    .join(", ")
            );
            println!(
                "Glycomeatlas: {}",
                glycomeatlas
                    .iter()
                    .map(|(species, places)| format!(
                        "{species}:{}",
                        places
                            .iter()
                            .map(|(place, id)| format!("{place}({id})"))
                            .join(", ")
                    ))
                    .join("\n")
            );
            match composition {
                GnoComposition::Weight(mass) => {
                    println!(
                        "Average weight: {}",
                        display_mass(mass.into_inner(), true, precision)
                    )
                }
                GnoComposition::Composition(composition) => {
                    println!(
                        "Composition: {}",
                        composition
                            .iter()
                            .map(|(sug, amount)| format!("{}{amount}", sug.to_string().green()))
                            .join("")
                    )
                }
                GnoComposition::Topology(structure) => {
                    println!("Structure: {}", structure.to_string().green())
                }
            }
        }
        _ => (),
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
        "Ontology: {}, name: {}{}",
        id.ontology.to_string().purple(),
        id.name.green(),
        id.id.map_or(String::new(), |id| format!(
            ", index: {}",
            id.to_string().blue()
        ))
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
    let scoring = AlignScoring::<'static> {
        matrix: rustyms::align::matrix::BLOSUM90,
        ..Default::default()
    };
    let alignment = rustyms::align::align::<1, UnAmbiguous, UnAmbiguous>(
        allele.sequence,
        allele.sequence,
        scoring,
        rustyms::align::AlignType::GLOBAL,
    );
    if args.display_fasta {
        println!(
            ">{} {} {}",
            allele.name().purple(),
            allele.species.scientific_name(),
            allele.species.common_name().purple(),
        );
    } else {
        println!(
            "{} {} {}",
            allele.species.scientific_name().to_string().purple(),
            allele.species.common_name(),
            format!("{} / {}", allele.name(), allele.fancy_name()).purple(),
        );
    }
    show_annotated_mass_alignment(
        &alignment,
        Some(&allele),
        true,
        args.display_fasta,
        ("", ""),
        args,
    );
}

fn align<'a, A: AtMax<SimpleLinear>, B: AtMax<SimpleLinear>>(
    seq_a: &'a Peptidoform<A>,
    seq_b: &'a Peptidoform<B>,
    scoring: AlignScoring<'a>,
    ty: AlignType,
    kind: AlignmentKind,
) -> Alignment<'a, A, B> {
    if kind.normal {
        rustyms::align::align::<1, A, B>(seq_a, seq_b, scoring, ty)
    } else if kind.mass_based_huge {
        rustyms::align::align::<{ u16::MAX }, A, B>(seq_a, seq_b, scoring, ty)
    } else if kind.mass_based_long {
        rustyms::align::align::<8, A, B>(seq_a, seq_b, scoring, ty)
    } else {
        rustyms::align::align::<4, A, B>(seq_a, seq_b, scoring, ty)
    }
}

fn consecutive_align(
    seq: &Peptidoform<SimpleLinear>,
    species: Option<HashSet<imgt::Species>>,
    chains: Option<HashSet<imgt::ChainType>>,
    allele: imgt::AlleleSelection,
    scoring: AlignScoring,
    return_number: usize,
    kind: AlignmentKind,
) -> ConsecutiveAlignment<'static, SimpleLinear> {
    if kind.normal {
        par_consecutive_align::<1, SimpleLinear>(
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
            allele,
            scoring,
            return_number,
        )
    } else if kind.mass_based_huge {
        par_consecutive_align::<{ u16::MAX }, SimpleLinear>(
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
            allele,
            scoring,
            return_number,
        )
    } else if kind.mass_based_long {
        par_consecutive_align::<8, SimpleLinear>(
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
            allele,
            scoring,
            return_number,
        )
    } else {
        par_consecutive_align::<4, SimpleLinear>(
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
            allele,
            scoring,
            return_number,
        )
    }
}
