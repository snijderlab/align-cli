use clap::{Args, Parser};
use imgt_germlines::{AlleleSelection, Gene, GeneType, ChainType, Selection, Species};
use rustyms::{
    modification::ReturnModification,
    placement_rule::*,
    AminoAcid,  ComplexPeptide, LinearPeptide, MassTolerance, Modification,
};
use std::{collections::HashSet, fmt::Display};

#[derive(Parser, Debug)]
#[command(author, version, about)]
#[command(long_about = "It supports multiple use cases:

1) Align two sequences `align <X> <Y>`, this shows the best alignment for these two sequences.

2) Align a single peptide to the IMGT database, `align <X> --imgt <SELECTION>` in this way you can select either a specific germline or a search across the whole database.

3) Align a single peptide to a database `align <X> --file <FILE.fasta>`, this shows the scores for the best matches for this peptide alongside the alignment for the best match.

4) Get information about a single sequence `align <sequence>`, this shows many basic properties (like mass) and generates isobaric sequences to this sequence.

5) Get information about a single modification, `align --modification <MODIFICATION>`, this shows basic properties, and if it is a mass shift, eg `+58.01`, it shows all modifications that have the same mass within the tolerance.

6) Get the sequence of one or more germlines `align --imgt <SELECTION>`.")]
pub struct Cli {
    /// First sequence
    #[arg()]
    pub x: Option<String>,

    /// The selection of second sequence, can only be one of these
    #[command(flatten)]
    pub second: SecondSelection,

    /// The type of alignment, can only be one of these
    #[command(flatten)]
    pub alignment_type: AlignmentType,

    /// Use normal alignment (instead of the default of Mass alignment) this uses Smith Waterman or Needleman Wunsch algorithms (based on the alignment mode)
    /// using the BLOSUM62 scoring table.
    #[arg(long)]
    pub normal: bool,

    /// The number of characters to show on a single line in the alignment
    #[arg(short = 'n', long, default_value_t = 50)]
    pub line_width: usize,

    /// Toggles the showing of additional context for the match (only useful on non global alignments)
    #[arg(short = 'c', long)]
    pub context: bool,

    /// The number of hits to show in the tables for file and IMGT alignment
    #[arg(short = 'N', long, default_value_t = 10)]
    pub number_of_hits: usize,

    /// The maximal number of isobaric sets the generate, use `all` to generate all options
    #[arg(short, long, default_value_t = IsobaricNumber::Limited(25), value_parser=options_parse)]
    pub isobaric: IsobaricNumber,

    /// All possible fixed modifications that will be used in the isobaric sets generation, separated by commas `,`, commas can be
    /// escaped by wrapping the entire modification in square brackets `[..]`.
    /// You can overwrite the default placement rules in the same was as for variable modifications.
    #[arg(short = 'F', long, default_value_t = Modifications::None, value_parser=modifications_parse, allow_hyphen_values=true)]
    pub fixed: Modifications,

    /// All possible variable modifications that will be used in the isobaric sets generation, separated by commas `,`, commas can be
    /// escaped by wrapping the entire modification in square brackets `[..]`.
    /// You can overwrite the default placement rules in the following way: `@AA-pos` where `AA` is the list of all amino acids it can be
    /// placed or a star to indicate it can be placed on all locations, and `pos` is the position: * -> Anywhere,
    /// N/n -> N terminal (protein/peptide), C/c -> C terminal (protein/peptide). The position can be left out which defaults to Anywhere.
    /// Examples for the rules: `Carboxymethyl@C`, `Oxidation@WFH`, `Amidated@*-C`.
    #[arg(short, long, default_value_t = Modifications::None, value_parser=modifications_parse, allow_hyphen_values=true)]
    pub variable: Modifications,

    /// The base to always include in generating isobaric sets. This is assumed to be a simple sequence (for details see rustyms::LinearPeptide::assume_simple).
    #[arg(long, value_parser=peptide_parser)]
    pub include: Option<LinearPeptide>,

    /// Overrule the default set of amino acids used in the isobaric sequences generation. The default set has all amino acids with a defined mass (no I/L in favour of J, no B/Z/X, but with U/O included).
    #[arg(long, value_parser=amino_acids_parser)]
    pub amino_acids: Option<AminoAcids>,

    /// The tolerance for the isobaric set search and the definition for isobaric sets in the alignment, use `<x>ppm` or `<x>da` to control the unit, eg `10.0ppm` or `2.3da`
    #[arg(short, long, default_value_t = MassTolerance::Ppm(10.0), value_parser=mass_tolerance_parse)]
    pub tolerance: MassTolerance,

    /// A modification you want details on, if it is a mass shift modification eg `+58.01` it will show all predefined modifications that are within the tolerance of this mass
    #[arg(short, long, value_parser=modification_parse, allow_hyphen_values=true)]
    pub modification: Option<Modification>,

    /// To turn on fancy display of germlines and compositions
    #[arg(long)]
    pub fancy: bool,
}

#[test]
fn verify_cli() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}

#[derive(Args, Debug)]
#[group(multiple = false)]
pub struct AlignmentType {
    /// Use global alignment [default]
    #[arg(short, long)]
    pub global: bool,

    /// Use semi-global alignment, meaning that the second sequence has to match fully, while the first sequence can be longer then the alignment.
    /// When the `--file` or `--imgt` mode is used this flag indicates that the given sequence can align semi globally to the provided database sequences.
    #[arg(short, long)]
    pub semi_global: bool,

    /// Use local alignment
    #[arg(short, long)]
    pub local: bool,
}

impl AlignmentType {
    pub fn ty(&self) -> rustyms::align::Type {
        if self.local {
            rustyms::align::Type::Local
        } else if self.semi_global {
            rustyms::align::Type::GlobalForB
        } else {
            rustyms::align::Type::Global
        }
    }
}

#[derive(Args, Debug)]
#[group(multiple = false)]
pub struct SecondSelection {
    /// Second sequence
    #[arg()]
    pub y: Option<String>,

    /// A fasta database file to open to align the sequence to, only provide a single sequence for this mode
    #[arg(short, long)]
    pub file: Option<String>,

    /// Align against IMGT germline sequences. 
    /// 
    /// You can select either a specific germline using `species:<SPECIES>&<NAME>` with if needed the allele specified using `allele:<NUMBER>`, otherwise it defaults to the first allele. 
    /// 
    /// Or you can give a selection using `species:<SPECIES>`, `kind:<KIND>`, `segment:<SEGMENT>`, `allele:<all|first>`. 
    /// Any of these criteria can be left out to select all of them, except for the allele, leaving that out will select only the first. Combine criteria using `&`.
    #[arg(long, value_parser=imgt_selection_parse)]
    pub imgt: Option<IMGTSelection>,
}

#[derive(Debug, Clone)]
pub enum IMGTSelection {
    Gene(Species, Gene, Option<usize>),
    Search(Selection),
}

fn imgt_selection_parse(s: &str) -> Result<IMGTSelection, String> {
    let species = s.split('&')
        .find_map(|p| p.strip_prefix("species:"))
        .map(|s| {
            s.split(',')
                .map(|s| {
                    s.parse()
                        .map_err(|()| format!("Not a recognised species: {s}"))
                })
                .collect::<Result<HashSet<Species>, String>>()
        })
        .transpose()?;
    let chains = s.split('&')
        .find_map(|p| p.strip_prefix("chain:"))
        .map(|s| {
            s.split(',')
                .map(|s| {
                    s.parse()
                        .map_err(|()| format!("Not a recognised chain: {s}"))
                })
                .collect::<Result<HashSet<ChainType>, String>>()
        })
        .transpose()?;
    let genes = s.split('&')
        .find_map(|p| p.strip_prefix("gene:"))
        .map(|s| {
            s.split(',')
                .map(|s| {
                    s.parse()
                        .map_err(|()| format!("Not a recognised gene: {s}"))
                })
                .collect::<Result<HashSet<GeneType>, String>>()
        })
        .transpose()?;
    let name = s.split('&')
        .find_map(|p| p.strip_prefix("name:"))
        .or(s.split('&').find(|p| p.starts_with("IG")));

    if let Some(name) = name {
        let species = species
            .ok_or("No species specified")?
            .into_iter()
            .collect::<Vec<_>>();
        let allele = s.split('&')
            .find_map(|p| p.strip_prefix("allele:"))
            .and_then(|s| match s {
                "first" => None,
                other => Some(
                    other
                        .parse::<usize>()
                        .map_err(|_| format!("Not a valid number for allele selection: {other}")),
                ),
            }).transpose()?;
        if species.len() != 1 {
            Err("You have to specify a single species for IMGT gene selection")?
        } else {
            Ok(IMGTSelection::Gene(
                species[0],
                Gene::from_imgt_name(name)?,
                allele,
            ))
        }
    } else {
        let allele = s.split('&')
            .find_map(|p| p.strip_prefix("allele:"))
            .map(|s| match s {
                "all" => Ok(AlleleSelection::All),
                "first" => Ok(AlleleSelection::First),
                err => Err(format!("Not a valid allele specification: {err}")),
            })
            .transpose()?
            .unwrap_or(AlleleSelection::First);
        Ok(IMGTSelection::Search(Selection {
            species,
            chains,
            genes,
            allele,
        }))
    }
}

#[derive(Debug, Clone)]
pub enum IsobaricNumber {
    All,
    Limited(usize),
}
impl Display for IsobaricNumber {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::All => write!(f, "all"),
            Self::Limited(limit) => write!(f, "{limit}"),
        }
    }
}
fn mass_tolerance_parse(input: &str) -> Result<MassTolerance, &'static str> {
    input.parse().map_err(|()| "Invalid tolerance parameter")
}
fn options_parse(input: &str) -> Result<IsobaricNumber, &'static str> {
    if input.to_lowercase() == "all" {
        Ok(IsobaricNumber::All)
    } else {
        input
            .parse::<usize>()
            .map(IsobaricNumber::Limited)
            .map_err(|_| "Invalid options parameter")
    }
}
fn peptide_parser(input: &str) -> Result<LinearPeptide, String> {
    Ok(ComplexPeptide::pro_forma(input)
        .map_err(|e| e.to_string())?
        .assume_linear()
        .assume_simple())
}
fn amino_acids_parser(input: &str) -> Result<AminoAcids, String> {
    input
        .chars()
        .map(|c| AminoAcid::try_from(c).map_err(|()| format!("`{c}` is not a valid amino acid")))
        .collect()
}
type AminoAcids = Vec<AminoAcid>;
#[derive(Debug, Clone)]
pub enum Modifications {
    None,
    Some(Vec<(Modification, Option<PlacementRule>)>),
}
impl Modifications {
    pub fn mods(&self) -> &[(Modification, Option<PlacementRule>)] {
        match self {
            Self::None => &[],
            Self::Some(m) => m,
        }
    }
}
impl Display for Modifications {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::None => write!(f, ""),
            Self::Some(mods) => {
                let mut start = true;
                for m in mods {
                    write!(f, "{}{}", if start { "" } else { "," }, m.0).unwrap();
                    start = false;
                }
                Ok(())
            }
        }
    }
}
fn modifications_parse(input: &str) -> Result<Modifications, String> {
    fn parse_position(pos: &str) -> Result<Position, String> {
        match pos {
            "*" => Ok(Position::Anywhere),
            "C" => Ok(Position::ProteinCTerm),
            "c" => Ok(Position::AnyCTerm),
            "N" => Ok(Position::ProteinNTerm),
            "n" => Ok(Position::AnyNTerm),
            _ => Err(format!(
                "'{pos}' is not a valid modification placement position use any of: */N/n/C/c"
            )),
        }
    }
    fn parse_aa(aa: &str) -> Result<Option<Vec<AminoAcid>>, String> {
        if aa == "*" {
            Ok(None)
        } else {
            Ok(Some(
                aa.chars()
                    .map(|c| {
                        AminoAcid::try_from(c)
                            .map_err(|_| format!("'{c}' is not a valid amino acid"))
                    })
                    .collect::<Result<Vec<_>, _>>()?,
            ))
        }
    }
    fn split(input: &str) -> Vec<&str> {
        let input = input.trim_end_matches(',');
        let mut index = None;
        let mut depth = 0;
        let mut res = Vec::new();
        for (i, c) in input.as_bytes().iter().enumerate() {
            match c {
                b'[' => {
                    if index.is_none() || depth == 0 {
                        index = Some(i + 1);
                        depth = 1;
                    } else {
                        depth += 1;
                    }
                }
                b']' if index.is_some() && depth > 0 => {
                    if depth == 1 {
                        res.push(&input[index.unwrap()..i]);
                        index = None;
                        depth = 0;
                    } else {
                        depth -= 1;
                    }
                }
                b',' if depth == 0 => {
                    if let Some(ind) = index {
                        res.push(&input[ind..i]);
                    }
                    index = Some(i + 1);
                }
                _ if index.is_none() && !c.is_ascii_whitespace() => index = Some(i),
                _ => (),
            }
        }
        if let Some(ind) = index {
            if ind != input.len() {
                res.push(&input[ind..]);
            }
        }
        res
    }

    if input.is_empty() {
        Ok(Modifications::None)
    } else {
        split(input).into_iter()
            .map(|m| {
                if let Some((head, tail)) = m.split_once('@') {
                    let modification = 
                    Modification::try_from(head, 0..head.len(), &mut Vec::new()).map_err(|e| e.to_string()).and_then(|m| if let Some(d) = m.defined() {
                        Ok(d) } else {
                            Err("Can not define ambiguous modifications for the modifications parameter".to_string())
                        }
                    )?;
                    let rule = if let Some((aa, position)) = tail.split_once('-') {
                        if let Some(aa) = parse_aa(aa)? {
                            PlacementRule::AminoAcid(aa, parse_position(position)?)
                        } else {
                            PlacementRule::Terminal(parse_position(position)?)
                        }
                    } else if let Some(aa) = parse_aa(tail)? {
                            PlacementRule::AminoAcid(aa, Position::Anywhere)
                        } else {
                            return Err("Cannot have a modification rule that allows a modification on all position on all amino acids".to_string())
                        };
                    Ok((modification, Some(rule)))
                } else {
                    Modification::try_from(m, 0..m.len(), &mut Vec::new()).map_err(|e| e.to_string()).and_then(|m| if let Some(d) = m.defined() {
                        Ok((d, None)) } else {
                            Err("Can not define ambiguous modifications for the modifications parameter".to_string())
                        }
                    )
                }
                })
            .collect::<Result<Vec<(Modification, Option<PlacementRule>)>, String>>()
            .map(Modifications::Some)
    }
}

fn modification_parse(input: &str) -> Result<Modification, String> {
    if input.is_empty() {
        Err("Empty".to_string())
    } else {
        Modification::try_from(input, 0..input.len(), &mut Vec::new())
            .map(|m| match m {
                ReturnModification::Defined(d) => d,
                _ => {
                    panic!("Can not define ambiguous modifications for the modifications parameter")
                }
            })
            .map_err(|err| err.to_string())
    }
}