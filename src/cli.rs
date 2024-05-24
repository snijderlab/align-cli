use clap::{Args, Parser};
use rustyms::imgt::{AlleleSelection, Gene, GeneType, ChainType, Species};
use rustyms::modification::SimpleModification;
use rustyms::system::Mass;
use rustyms::{ReturnModification, Simple};
use rustyms::{
    placement_rule::*,
    AminoAcid,  LinearPeptide, Tolerance, align::{AlignType, self},
};
use std::str::FromStr;
use std::{collections::HashSet, fmt::Display};

#[derive(Parser, Debug)]
#[command(author, version, about)]
#[command(long_about = "It supports multiple use cases:

1. Pairwise alignment
   - Align two sequences `align <A> <B>`, this shows the best alignment for these two sequences.
   - Align a single peptide to a database `align <A> --file <FILE.fasta>`.
   - Align a single peptide to the IMGT database `align <A> --imgt`.
   - Align a single peptide to the V-J-C domains in the IMGT database `align <A> --domain`.
   - Align a single peptide to a specific gene in IMGT database `align <A> --specific-gene <GENE>`.

2. Get information about a single sequence `align <sequence>`, this shows many basic properties (like mass) and generates isobaric sequences to this sequence.
   - Use `--fixed <MODIFICATIONS>` and `--variable <MODIFICATIONS>` to fine tune the generated isobaric sequences.

3. Get information about a single modification `align --modification <MODIFICATION>`.
   - Use a full name to list its properties eg `--modification Oxidation`
   - Use a formula to find all modifications with that formula eg `--modification Formula:O`
   - Use a mass to find all modifications with that mass eg `--modification +15.995`

4. List IMGT genes `align --imgt` or `align --specific-gene <GENE>`.")]
pub struct Cli {
    /// First sequence
    #[arg()]
    pub a: Option<String>,

    /// The selection of second sequence, can only be one of these
    #[command(flatten)]
    pub second: SecondSelection,

    /// The kind of alignment (normal/mass based/etc), can only be one of these
    #[command(flatten)]
    pub alignment_kind: AlignmentKind,

    /// The type of alignment, can only be one of these
    #[command(flatten)]
    pub alignment_type: AlignmentType,

    /// The scoring matrix used, forced to be only one of the possible values
    #[command(flatten)]
    pub scoring_matrix: ScoringMatrix,

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
    /// You can overwrite the default placement rules in the same way as for variable modifications.
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
    pub include: Option<LinearPeptide<Simple>>,

    /// Overrule the default set of amino acids used in the isobaric sequences generation. The default set has all amino acids with a defined mass (no I/L in favour of J, no B/Z/X, but with U/O included).
    #[arg(long, value_parser=amino_acids_parser)]
    pub amino_acids: Option<AminoAcids>,

    /// The tolerance for the isobaric set search and the definition for isobaric sets in the alignment, use `<x>ppm` or `<x>da` to control the unit, eg `10.0ppm` or `2.3da`
    #[arg(short, long, default_value_t = Tolerance::new_ppm(10.0.into()), value_parser=mass_tolerance_parse)]
    pub tolerance: Tolerance<Mass>,

    /// A modification you want details on, if it is a mass shift modification eg `+58.01` it will show all predefined modifications that are within the tolerance of this mass
    #[arg(short, long, value_parser=modification_parse, allow_hyphen_values=true)]
    pub modification: Option<SimpleModification>,

    /// The species selected for any IMGT based alignments, you can use either the common or scientific name for the species.
    #[arg(long)]
    pub species: Option<Species>,

    /// The chains selected for any IMGT based alignments, you can use any number of H, K, L, and I.
    #[arg(long, value_parser=chains_parser)]
    pub chains: Option<HashSet<ChainType>>,

    /// The genes selected for any IMGT based alignments, you can use any number of V, J, C, A, D, E, G, M, O, and T.
    #[arg(long, value_parser=genes_parser)]
    pub genes: Option<HashSet<GeneType>>,

    /// The genes selected for any IMGT based alignments, you can use either 'all' or 'first'.
    #[arg(long, value_parser=allele_parser, default_value = "first")]
    pub allele: AlleleSelection,

    /// Show full mass precision according to floating point math instead of the normal capped number of digits
    #[arg(long)]
    pub full_number: bool,
}

fn chains_parser(value: &str) -> Result<HashSet<ChainType>, String> {
    let mut set = HashSet::new();
    for c in value.chars() {
        set.insert(ChainType::from_str(c.to_string().as_str()).map_err(|()| format!("Not a valid chain type: {c}"))?);
    }
    Ok(set)
}

fn genes_parser(value: &str) -> Result<HashSet<GeneType>, String> {
    let mut set = HashSet::new();
    for c in value.chars() {
        set.insert(GeneType::from_str(c.to_string().as_str()).map_err(|()| format!("Not a valid gene type: {c}"))?);
    }
    Ok(set)
}

fn allele_parser(value: &str) -> Result<AlleleSelection, String> {
    match value.trim().to_lowercase().as_str() {
        "all" => Ok(AlleleSelection::All),
        "first" => Ok(AlleleSelection::First),
        _ => Err(format!("Not a valid allele selection: {value}, use 'all' or 'first'."))
    }   
}

#[test]
fn verify_cli() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}

#[derive(Args, Debug, Clone, Copy)]
#[group(multiple = false)]
pub struct AlignmentKind {
    /// Use normal alignment (instead of the default of Mass alignment) this uses Smith Waterman or Needleman Wunsch algorithms (based on the alignment mode)
    /// using the same modified BLOSUM62 scoring table as used in mass based alignment. Note: this is the same mass based alignment algorithm but set to a
    /// maximal length of isobaric sets of 1, meaning it will still handle modifications and show I/L as isobaric.
    #[arg(long)]
    pub normal: bool,

    /// Do mass based alignment but allow for a maximal isobaric set length of 8 instead of the default 4.
    #[arg(long)]
    pub mass_based_long: bool,

    /// Do mass based alignment but allow for an unbounded maximal isobaric set instead of the default 4.
    #[arg(long)]
    pub mass_based_huge: bool,
}

#[derive(Args, Debug)]
#[group(multiple = false)]
pub struct ScoringMatrix {
    /// BLOSUM45 matrix
    #[arg(long)]
    pub blosum45: bool,
    /// BLOSUM50 matrix
    #[arg(long)]
    pub blosum50: bool,
    /// BLOSUM62 matrix [default]
    #[arg(long)]
    pub blosum62: bool,
    /// BLOSUM80 matrix
    #[arg(long)]
    pub blosum80: bool,
    /// BLOSUM90 matrix
    #[arg(long)]
    pub blosum90: bool,
    /// Identity matrix
    #[arg(long)]
    pub identity: bool,
    /// PAM30 matrix
    #[arg(long)]
    pub pam30: bool,
    /// PAM70 matrix
    #[arg(long)]
    pub pam70: bool,
    /// PAM250 matrix
    #[arg(long)]
    pub pam250: bool,
}

impl ScoringMatrix {
    pub fn matrix(&self) -> &[[i8; AminoAcid::TOTAL_NUMBER]; AminoAcid::TOTAL_NUMBER] {
        if self.blosum45 {
            align::matrix::BLOSUM45
        } else if self.blosum50 {
            align::matrix::BLOSUM50
        } else if self.blosum62 {
            align::matrix::BLOSUM62
        } else if self.blosum80 {
            align::matrix::BLOSUM80
        } else if self.blosum90 {
            align::matrix::BLOSUM90
        } else if self.identity {
            align::matrix::IDENTITY
        } else if self.pam30 {
            align::matrix::PAM30
        } else if self.pam70 {
            align::matrix::PAM70
        } else if self.pam250 {
            align::matrix::PAM250
        } else {
            align::matrix::BLOSUM62
        }
    }
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

    /// Use semi-global alignment, meaning that the first sequence has to match fully, while the second sequence can be longer then the alignment.
    /// When the `--file` or `--imgt` mode is used this flag indicates that the database sequences can align semi globally to the provided sequence.
    #[arg(short='S', long)]
    pub semi_global_a: bool,

    /// Use local alignment
    #[arg(short, long)]
    pub local: bool,

    /// Specify the type fully. Specify each position as local `0` or global `1` in the following order: left A, left B, right A, right B. 
    /// An either global can be specified by putting a hyphen on the left or right side, `-xx` is either global left, `xx-` is either global right, `--` is either global.
    /// For example `1001` means global on left A and right B which will make peptide A extend peptide B.
    #[arg(long, value_parser=type_parser, allow_hyphen_values=true)]
    pub r#type: Option<AlignType>,
}

impl AlignmentType {
    pub fn ty(&self) -> rustyms::align::AlignType {
        if let Some(ty) = self.r#type {
            ty
        } else if self.local {
            rustyms::align::AlignType::LOCAL
        } else if self.semi_global {
            rustyms::align::AlignType::GLOBAL_B
        } else if self.semi_global_a {
            rustyms::align::AlignType::GLOBAL_A
        } else {
            rustyms::align::AlignType::GLOBAL
        }
    }
}

#[derive(Args, Debug)]
#[group(multiple = false)]
pub struct SecondSelection {
    /// Second sequence
    #[arg()]
    pub b: Option<String>,

    /// A fasta database file to open to align the sequence to, only provide a single sequence for this mode
    #[arg(short, long)]
    pub file: Option<String>,

    /// A csv file of pairs of sequences to score it returns a csv file with statistics added as last columns. 
    /// The requirement is that the pair columns have to be called "a" and "b".
    #[arg(long)]
    pub csv: Option<String>,

    /// Align against IMGT germline sequences. Use species/chains/genes/allele to further specify the IMGT selection.
    #[arg(long)]
    pub imgt: bool,

    /// Align against one specific IMGT gene, using species is required if this is used.
    #[arg(long, value_parser=parse_specific_gene)]
    pub specific_gene: Option<(Gene, Option<usize>)>,

    /// Do a consecutive alignment against V-J-D (in that order) of the IMGT database. Use species/chains/genes/allele to further specify the IMGT selection.
    #[arg(long)]
    pub domain: bool,
}

fn parse_specific_gene(value: &str) -> Result<(Gene, Option<usize>), String> {
    Gene::from_imgt_name_with_allele(value)
        .map(|(g, a)| (g, Some(a)))
        .or_else(|_| Gene::from_imgt_name(value).map(|g| (g, None)))
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
fn mass_tolerance_parse(input: &str) -> Result<Tolerance<Mass>, &'static str> {
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
fn peptide_parser(input: &str) -> Result<LinearPeptide<Simple>, String> {
    LinearPeptide::pro_forma(input, None).map_err(|e| e.to_string())?.simple().ok_or("Not a simple peptide".to_string())
}
fn amino_acids_parser(input: &str) -> Result<AminoAcids, String> {
    input
        .chars()
        .map(|c| AminoAcid::try_from(c).map_err(|()| format!("`{c}` is not a valid amino acid")))
        .collect()
}
type AminoAcids = Vec<AminoAcid>;

fn type_parser(input: &str) -> Result<AlignType, String> {
    input.parse().map_err(|()| format!("Not a valid alignment type: '{input}'"))
}

#[derive(Debug, Clone)]
pub enum Modifications {
    None,
    Some(Vec<(SimpleModification, Option<PlacementRule>)>),
}
impl Modifications {
    pub fn mods(&self) -> &[(SimpleModification, Option<PlacementRule>)] {
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
                    SimpleModification::try_from(head, 0..head.len(), &mut Vec::new(), &mut Vec::new(), None).map_err(|e| e.to_string()).and_then(|m| if let Some(d) = m.defined() {
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
                    SimpleModification::try_from(m, 0..m.len(), &mut Vec::new(), &mut Vec::new(), None).map_err(|e| e.to_string()).and_then(|m| if let Some(d) = m.defined() {
                        Ok((d, None)) } else {
                            Err("Can not define ambiguous modifications for the modifications parameter".to_string())
                        }
                    )
                }
                })
            .collect::<Result<Vec<(SimpleModification, Option<PlacementRule>)>, String>>()
            .map(Modifications::Some)
    }
}

fn modification_parse(input: &str) -> Result<SimpleModification, String> {
    if input.is_empty() {
        Err("Empty".to_string())
    } else {
        SimpleModification::try_from(input, 0..input.len(), &mut Vec::new(), &mut Vec::new(), None)
            .map(|m| match m {
                ReturnModification::Defined(d) => d,
                _ => {
                    panic!("Can not define ambiguous modifications for the modifications parameter")
                }
            })
            .map_err(|err| err.to_string())
    }
}