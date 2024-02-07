use colored::Color;

pub trait Legend {
    fn fg_color(&self) -> Option<Color>;
    fn bg_color(&self) -> Option<Color>;
}

impl Legend for imgt_germlines::Annotation {
    fn fg_color(&self) -> Option<Color> {
        Some(match self {
            Self::Cysteine1 => Color::Blue,
            Self::Cysteine2 => Color::Blue,
            Self::Tryptophan => Color::Blue,
            Self::Phenylalanine => Color::Blue,
            Self::Glycine => Color::Blue,
            Self::NGlycan => Color::Green, // TODO: If on CDR2 not visible
        })
    }
    fn bg_color(&self) -> Option<Color> {
        None
    }
}

impl Legend for imgt_germlines::Region {
    fn fg_color(&self) -> Option<Color> {
        match self {
            Self::CDR1 | Self::CDR2 | Self::CDR3 => Some(Color::Black),
            _ => None,
        }
    }
    fn bg_color(&self) -> Option<Color> {
        match self {
            Self::CDR1 => Some(Color::Red),
            Self::CDR2 => Some(Color::Green),
            Self::CDR3 => Some(Color::Blue),
            _ => None,
        }
    }
}
