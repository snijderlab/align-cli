use colored::Color;

pub trait Legend {
    fn fg_color(&self) -> Option<Color>;
    fn bg_color(&self) -> Option<Color>;
}

impl Legend for mzcore::sequence::Annotation {
    fn fg_color(&self) -> Option<Color> {
        match self {
            Self::Conserved => Some(Color::Blue),
            Self::NGlycan => Some(Color::Green), // TODO: If on CDR2 not visible
            Self::Other(_) => None,
        }
    }
    fn bg_color(&self) -> Option<Color> {
        None
    }
}

impl Legend for mzcore::sequence::Region {
    fn fg_color(&self) -> Option<Color> {
        match self {
            Self::ComplementarityDetermining(_) => Some(Color::Black),
            _ => None,
        }
    }
    fn bg_color(&self) -> Option<Color> {
        match self {
            Self::ComplementarityDetermining(1) => Some(Color::Red),
            Self::ComplementarityDetermining(2) => Some(Color::Green),
            Self::ComplementarityDetermining(3) => Some(Color::Blue),
            _ => None,
        }
    }
}
