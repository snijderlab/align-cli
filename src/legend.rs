use colored::Color;
use rustyms::identification;

pub trait Legend {
    fn fg_color(&self) -> Option<Color>;
    fn bg_color(&self) -> Option<Color>;
}

impl Legend for identification::Annotation {
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

impl Legend for identification::Region {
    fn fg_color(&self) -> Option<Color> {
        match self {
            Self::ComplementarityDeterminingRegion(_) => Some(Color::Black),
            _ => None,
        }
    }
    fn bg_color(&self) -> Option<Color> {
        match self {
            Self::ComplementarityDeterminingRegion(1) => Some(Color::Red),
            Self::ComplementarityDeterminingRegion(2) => Some(Color::Green),
            Self::ComplementarityDeterminingRegion(3) => Some(Color::Blue),
            _ => None,
        }
    }
}
