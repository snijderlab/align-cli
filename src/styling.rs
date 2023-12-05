use colored::{Color, ColoredString, Colorize, Styles};

#[derive(Default, Clone)]
pub struct Styling {
    fg: Option<Color>,
    bg: Option<Color>,
    styles: Vec<Styles>,
}

#[allow(dead_code)]
impl Styling {
    pub const fn none() -> Self {
        Self {
            fg: None,
            bg: None,
            styles: Vec::new(),
        }
    }
    pub fn with_fg(color: Option<Color>) -> Self {
        Self {
            fg: color,
            ..Self::none()
        }
    }
    pub fn with_bg(color: Option<Color>) -> Self {
        Self {
            bg: color,
            ..Self::none()
        }
    }
    pub fn with_style(style: Styles) -> Self {
        Self {
            styles: vec![style],
            ..Self::none()
        }
    }
    pub fn fg(self, color: Option<Color>) -> Self {
        Self { fg: color, ..self }
    }
    pub fn bg(self, color: Option<Color>) -> Self {
        Self { bg: color, ..self }
    }
    /// If there is no color set already, set the given color
    pub fn or_fg(self, color: Option<Color>) -> Self {
        Self {
            fg: self.fg.or(color),
            ..self
        }
    }
    /// If there is no color set already, set the given color
    pub fn or_bg(self, color: Option<Color>) -> Self {
        Self {
            bg: self.bg.or(color),
            ..self
        }
    }
    pub fn style(mut self, style: Styles) -> Self {
        self.styles.push(style);
        self
    }
    pub fn maybe_style(self, style: Option<Styles>) -> Self {
        if let Some(style) = style {
            self.style(style)
        } else {
            self
        }
    }
}

pub trait ExtendedColorize {
    type Output: ExtendedColorize<Output = Self::Output> + Colorize;
    fn apply_style(&self, style: Option<Styles>) -> Self::Output;
    fn on_color_e(&self, color: Option<Color>) -> Self::Output;
    fn color_e(&self, color: Option<Color>) -> Self::Output;
    fn apply(&self, styles: &Styling) -> Self::Output {
        styles.styles.iter().fold(
            self.color_e(styles.fg).on_color_e(styles.bg),
            |output, style| output.apply_style(Some(*style)),
        )
    }
}

impl ExtendedColorize for ColoredString {
    type Output = Self;
    fn apply_style(&self, style: Option<Styles>) -> Self {
        if let Some(style) = style {
            match style {
                Styles::Clear => self.clone().clear(),
                Styles::Bold => self.clone().bold(),
                Styles::Dimmed => self.clone().dimmed(),
                Styles::Underline => self.clone().underline(),
                Styles::Reversed => self.clone().reversed(),
                Styles::Italic => self.clone().italic(),
                Styles::Blink => self.clone().blink(),
                Styles::Hidden => self.clone().hidden(),
                Styles::Strikethrough => self.clone().strikethrough(),
            }
        } else {
            self.clone()
        }
    }
    fn on_color_e(&self, color: Option<Color>) -> Self {
        match color {
            Some(clr) => self.clone().on_color(clr),
            None => self.clone(),
        }
    }
    fn color_e(&self, color: Option<Color>) -> Self {
        match color {
            Some(clr) => self.clone().color(clr),
            None => self.clone(),
        }
    }
}

impl ExtendedColorize for &str {
    type Output = ColoredString;
    fn apply_style(&self, style: Option<Styles>) -> Self::Output {
        if let Some(style) = style {
            match style {
                Styles::Clear => self.clear(),
                Styles::Bold => self.bold(),
                Styles::Dimmed => self.dimmed(),
                Styles::Underline => self.underline(),
                Styles::Reversed => self.reversed(),
                Styles::Italic => self.italic(),
                Styles::Blink => self.blink(),
                Styles::Hidden => self.hidden(),
                Styles::Strikethrough => self.strikethrough(),
            }
        } else {
            self.normal()
        }
    }
    fn on_color_e(&self, color: Option<Color>) -> Self::Output {
        match color {
            Some(clr) => self.on_color(clr),
            None => self.normal(),
        }
    }
    fn color_e(&self, color: Option<Color>) -> Self::Output {
        match color {
            Some(clr) => self.color(clr),
            None => self.normal(),
        }
    }
}

impl ExtendedColorize for String {
    type Output = ColoredString;
    fn apply_style(&self, style: Option<Styles>) -> Self::Output {
        self.as_str().apply_style(style)
    }
    fn on_color_e(&self, color: Option<Color>) -> Self::Output {
        self.as_str().on_color_e(color)
    }
    fn color_e(&self, color: Option<Color>) -> Self::Output {
        self.as_str().color_e(color)
    }
}

impl ExtendedColorize for char {
    type Output = ColoredString;
    fn apply_style(&self, style: Option<Styles>) -> Self::Output {
        self.to_string().apply_style(style)
    }
    fn on_color_e(&self, color: Option<Color>) -> Self::Output {
        self.to_string().on_color_e(color)
    }
    fn color_e(&self, color: Option<Color>) -> Self::Output {
        self.to_string().color_e(color)
    }
}
