//! Project files: an input OBJ, an output directory and the job settings, as
//! TOML. The GUI saves one with every export, and `mesher run project.toml`
//! replays it, so anything done in the GUI is reproducible from the CLI.
//!
//! ```toml
//! input = "model.obj"
//! output = "model_out"
//!
//! [settings]
//! parts = ["body", "lid"]
//! element_size = 0.03
//! ```
//!
//! Relative paths are relative to the project file's directory, and are
//! written that way (with `/`) when the path lies under it, so a project and
//! its model can move together.

use std::path::{Component, Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::protocol::Job;
use crate::report::JobConfig;

/// A saved job.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Project {
    /// The input OBJ file.
    pub input: PathBuf,
    /// Where the worker writes the report and meshes.
    pub output: PathBuf,
    #[serde(default)]
    pub settings: JobConfig,
}

/// Why a project file could not be used.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProjectError {
    /// The file could not be read or written.
    Io(String),
    /// The text is not a valid project.
    Invalid(String),
}

impl std::fmt::Display for ProjectError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ProjectError::Io(e) => write!(f, "{e}"),
            ProjectError::Invalid(e) => write!(f, "not a valid project: {e}"),
        }
    }
}

impl Project {
    /// Parse project text; relative paths are resolved against `base`.
    ///
    /// # Errors
    /// `Invalid` when the text is not a valid project.
    pub fn from_toml(text: &str, base: &Path) -> Result<Project, ProjectError> {
        let mut p: Project =
            toml::from_str(text).map_err(|e| ProjectError::Invalid(e.message().to_owned()))?;
        if p.input.as_os_str().is_empty() {
            return Err(ProjectError::Invalid("`input` is empty".to_owned()));
        }
        if p.output.as_os_str().is_empty() {
            return Err(ProjectError::Invalid("`output` is empty".to_owned()));
        }
        let s = &p.settings;
        if !(s.max_deviation > 0.0 && s.max_deviation.is_finite())
            || !(s.element_size > 0.0 && s.element_size.is_finite())
        {
            return Err(ProjectError::Invalid(
                "`max_deviation` and `element_size` must be positive".to_owned(),
            ));
        }
        p.input = base.join(&p.input);
        p.output = base.join(&p.output);
        Ok(p)
    }

    /// The project as TOML, with paths relative to `base` where they lie
    /// under it.
    #[must_use]
    pub fn to_toml(&self, base: &Path) -> String {
        let rel = Project {
            input: relative(&self.input, base),
            output: relative(&self.output, base),
            settings: self.settings.clone(),
        };
        toml::to_string(&rel).unwrap_or_default()
    }

    /// Read a project file.
    ///
    /// # Errors
    /// `Io` when the file cannot be read, `Invalid` when it is not a project.
    pub fn read(path: &Path) -> Result<Project, ProjectError> {
        let text = std::fs::read_to_string(path)
            .map_err(|e| ProjectError::Io(format!("cannot read {}: {e}", path.display())))?;
        Project::from_toml(&text, path.parent().unwrap_or(Path::new("")))
    }

    /// Write the project file.
    ///
    /// # Errors
    /// `Io` when the file cannot be written.
    pub fn write(&self, path: &Path) -> Result<(), ProjectError> {
        let text = self.to_toml(path.parent().unwrap_or(Path::new("")));
        std::fs::write(path, text)
            .map_err(|e| ProjectError::Io(format!("cannot write {}: {e}", path.display())))
    }

    /// The worker job this project describes.
    #[must_use]
    pub fn job(&self) -> Job {
        Job {
            input: self.input.display().to_string(),
            output_dir: self.output.display().to_string(),
            config: self.settings.clone(),
        }
    }
}

/// `path` relative to `base` when it lies under it (with `/` separators),
/// else unchanged.
fn relative(path: &Path, base: &Path) -> PathBuf {
    let clean = |p: &Path| -> PathBuf {
        p.components()
            .filter(|c| !matches!(c, Component::CurDir))
            .collect()
    };
    let (p, b) = (clean(path), clean(base));
    match p.strip_prefix(&b) {
        Ok(rest) if !rest.as_os_str().is_empty() => {
            let parts: Vec<String> = rest
                .components()
                .map(|c| c.as_os_str().to_string_lossy().into_owned())
                .collect();
            PathBuf::from(parts.join("/"))
        }
        _ => path.to_path_buf(),
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;

    #[test]
    fn round_trip_with_relative_paths() {
        let base = Path::new("projects/a");
        let p = Project {
            input: base.join("models/two.obj"),
            output: base.join("out"),
            settings: JobConfig {
                parts: vec!["lid".to_owned(), "body".to_owned()],
                strict: true,
                ..JobConfig::default()
            },
        };
        let text = p.to_toml(base);
        assert!(text.contains("input = \"models/two.obj\""), "{text}");
        assert!(text.contains("output = \"out\""), "{text}");
        assert_eq!(Project::from_toml(&text, base).unwrap(), p);
    }

    #[test]
    fn defaults_and_errors() {
        let p = Project::from_toml("input = \"m.obj\"\noutput = \"o\"\n", Path::new("")).unwrap();
        assert_eq!(p.settings, JobConfig::default());
        assert_eq!(p.job().input, "m.obj");
        for bad in [
            "input = \"m.obj\"\n",
            "input = \"m.obj\"\noutput = \"o\"\nunknown = 1\n",
            "input = \"m.obj\"\noutput = \"o\"\n[settings]\nelement_size = -1.0\n",
            "input = \"\"\noutput = \"o\"\n",
            "not toml",
        ] {
            assert!(
                matches!(
                    Project::from_toml(bad, Path::new("")),
                    Err(ProjectError::Invalid(_))
                ),
                "{bad}"
            );
        }
        assert!(matches!(
            Project::read(Path::new("does/not/exist.toml")),
            Err(ProjectError::Io(_))
        ));
    }
}
