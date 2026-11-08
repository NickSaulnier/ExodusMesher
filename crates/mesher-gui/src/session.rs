//! What the user has chosen: the loaded scene, which parts are selected, the
//! job settings and the output directory. It turns into a [`Project`], which
//! the GUI runs and saves, and which `mesher run project.toml` replays.

use std::path::PathBuf;

use mesher_core::JobConfig;
use mesher_core::project::Project;

use crate::scene::Scene;

/// File name of the project saved beside the outputs.
pub const PROJECT_FILE: &str = "project.toml";

#[derive(Debug, Clone)]
pub struct Session {
    pub scene: Scene,
    pub selected: Vec<bool>,
    pub settings: JobConfig,
    pub output: PathBuf,
}

impl Session {
    /// A fresh session: every part selected, default settings, output in
    /// `<stem>_out` beside the input.
    #[must_use]
    pub fn new(scene: Scene) -> Session {
        let stem = scene
            .path
            .file_stem()
            .map_or_else(|| "mesh".to_owned(), |s| s.to_string_lossy().into_owned());
        let output = scene
            .path
            .parent()
            .unwrap_or(std::path::Path::new(""))
            .join(format!("{stem}_out"));
        Session {
            selected: vec![true; scene.parts.len()],
            settings: JobConfig::default(),
            output,
            scene,
        }
    }

    /// Adopt a saved project for the same input: its selection, settings and
    /// output directory.
    pub fn apply(&mut self, project: &Project) {
        self.settings = project.settings.clone();
        self.output.clone_from(&project.output);
        let wanted = &project.settings.parts;
        for (i, (part, sel)) in self.scene.parts.iter().zip(&mut self.selected).enumerate() {
            *sel = wanted.is_empty()
                || wanted
                    .iter()
                    .any(|w| *w == part.name || *w == (i + 1).to_string());
        }
    }

    pub fn toggle(&mut self, part: usize) {
        if let Some(s) = self.selected.get_mut(part) {
            *s = !*s;
        }
    }

    #[must_use]
    pub fn any_selected(&self) -> bool {
        self.selected.iter().any(|&s| s)
    }

    /// The selection as `--parts` values: empty when every part is selected;
    /// names when they are unique, 1-based indices otherwise.
    #[must_use]
    pub fn part_args(&self) -> Vec<String> {
        if self.selected.iter().all(|&s| s) {
            return Vec::new();
        }
        let parts = &self.scene.parts;
        self.selected
            .iter()
            .enumerate()
            .filter(|(_, s)| **s)
            .map(|(i, _)| {
                let name = &parts[i].name;
                if parts.iter().filter(|p| p.name == *name).count() == 1 {
                    name.clone()
                } else {
                    (i + 1).to_string()
                }
            })
            .collect()
    }

    /// The job as a project.
    #[must_use]
    pub fn project(&self) -> Project {
        Project {
            input: self.scene.path.clone(),
            output: self.output.clone(),
            settings: JobConfig {
                parts: self.part_args(),
                ..self.settings.clone()
            },
        }
    }

    /// Where the project is saved: beside the outputs.
    #[must_use]
    pub fn project_path(&self) -> PathBuf {
        self.output.join(PROJECT_FILE)
    }
}
