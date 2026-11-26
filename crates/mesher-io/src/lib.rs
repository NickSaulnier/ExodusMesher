//! Readers and writers: the OBJ reader and writer, the defect overlay (.vtp),
//! and the volume-mesh writers: Exodus II, Gmsh .msh 4.1, VTK .vtu, an
//! OpenFOAM polyMesh case, a Fluent mesh file and libMeshb (`.meshb`, also
//! read as an input surface).

pub mod exodus;
pub mod fluent;
pub mod foam;
pub mod meshb;
pub mod msh;
pub mod obj;
mod shared;
pub mod vtp;
pub mod vtu;

pub use mesher_core::model::{Polygon, SurfaceModel, Triangle};
pub use obj::{NamedSurface, ObjModel, parse_obj, read_obj, write_obj};

use std::path::Path;

use mesher_core::{Action, DefectCode, Outcome, Record, RecordLog, Stage};

/// Whether `path` is read as a libMeshb binary file (extension `.meshb`,
/// any case) rather than as OBJ.
#[must_use]
pub fn is_meshb(path: &Path) -> bool {
    path.extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("meshb"))
}

/// Read an input surface: a `.meshb` file by its extension, anything else as
/// OBJ. On an I/O error, records `PARSE_READ_FAILED` and returns `None`.
pub fn read_surface(
    path: &Path,
    log: &mut RecordLog,
    progress: &mut dyn FnMut(f64),
) -> Option<SurfaceModel> {
    if !is_meshb(path) {
        return read_obj(path, log, progress);
    }
    match std::fs::read(path) {
        Ok(bytes) => {
            let model = meshb::parse_meshb(&bytes, log);
            progress(1.0);
            Some(model)
        }
        Err(err) => {
            log.push(
                Record::new(DefectCode::ParseReadFailed, Stage::Parse)
                    .measured("io_error", err.kind().to_string())
                    .action(Action::Abort, Outcome::Applied)
                    .hint(format!("could not read {}: {err}", path.display())),
            );
            None
        }
    }
}
