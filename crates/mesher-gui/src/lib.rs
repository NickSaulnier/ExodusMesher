//! The ExodusMesher GUI's logic, kept apart from rendering so it can be tested
//! headless: the orbit camera, the loaded scene (read and diagnosed exactly as
//! the worker does), and the session that becomes a project file. The GUI
//! holds no meshing logic; it runs the worker through `mesher-client`.

// Camera and display maths only: nothing here reaches a mesh or a report, so
// the platform's maths library is fine (see `clippy.toml`).
#![allow(clippy::disallowed_methods)]

pub mod camera;
pub mod scene;
pub mod session;
