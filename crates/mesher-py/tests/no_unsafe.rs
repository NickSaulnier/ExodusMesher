//! The crate allows `unsafe_code` in its Python module only because PyO3's
//! macros expand to it; nothing in the source may use `unsafe` itself.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

#[test]
fn no_unsafe_is_written_by_hand() {
    let src = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let mut checked = 0;
    for entry in std::fs::read_dir(&src).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().is_some_and(|e| e == "rs") {
            let text = std::fs::read_to_string(&path).unwrap();
            for (i, line) in text.lines().enumerate() {
                let code = line.split("//").next().unwrap_or("");
                assert!(
                    !code
                        .split(|c: char| !c.is_alphanumeric() && c != '_')
                        .any(|w| w == "unsafe"),
                    "{}:{}: `unsafe` in the source",
                    path.display(),
                    i + 1
                );
            }
            checked += 1;
        }
    }
    assert!(checked >= 2, "no source files found in {}", src.display());
}
