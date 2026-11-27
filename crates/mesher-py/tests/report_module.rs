//! Generates the Python package's report types (`python/mesher/_schema.py`)
//! from the committed report schema, and checks the committed module is
//! current: run `UPDATE_REPORT_MODULE=1 cargo test -p mesher-py --test
//! report_module` after the schema changes.
//!
//! Each object in the schema becomes a frozen dataclass with a `_parse`
//! classmethod that builds it from the report's parsed JSON, and each set
//! of string constants a `str` enum. A field the worker may leave out is
//! its schema default when it has one, else empty for a list or a map, else
//! `None`.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::fmt::Write as _;
use std::path::PathBuf;

use serde_json::{Map, Value};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

#[test]
fn generated_report_module_is_current() {
    let schema_text = std::fs::read_to_string(repo_root().join("schema/report.schema.json"))
        .unwrap()
        .replace("\r\n", "\n");
    let schema: Value = serde_json::from_str(&schema_text).unwrap();
    let generated = generate(&schema);
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("python/mesher/_schema.py");
    if std::env::var_os("UPDATE_REPORT_MODULE").is_some() {
        std::fs::write(&path, &generated).unwrap();
        return;
    }
    let committed = std::fs::read_to_string(&path)
        .unwrap_or_default()
        .replace("\r\n", "\n");
    assert!(
        committed == generated,
        "crates/mesher-py/python/mesher/_schema.py is stale; run \
         `UPDATE_REPORT_MODULE=1 cargo test -p mesher-py --test report_module`"
    );
}

/// Every object and enum in the schema is in the module, and every
/// property of each object is a field of its class.
#[test]
fn every_schema_type_and_property_is_generated() {
    let text = std::fs::read_to_string(repo_root().join("schema/report.schema.json")).unwrap();
    let schema: Value = serde_json::from_str(&text).unwrap();
    let module = generate(&schema);
    let defs = schema["$defs"].as_object().unwrap();
    for (name, def) in defs.iter().chain([(&"Report".to_owned(), &schema)]) {
        assert!(
            module.contains(&format!("\nclass {name}:"))
                || module.contains(&format!("\nclass {name}(str, Enum):")),
            "{name} missing"
        );
        for prop in def
            .get("properties")
            .and_then(Value::as_object)
            .into_iter()
            .flat_map(Map::keys)
        {
            assert!(
                module.contains(&format!("\n    {}: ", field_name(prop))),
                "{name}.{prop} missing"
            );
        }
    }
}

const KEYWORDS: &[&str] = &[
    "False", "None", "True", "and", "as", "assert", "async", "await", "break", "class", "continue",
    "def", "del", "elif", "else", "except", "finally", "for", "from", "global", "if", "import",
    "in", "is", "lambda", "nonlocal", "not", "or", "pass", "raise", "return", "try", "while",
    "with", "yield",
];

fn field_name(prop: &str) -> String {
    if KEYWORDS.contains(&prop) {
        format!("{prop}_")
    } else {
        prop.to_owned()
    }
}

/// A docstring's body: backslashes and triple quotes escaped, lines kept.
fn doc(text: &str, indent: &str) -> String {
    let body = text
        .trim()
        .replace('\\', "\\\\")
        .replace("\"\"\"", "\\\"\\\"\\\"");
    let mut lines = body.lines();
    let mut out = format!("{indent}\"\"\"{}", lines.next().unwrap_or(""));
    let rest: Vec<&str> = lines.collect();
    if rest.is_empty() {
        out.push_str("\"\"\"\n");
    } else {
        out.push('\n');
        for line in rest {
            if line.is_empty() {
                out.push('\n');
            } else {
                let _ = writeln!(out, "{indent}{line}");
            }
        }
        let _ = writeln!(out, "{indent}\"\"\"");
    }
    out
}

/// The enum's values, each with its description, if the def is a set of
/// string constants.
fn enum_values(def: &Value) -> Option<Vec<(String, Option<String>)>> {
    if let Some(values) = def.get("enum").and_then(Value::as_array) {
        return Some(
            values
                .iter()
                .map(|v| (v.as_str().unwrap().to_owned(), None))
                .collect(),
        );
    }
    let variants = def.get("oneOf")?.as_array()?;
    let mut out = Vec::new();
    for variant in variants {
        let description = variant
            .get("description")
            .and_then(Value::as_str)
            .map(str::to_owned);
        if let Some(c) = variant.get("const") {
            out.push((c.as_str()?.to_owned(), description));
            continue;
        }
        for v in variant.get("enum").and_then(Value::as_array)? {
            out.push((v.as_str()?.to_owned(), description.clone()));
        }
    }
    Some(out)
}

fn ref_name(v: &Value) -> Option<&str> {
    v.get("$ref")?.as_str()?.strip_prefix("#/$defs/")
}

/// The Python expression converting the JSON value `x` (an expression),
/// at a nesting depth (for comprehension variables' names).
type Convert = Box<dyn Fn(&str, usize) -> String>;

struct Ty {
    /// The annotation.
    annotation: String,
    /// The conversion of the JSON value to it.
    convert: Convert,
    /// The value for a field the worker left out, absent a schema default.
    missing: &'static str,
}

fn ty(v: &Value, defs: &Map<String, Value>) -> Ty {
    if let Some(name) = ref_name(v) {
        let name = name.to_owned();
        let is_enum = enum_values(&defs[&name]).is_some();
        let annotation = name.clone();
        return Ty {
            annotation,
            convert: Box::new(move |x, _| {
                if is_enum {
                    format!("{name}({x})")
                } else {
                    format!("{name}._parse({x})")
                }
            }),
            missing: "None",
        };
    }
    if let Some(any) = v.get("anyOf").and_then(Value::as_array) {
        let inner: Vec<&Value> = any.iter().filter(|a| a["type"] != "null").collect();
        assert!(inner.len() == 1 && any.len() == 2, "unsupported anyOf: {v}");
        return nullable(ty(inner[0], defs));
    }
    let types: Vec<&str> = match &v["type"] {
        Value::String(t) => vec![t.as_str()],
        Value::Array(ts) => ts.iter().map(|t| t.as_str().unwrap()).collect(),
        other => panic!("unsupported type {other}: {v}"),
    };
    let nullable_type = types.contains(&"null");
    let base: Vec<&str> = types.into_iter().filter(|&t| t != "null").collect();
    assert!(base.len() == 1, "unsupported type: {v}");
    let t = match base[0] {
        "string" => scalar("str", false),
        "integer" => scalar("int", false),
        "number" => scalar("float", true),
        "boolean" => scalar("bool", false),
        "array" => {
            let item = ty(&v["items"], defs);
            let n = v.get("minItems").and_then(Value::as_u64);
            let annotation = match n {
                Some(n) if v.get("maxItems").and_then(Value::as_u64) == Some(n) => {
                    format!(
                        "Tuple[{}]",
                        vec![item.annotation.as_str(); n as usize].join(", ")
                    )
                }
                _ => format!("Tuple[{}, ...]", item.annotation),
            };
            let convert = item.convert;
            Ty {
                annotation,
                convert: Box::new(move |x, depth| {
                    let e = format!("e{depth}");
                    let item = convert(&e, depth + 1);
                    if item == e {
                        format!("tuple({x})")
                    } else {
                        format!("tuple({item} for {e} in {x})")
                    }
                }),
                missing: "()",
            }
        }
        "object" => match v.get("additionalProperties") {
            Some(Value::Bool(true)) | None => Ty {
                annotation: "Dict[str, Any]".to_owned(),
                convert: Box::new(|x, _| format!("dict({x})")),
                missing: "{}",
            },
            Some(item) => {
                let item = ty(item, defs);
                let annotation = format!("Dict[str, {}]", item.annotation);
                let convert = item.convert;
                Ty {
                    annotation,
                    convert: Box::new(move |x, depth| {
                        let (k, e) = (format!("k{depth}"), format!("e{depth}"));
                        format!(
                            "{{{k}: {} for {k}, {e} in {x}.items()}}",
                            convert(&e, depth + 1)
                        )
                    }),
                    missing: "{}",
                }
            }
        },
        other => panic!("unsupported type {other}: {v}"),
    };
    if nullable_type { nullable(t) } else { t }
}

fn scalar(name: &str, float: bool) -> Ty {
    Ty {
        annotation: name.to_owned(),
        convert: Box::new(move |x, _| {
            if float {
                format!("float({x})")
            } else {
                x.to_owned()
            }
        }),
        missing: "None",
    }
}

fn nullable(t: Ty) -> Ty {
    let convert = t.convert;
    Ty {
        annotation: format!("Optional[{}]", t.annotation),
        convert: Box::new(move |x, depth| {
            let inner = convert(x, depth);
            if inner == x {
                x.to_owned()
            } else {
                format!("None if {x} is None else {inner}")
            }
        }),
        missing: "None",
    }
}

/// A Python literal for a schema default.
fn literal(v: &Value) -> String {
    match v {
        Value::Null => "None".to_owned(),
        Value::Bool(b) => if *b { "True" } else { "False" }.to_owned(),
        Value::Number(n) => n.to_string(),
        Value::String(s) => format!("{s:?}"),
        Value::Array(a) if a.is_empty() => "()".to_owned(),
        other => panic!("unsupported default {other}"),
    }
}

fn class(name: &str, def: &Value, defs: &Map<String, Value>) -> String {
    let mut out = String::new();
    let props = def["properties"].as_object().unwrap();
    let required: Vec<&str> = def
        .get("required")
        .and_then(Value::as_array)
        .map(|r| r.iter().map(|x| x.as_str().unwrap()).collect())
        .unwrap_or_default();
    let _ = writeln!(out, "@dataclass(frozen=True)\nclass {name}:");
    out.push_str(&doc(
        def.get("description")
            .and_then(Value::as_str)
            .unwrap_or(name),
        "    ",
    ));
    // Required fields first: dataclasses put fields with defaults last.
    let mut order: Vec<(&String, &Value)> = props.iter().collect();
    order.sort_by_key(|(p, _)| !required.contains(&p.as_str()));
    let mut parse = Vec::new();
    for (prop, schema) in order {
        let t = ty(schema, defs);
        let field = field_name(prop);
        let access = format!("v[{prop:?}]");
        out.push('\n');
        if required.contains(&prop.as_str()) {
            let _ = writeln!(out, "    {field}: {}", t.annotation);
            parse.push(format!("{field}={}", (t.convert)(&access, 0)));
        } else {
            let default = schema.get("default").map_or(t.missing.to_owned(), literal);
            let declared = if default == "{}" {
                "field(default_factory=dict)".to_owned()
            } else {
                default.clone()
            };
            let annotation = if default == "None" && !t.annotation.starts_with("Optional[") {
                format!("Optional[{}]", t.annotation)
            } else {
                t.annotation.clone()
            };
            let _ = writeln!(out, "    {field}: {annotation} = {declared}");
            let converted = (t.convert)(&access, 0);
            parse.push(if converted == access {
                format!("{field}=v.get({prop:?}, {default})")
            } else {
                format!("{field}=({converted}) if {prop:?} in v else {default}")
            });
        }
        if let Some(d) = schema.get("description").and_then(Value::as_str) {
            out.push_str(&doc(d, "    "));
        }
    }
    let _ = write!(
        out,
        "\n    @classmethod\n    def _parse(cls, v: Mapping[str, Any]) -> {name}:\n        return cls(\n"
    );
    for p in parse {
        let _ = writeln!(out, "            {p},");
    }
    out.push_str("        )\n");
    out
}

fn enumeration(name: &str, def: &Value, values: &[(String, Option<String>)]) -> String {
    let mut out = format!("class {name}(str, Enum):\n");
    out.push_str(&doc(
        def.get("description")
            .and_then(Value::as_str)
            .unwrap_or(name),
        "    ",
    ));
    out.push('\n');
    for (value, description) in values {
        if let Some(d) = description {
            for line in d.trim().lines() {
                let _ = writeln!(out, "{}", format!("    #: {line}").trim_end());
            }
        }
        let member = value.to_uppercase().replace(['-', '.', ' '], "_");
        let _ = writeln!(out, "    {member} = {value:?}");
    }
    out
}

fn generate(schema: &Value) -> String {
    let defs = schema["$defs"].as_object().unwrap();
    let mut out = String::from(
        r#""""The report's types, generated from ``schema/report.schema.json`` by
``crates/mesher-py/tests/report_module.rs``: do not edit; run
``UPDATE_REPORT_MODULE=1 cargo test -p mesher-py --test report_module``.

Each object in the schema is a frozen dataclass (built from the report's
JSON by its ``_parse``), each set of string constants a ``str`` enum. A field
the worker may leave out is its schema default when it has one, else empty
for a list or a map, else ``None``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, Mapping, Optional, Tuple

"#,
    );
    let mut names: Vec<&str> = defs.keys().map(String::as_str).collect();
    names.push("Report");
    names.sort_unstable();
    let _ = writeln!(out, "__all__ = [");
    for n in &names {
        let _ = writeln!(out, "    {n:?},");
    }
    out.push_str("]\n");
    for (name, def) in defs {
        out.push_str("\n\n");
        match enum_values(def) {
            Some(values) => out.push_str(&enumeration(name, def, &values)),
            None => out.push_str(&class(name, def, defs)),
        }
    }
    out.push_str("\n\n");
    out.push_str(&class("Report", schema, defs));
    out
}
