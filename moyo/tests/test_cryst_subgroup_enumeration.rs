//! Optional cross-check of subgroup enumeration against GAP/Cryst.
//!
//! Cryst is licensed under GPL-2.0-or-later, while Moyo is licensed under
//! MIT OR Apache-2.0. Keep the two programs separate: this test neither links
//! to Cryst nor distributes Cryst code or data. It invokes a user-installed
//! `gap` executable as an independent process and compares only numeric
//! invariants printed on standard output.
//!
//! The test is ignored by default, so GAP/Cryst is not a build, test, or
//! distribution dependency of Moyo. Run it explicitly with:
//!
//! ```text
//! pixi run test-cryst
//! ```

use std::collections::BTreeMap;
use std::io::Write;
use std::process::{Command, Stdio};

use moyo::data::{Setting, operations_from_number};
use moyo::subgroup::{
    enumerate_klassengleiche_subgroups_by_index, enumerate_translationengleiche_subgroups,
};

const SPACE_GROUP_NUMBERS: std::ops::RangeInclusive<i32> = 1..=230;
const KLASSENGLEICHE_SPACE_GROUP_NUMBERS: std::ops::RangeInclusive<i32> = 1..=230;
const KLASSENGLEICHE_PRIMES: [usize; 2] = [2, 3];

#[derive(Debug, Default)]
struct CrystInvariants {
    translationengleiche_indices: BTreeMap<i32, Vec<usize>>,
    klassengleiche_counts: BTreeMap<(i32, usize), usize>,
}

#[test]
#[ignore = "requires a separately installed GAP with the GPL-licensed Cryst package"]
fn test_subgroup_enumeration_against_cryst() {
    let cryst = run_cryst();

    for number in SPACE_GROUP_NUMBERS {
        let operations = operations_from_number(number, Setting::Spglib, true).unwrap();
        let mut moyo_indices = enumerate_translationengleiche_subgroups(&operations, 1e-8)
            .unwrap()
            .into_iter()
            .filter(|class| class.representative.depth == 1)
            .map(|class| class.representative.translationengleiche_index)
            .collect::<Vec<_>>();
        moyo_indices.sort_unstable();

        assert_eq!(
            moyo_indices, cryst.translationengleiche_indices[&number],
            "maximal Translationengleiche subgroup indices differ for space group {number}"
        );
    }

    for number in KLASSENGLEICHE_SPACE_GROUP_NUMBERS {
        let operations = operations_from_number(number, Setting::Spglib, true).unwrap();
        for prime in KLASSENGLEICHE_PRIMES {
            let moyo_subgroups =
                enumerate_klassengleiche_subgroups_by_index(&operations, prime, 1e-8).unwrap();

            assert!(
                moyo_subgroups
                    .iter()
                    .all(|class| class.representative.klassengleiche_index == prime),
                "Moyo returned a Klassengleiche subgroup with the wrong index for space group {number} at index {prime}"
            );

            assert_eq!(
                moyo_subgroups.len(),
                cryst.klassengleiche_counts[&(number, prime)],
                "prime-index Klassengleiche subgroup counts differ for space group {number} at index {prime}"
            );
        }
    }
}

fn run_cryst() -> CrystInvariants {
    let mut child = Command::new("gap")
        .arg("-q")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to start `gap`; install GAP with the Cryst package");

    child
        .stdin
        .take()
        .unwrap()
        .write_all(cryst_program().as_bytes())
        .expect("failed to send the validation program to GAP");

    let output = child
        .wait_with_output()
        .expect("failed to wait for GAP/Cryst");
    assert!(
        output.status.success(),
        "GAP/Cryst failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    parse_cryst_output(&String::from_utf8(output.stdout).expect("GAP emitted non-UTF-8 output"))
}

fn cryst_program() -> String {
    let klassengleiche_numbers = KLASSENGLEICHE_SPACE_GROUP_NUMBERS
        .map(|number| number.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let klassengleiche_primes = KLASSENGLEICHE_PRIMES
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(",");

    format!(
        r#"
if LoadPackage("cryst") = fail then
    PrintTo("*errout*", "the GAP Cryst package is unavailable\n");
    QUIT_GAP(1);
fi;

for number in [1..230] do
    space_group := SpaceGroupIT(3, number);
    classes := MaximalSubgroupClassReps(
        space_group,
        rec(latticeequal := true)
    );
    indices := List(classes, IndexInParent);
    Sort(indices);
    Print("MOYO_CRYST|t|", number, "|", Length(indices));
    for index in indices do
        Print("|", index);
    od;
    Print("\n");
od;

for number in [{klassengleiche_numbers}] do
    space_group := SpaceGroupIT(3, number);
    for prime in [{klassengleiche_primes}] do
        classes := MaximalSubgroupClassReps(
            space_group,
            rec(classequal := true, primes := [prime])
        );
        count := Length(
            Filtered(classes, subgroup -> IndexInParent(subgroup) = prime)
        );
        Print("MOYO_CRYST|k|", number, "|", prime, "|", count, "\n");
    od;
od;

QUIT_GAP(0);
"#
    )
}

fn parse_cryst_output(output: &str) -> CrystInvariants {
    let mut invariants = CrystInvariants::default();

    for line in output.lines().map(|line| line.trim_start_matches('\r')) {
        if !line.starts_with("MOYO_CRYST|") {
            continue;
        }
        let fields = line.split('|').collect::<Vec<_>>();
        match fields[1] {
            "t" => {
                let number = parse_field(fields[2], line);
                let count = parse_field(fields[3], line);
                let indices = parse_indices(&fields[4..], count, line);
                assert!(
                    invariants
                        .translationengleiche_indices
                        .insert(number, indices)
                        .is_none(),
                    "duplicate Cryst result: {line}"
                );
            }
            "k" => {
                assert_eq!(fields.len(), 5, "invalid Cryst result: {line}");
                let number = parse_field(fields[2], line);
                let prime = parse_field(fields[3], line);
                let count = parse_field(fields[4], line);
                assert!(
                    invariants
                        .klassengleiche_counts
                        .insert((number, prime), count)
                        .is_none(),
                    "duplicate Cryst result: {line}"
                );
            }
            kind => panic!("unknown Cryst result kind `{kind}`: {line}"),
        }
    }

    assert_eq!(
        invariants.translationengleiche_indices.len(),
        SPACE_GROUP_NUMBERS.count(),
        "GAP/Cryst did not report every space group; output:\n{output}"
    );
    assert_eq!(
        invariants.klassengleiche_counts.len(),
        KLASSENGLEICHE_SPACE_GROUP_NUMBERS.count() * KLASSENGLEICHE_PRIMES.len(),
        "GAP/Cryst did not report every Klassengleiche case; output:\n{output}"
    );

    invariants
}

fn parse_field<T>(field: &str, line: &str) -> T
where
    T: std::str::FromStr,
    T::Err: std::fmt::Debug,
{
    field
        .parse()
        .unwrap_or_else(|error| panic!("invalid Cryst result `{line}`: {error:?}"))
}

fn parse_indices(fields: &[&str], count: usize, line: &str) -> Vec<usize> {
    assert_eq!(fields.len(), count, "invalid Cryst result: {line}");
    fields
        .iter()
        .map(|field| parse_field(field, line))
        .collect()
}
