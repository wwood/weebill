use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::predicate;
use serial_test::serial;
use std::fs;
use std::path::Path;
use std::process::Command;
use std::str; // Run programs

fn fresh() {
    let dir = "./tests/results/test_sketch_dir";
    if Path::new(dir).exists() {
        let _ = fs::remove_dir_all(dir);
    }
}

#[serial]
#[test]
fn test_sketch_commands() {
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("test_files/e.coli-EC590.fasta.gz")
        .arg("test_files/e.coli-K12.fasta.gz")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-o")
        .arg("tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("-l")
        .arg("./test_files/list.txt")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("-i")
        .arg("-m")
        .arg("90")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("./test_files/t1.fq")
        .arg("-2")
        .arg("./test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);
    assert!(
        Path::new("./tests/results/test_sketch_dir/t1.fq.paired.sylsp").exists(),
        "Output file was not created"
    );
    fresh();

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("--l1")
        .arg("./test_files/pair_list1.txt")
        .arg("--l2")
        .arg("./test_files/pair_list2.txt")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);
    assert!(
        Path::new("./tests/results/test_sketch_dir/t1.fq.paired.sylsp").exists(),
        "Output file was not created"
    );

    fresh();
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-g")
        .arg("./test_files/t1.fq")
        .arg("-r")
        .arg("./test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/testdb")
        .assert();
    assert.success().code(0);
    assert!(
        Path::new("./tests/results/test_sketch_dir/t2.fq.sylsp").exists(),
        "Output file was not created"
    );
    assert!(
        Path::new("./tests/results/test_sketch_dir/testdb.syldb").exists(),
        "Output file was not created"
    );
}

#[serial]
#[test]
fn test_profile_vs_query() {
    fresh();

    let mut output = Command::cargo_bin("weebill").unwrap();
    let output = output
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .output()
        .expect("Output failed");
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    dbg!(stdout.matches('\n').count());
    assert!(stdout.matches('\n').count() == 2);

    let mut output = Command::cargo_bin("weebill").unwrap();
    let output = output
        .arg("query")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .output()
        .expect("Output failed");
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    dbg!(stdout.matches('\n').count());
    println!("{}", stdout);
    assert!(stdout.matches('\n').count() == 4);
}

#[serial]
#[test]
fn test_sketch_list() {
    fresh();
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-r")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);
    assert!(
        Path::new("./tests/results/test_sketch_dir/e.coli-EC590.fasta.gz.sylsp").exists(),
        "Output file was not created"
    );
    assert!(
        Path::new("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp").exists(),
        "Output file was not created"
    );
    assert!(
        !Path::new("./tests/results/test_sketch_dir/db.syldb").exists(),
        "Output file was created"
    );
    fresh();

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-g")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);
    assert!(
        !Path::new("./tests/results/test_sketch_dir/e.coli-EC590.fasta.gz.sylsp").exists(),
        "Output file was created"
    );
    assert!(
        !Path::new("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp").exists(),
        "Output file was created"
    );
    assert!(
        Path::new("./tests/results/test_sketch_dir/db.syldb").exists(),
        "Output file was not created"
    );
    fresh();

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("--gl")
        .arg("test_files/list.txt")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/db")
        .assert();
    assert.success().code(0);
    assert!(
        Path::new("./tests/results/test_sketch_dir/db.syldb").exists(),
        "Output file was not created"
    );
    fresh();

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("--rl")
        .arg("test_files/list.txt")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);
    assert!(
        !Path::new("./tests/results/test_sketch_dir/db.syldb").exists(),
        "Output file was not created"
    );
    assert!(
        Path::new("./tests/results/test_sketch_dir/e.coli-EC590.fasta.gz.sylsp").exists(),
        "Output file was not created"
    );
    assert!(
        Path::new("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp").exists(),
        "Output file was not created"
    );
    fresh();
}
#[serial]
#[test]
fn test_profile_disabling() {
    fresh();

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-g")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("--disable-profiling")
        .assert();
    assert.success().code(0);

    let mut output = Command::cargo_bin("weebill").unwrap();
    let assert = output
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .assert();
    assert.failure().code(1);

    let mut output = Command::cargo_bin("weebill").unwrap();
    let assert = output
        .arg("query")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .assert();
    assert.success().code(0);

    fresh();
}
#[serial]
#[test]
fn test_sketch_fasta_fastq_concord() {
    fresh();
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);

    let mut output = Command::cargo_bin("weebill").unwrap();
    let out1 = output
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .output()
        .expect("Fail");

    let mut output = Command::cargo_bin("weebill").unwrap();
    let out2 = output
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .output()
        .expect("Fail");

    let mut output = Command::cargo_bin("weebill").unwrap();
    let out3 = output
        .arg("profile")
        .arg("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .output()
        .expect("Fail");

    let stdout1 = str::from_utf8(&out1.stdout).expect("Output was not valid UTF-8");
    let stdout2 = str::from_utf8(&out2.stdout).expect("Output was not valid UTF-8");
    let stdout3 = str::from_utf8(&out3.stdout).expect("Output was not valid UTF-8");

    assert!(stdout1 == stdout2);
    assert!(stdout1 == stdout3);
    assert!(stdout2 == stdout3);

    fresh();
}
#[serial]
#[test]
fn test_sample_names() {
    fresh();
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("--lS")
        .arg("./test_files/single_sample.txt")
        .assert();
    assert.success().code(0);
    assert!(
        Path::new("./tests/results/test_sketch_dir/SAMPLE_TEST.paired.sylsp").exists(),
        "Output file was not created"
    );
    fresh();

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("test_files/t1.fq")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("--lS")
        .arg("./test_files/sample_list.txt")
        .assert();
    assert.success().code(0);
    assert!(
        Path::new("./tests/results/test_sketch_dir/S1.sylsp").exists(),
        "Output file was not created"
    );
    assert!(
        Path::new("./tests/results/test_sketch_dir/S2.sylsp").exists(),
        "Output file was not created"
    );

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let output = cmd
        .arg("profile")
        .arg("./tests/results/test_sketch_dir/S2.sylsp")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .output()
        .unwrap();
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    dbg!(&stdout);
    assert!(stdout.contains("S2"));
    assert!(!stdout.contains("o157_reads"));

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("-S")
        .arg("SAMPLE_TEST_S")
        .assert();
    assert.success().code(0);
    assert!(
        Path::new("./tests/results/test_sketch_dir/SAMPLE_TEST_S.paired.sylsp").exists(),
        "Output file was not created, -S"
    );

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .arg("test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("-S")
        .arg("SAMPLE_TEST_S")
        .arg("SAMPLE_TEST_S1")
        .assert();
    assert.success().code(0);
    assert!(
        Path::new("./tests/results/test_sketch_dir/SAMPLE_TEST_S1.paired.sylsp").exists(),
        "Output file was not created, -S"
    );

    fresh();
}
#[serial]
#[test]
fn test_fpr() {
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("0")
        .assert();
    assert.success().code(0);
    fresh();

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("--fpr")
        .arg("0.001")
        .assert();
    assert.success().code(0);
    fresh();
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("--fpr")
        .arg("2")
        .assert();
    assert.failure().code(1);
    fresh();
}
#[serial]
#[test]
fn test_raw_inputs_profile_simple() {
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .assert();
    assert.success().code(0);
    fresh();

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("-1")
        .arg("test_files/t1.fq")
        .assert();
    assert.failure().code(1);
    fresh();

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .arg("test_files/t1.fq")
        .assert();
    assert.success().code(0);
    fresh();
}

#[serial]
#[test]
fn test_estimate_read_counts() {
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let output = cmd
        .arg("profile")
        .arg("--estimate-read-counts")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .output()
        .expect("output failed");
    let stdout_1 = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    let mut lines = stdout_1.lines();
    dbg!(stdout_1);
    lines.next();
    let output = lines.next().unwrap();
    let split: Vec<&str> = output.split('\t').collect();
    assert!(split[3].parse::<f64>().unwrap() > 1000.0);

    fresh();

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let output = cmd
        .arg("profile")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .output()
        .expect("output failed");
    let stdout_1 = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    let mut lines = stdout_1.lines();
    dbg!(stdout_1);
    lines.next();
    let output = lines.next().unwrap();
    let split: Vec<&str> = output.split('\t').collect();
    assert!(split[3].parse::<f64>().unwrap() < 101.00);

    fresh();
}

#[serial]
#[test]
fn test_raw_inputs_profile_with_sketch() {
    let mut output = Command::cargo_bin("weebill").unwrap();
    let output = output
        .arg("profile")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .output()
        .expect("Output failed");
    let stdout_1 = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);

    let mut output = Command::cargo_bin("weebill").unwrap();
    let output = output
        .arg("profile")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./tests/results/test_sketch_dir/k12_R1.fq.paired.sylsp")
        .output()
        .expect("Output failed");
    let stdout_2 = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");

    assert!(stdout_1 == stdout_2);
}

#[serial]
#[test]
fn test_inspect() {
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("test_files/e.coli-EC590.fasta.gz")
        .arg("test_files/e.coli-K12.fasta.gz")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-o")
        .arg("tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let output = cmd
        .arg("inspect")
        .arg("./tests/results/test_sketch_dir/k12_R1.fq.paired.sylsp")
        .output()
        .expect("Output failed");

    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    assert!(stdout.contains("k12_R1.fq"));

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let output = cmd
        .arg("inspect")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .output()
        .expect("Output failed");
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    assert!(stdout.contains("e.coli-EC590.fasta.gz"));
    assert!(stdout.contains("e.coli-K12.fasta.gz"));
}

#[serial]
#[test]
fn test_refdelta_query_with_reference() {
    fresh();
    let dir = "./tests/results/test_sketch_dir";

    // sketch a database and a sample
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("test_files/e.coli-K12.fasta.gz")
        .arg("test_files/e.coli-o157.fasta.gz")
        .arg("test_files/e.coli-EC590.fasta.gz")
        .arg("-o")
        .arg(format!("{}/db", dir))
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);

    // build a reference and compress the sample sketch against it
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("ref-build")
        .arg(format!("{}/db.syldb", dir))
        .arg("-o")
        .arg(format!("{}/ref", dir))
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("ref-compress")
        .arg(format!("{}/o157_reads.fastq.gz.sylsp", dir))
        .arg("-r")
        .arg(format!("{}/ref.sylref", dir))
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);
    assert!(
        Path::new(&format!("{}/o157_reads.fastq.gz.sylspr", dir)).exists(),
        "ref-compress did not produce a .sylspr"
    );

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let inspect = cmd
        .arg("ref-compress")
        .arg("--inspect")
        .arg(format!("{}/o157_reads.fastq.gz.sylspr", dir))
        .output()
        .expect("Output failed");
    assert!(inspect.status.success());
    let inspect_stdout = str::from_utf8(&inspect.stdout).expect("Output was not valid UTF-8");
    assert!(inspect_stdout.contains("reference_db"));
    assert!(inspect_stdout.contains("assigned_to_genomes"));
    assert!(inspect_stdout.contains(&format!("{}/ref.sylref", dir)));

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("ref-compress")
        .arg("--verify")
        .arg("-r")
        .arg(format!("{}/ref.sylref", dir))
        .arg(format!("{}/o157_reads.fastq.gz.sylspr", dir))
        .assert()
        .success()
        .code(0);

    // querying the .sylspr via --reference must match querying the original .sylsp
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let orig = cmd
        .arg("query")
        .arg(format!("{}/db.syldb", dir))
        .arg(format!("{}/o157_reads.fastq.gz.sylsp", dir))
        .output()
        .expect("Output failed");
    let orig = str::from_utf8(&orig.stdout).expect("Output was not valid UTF-8");

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let from_ref = cmd
        .arg("query")
        .arg(format!("{}/db.syldb", dir))
        .arg(format!("{}/o157_reads.fastq.gz.sylspr", dir))
        .arg("--reference")
        .arg(format!("{}/ref.sylref", dir))
        .output()
        .expect("Output failed");
    let from_ref = str::from_utf8(&from_ref.stdout).expect("Output was not valid UTF-8");

    assert!(orig.contains("e.coli-o157.fasta.gz"));
    assert_eq!(
        orig, from_ref,
        "query of .sylspr via --reference differs from query of original .sylsp"
    );

    // a .sylspr without --reference must fail
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("query")
        .arg(format!("{}/db.syldb", dir))
        .arg(format!("{}/o157_reads.fastq.gz.sylspr", dir))
        .assert()
        .failure();
}

/// `sketch --reference` accepts the same compression tunables as `ref-compress`, and they
/// change only how tightly the sample is packed: --no-error-kmer (or an unreachable
/// --min-dense-kmers-for-error) suppresses error-k-mer recoding (a bigger file), a near-100
/// --ref-screen-ani decodes fewer reference genomes (bigger still), and every setting still
/// round-trips losslessly.
#[serial]
#[test]
fn test_sketch_reference_compression_options() {
    fresh();
    let dir = "./tests/results/test_sketch_dir";

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("test_files/e.coli-K12.fasta.gz")
        .arg("test_files/e.coli-o157.fasta.gz")
        .arg("-o")
        .arg(format!("{}/db", dir))
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);
    // --store-genomes: error-k-mer encoding (and so --min-dense-kmers-for-error) only has
    // anything to do when the reference carries the genome sequences.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("ref-build")
        .arg(format!("{}/db.syldb", dir))
        .arg("--store-genomes")
        .arg("-o")
        .arg(format!("{}/ref", dir))
        .assert()
        .success()
        .code(0);
    let refdb = format!("{}/ref.sylref", dir);

    // sketch the reads straight to *.sylspr under three settings of the new flags
    let sketch_with = |subdir: &str, flags: &[&str]| -> u64 {
        let out = format!("{}/{}", dir, subdir);
        let mut cmd = Command::cargo_bin("weebill").unwrap();
        cmd.arg("sketch")
            .arg("test_files/o157_reads.fastq.gz")
            .arg("--reference")
            .arg(&refdb)
            .arg("-d")
            .arg(&out)
            .args(flags)
            .assert()
            .success()
            .code(0);
        let sample = format!("{}/o157_reads.fastq.gz.sylspr", out);
        fs::metadata(&sample)
            .unwrap_or_else(|_| panic!("sketch --reference did not produce {}", sample))
            .len()
    };
    let default_len = sketch_with("default", &[]);
    let no_error_len = sketch_with("no_error", &["--min-dense-kmers-for-error", "100000000"]);
    let strict_screen_len = sketch_with("strict_screen", &["--ref-screen-ani", "99.9"]);
    let no_error_kmer_len = sketch_with("no_error_kmer", &["--no-error-kmer"]);

    // --no-error-kmer and an unreachable --min-dense-kmers-for-error both suppress error-k-mer
    // recoding, by different routes, so they must land on the same bytes.
    assert_eq!(
        no_error_kmer_len, no_error_len,
        "--no-error-kmer gave {} bytes, but suppressing error k-mers via the dense threshold gave {}",
        no_error_kmer_len, no_error_len
    );
    assert!(
        no_error_len > default_len,
        "--min-dense-kmers-for-error was ignored: suppressing error-k-mer recoding gave {} bytes, same or less than the default {}",
        no_error_len,
        default_len
    );
    assert!(
        strict_screen_len > default_len,
        "--ref-screen-ani was ignored: screening out reference genomes gave {} bytes, same or less than the default {}",
        strict_screen_len,
        default_len
    );

    // whatever the settings, the sketch still decodes to the same thing
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let orig = cmd
        .arg("query")
        .arg(format!("{}/db.syldb", dir))
        .arg(format!("{}/o157_reads.fastq.gz.sylsp", dir))
        .output()
        .expect("Output failed");
    let orig = str::from_utf8(&orig.stdout).expect("Output was not valid UTF-8");
    assert!(orig.contains("e.coli-o157.fasta.gz"));

    for subdir in ["default", "no_error", "strict_screen", "no_error_kmer"] {
        let mut cmd = Command::cargo_bin("weebill").unwrap();
        let from_ref = cmd
            .arg("query")
            .arg(format!("{}/db.syldb", dir))
            .arg(format!("{}/{}/o157_reads.fastq.gz.sylspr", dir, subdir))
            .arg("--reference")
            .arg(&refdb)
            .output()
            .expect("Output failed");
        let from_ref = str::from_utf8(&from_ref.stdout).expect("Output was not valid UTF-8");
        // the sample path differs between the .sylsp and .sylspr runs; the rest must not
        let strip = |s: &str| {
            s.replace(&format!("{}/", subdir), "")
                .replace(".sylspr", "")
                .replace(".sylsp", "")
        };
        assert_eq!(
            strip(orig),
            strip(from_ref),
            "sketch --reference ({}) did not round-trip to the same query result",
            subdir
        );
    }
}

#[serial]
#[test]
fn test_two_stage_profile() {
    fresh();
    let dir = "./tests/results/two_stage";
    let _ = fs::remove_dir_all(dir);

    // Sparse (-c 200) database that retains the source fasta paths.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-c")
        .arg("200")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .arg("-o")
        .arg(format!("{}/db_c200", dir))
        .assert()
        .success()
        .code(0);

    // Dense (-c 50) read sample.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-c")
        .arg("50")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);

    let db = format!("{}/db_c200.syldb", dir);
    let sample = format!("{}/o157_reads.fastq.gz.sylsp", dir);
    let cache = format!("{}/cache", dir);

    // Two-stage profile: screen at c=200, densely profile the survivors at c=50
    // by re-sketching their source fastas, caching the dense sketches.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let output = cmd
        .arg("profile")
        .arg("--two-stage")
        .arg("--dense-c")
        .arg("50")
        .arg("--dense-cache")
        .arg(&cache)
        .arg(&db)
        .arg(&sample)
        .output()
        .expect("Output failed");
    assert!(output.status.success());
    let two_stage = str::from_utf8(&output.stdout)
        .expect("not UTF-8")
        .to_string();
    // The reads are E. coli O157 -> that genome must be profiled.
    assert!(two_stage.contains("e.coli-o157.fasta.gz"));
    // The dense stage caches a per-genome sketch for each screened survivor.
    assert!(Path::new(&cache).exists());
    assert!(
        fs::read_dir(&cache).unwrap().count() >= 1,
        "dense cache was not populated"
    );

    // Genomes detected by two-stage must equal those of a plain single-stage
    // profile of the same dense reads against a dense (-c 50) database.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-c")
        .arg("50")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .arg("-o")
        .arg(format!("{}/db_c50", dir))
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let output = cmd
        .arg("profile")
        .arg(format!("{}/db_c50.syldb", dir))
        .arg(&sample)
        .output()
        .expect("Output failed");
    let single = str::from_utf8(&output.stdout)
        .expect("not UTF-8")
        .to_string();

    let detected = |tsv: &str| -> Vec<String> {
        let mut v: Vec<String> = tsv
            .lines()
            .skip(1)
            .filter_map(|l| l.split('\t').nth(1).map(|s| s.to_string()))
            .collect();
        v.sort();
        v
    };
    assert_eq!(
        detected(&two_stage),
        detected(&single),
        "two-stage and single-stage detected different genome sets"
    );
}

#[serial]
#[test]
fn test_two_stage_db_convert_and_profile() {
    fresh();
    let dir = "./tests/results/two_stage_db";
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();

    // Dense (-c 50) database carrying profiling k-mers.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-c")
        .arg("50")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .arg("-o")
        .arg(format!("{}/db_c50", dir))
        .assert()
        .success()
        .code(0);

    // Dense (-c 50) read sample.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-c")
        .arg("50")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);

    let dense_db = format!("{}/db_c50.syldb", dir);
    let sample = format!("{}/o157_reads.fastq.gz.sylsp", dir);

    // Convert the dense db into a two-stage seekable database: dense blocks at
    // c=50, sparse stage-1 screen index at c=200.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("db-convert")
        .arg(&dense_db)
        .arg("--screen-c")
        .arg("200")
        .arg("-o")
        .arg(format!("{}/db2", dir))
        .assert()
        .success()
        .code(0);
    let two_stage_db = format!("{}/db2.syl2db", dir);
    assert!(
        Path::new(&two_stage_db).exists(),
        "db-convert did not produce a .syl2db"
    );
    // The two-stage db should be no larger than the dense .syldb it came from
    // (dense blocks are Golomb-Rice compressed; only the small sparse index adds).
    let dense_sz = fs::metadata(&dense_db).unwrap().len();
    let two_sz = fs::metadata(&two_stage_db).unwrap().len();
    assert!(
        two_sz < dense_sz,
        "compressed two-stage db ({} B) not smaller than dense db ({} B)",
        two_sz,
        dense_sz
    );

    // Profile --two-stage directly against the .syl2db: stage 1 screens via the
    // sparse index, stage 2 decodes only the screened genomes' dense blocks.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let output = cmd
        .arg("profile")
        .arg("--two-stage")
        .arg(&two_stage_db)
        .arg(&sample)
        .output()
        .expect("Output failed");
    assert!(output.status.success());
    let from_db2 = str::from_utf8(&output.stdout)
        .expect("not UTF-8")
        .to_string();
    assert!(from_db2.contains("e.coli-o157.fasta.gz"));

    // The detected genome set must equal a plain single-stage dense profile.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let output = cmd
        .arg("profile")
        .arg(&dense_db)
        .arg(&sample)
        .output()
        .expect("Output failed");
    let single = str::from_utf8(&output.stdout)
        .expect("not UTF-8")
        .to_string();

    let detected = |tsv: &str| -> Vec<String> {
        let mut v: Vec<String> = tsv
            .lines()
            .skip(1)
            .filter_map(|l| l.split('\t').nth(1).map(|s| s.to_string()))
            .collect();
        v.sort();
        v
    };
    assert_eq!(
        detected(&from_db2),
        detected(&single),
        "two-stage .syl2db and single-stage detected different genome sets"
    );

    // `query` must refuse a .syl2db (it is profile-only).
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("query")
        .arg(&two_stage_db)
        .arg(&sample)
        .assert()
        .failure();
}

#[serial]
#[test]
fn test_two_stage_individual_records() {
    fresh();
    let dir = "./tests/results/two_stage_indiv";
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();

    // Dense (-c 50) database built with --individual-records: e.coli-o157 has two
    // records, so multiple database entries share one file name -- the case that
    // must be preserved per record by db-convert (and rejected by the densify
    // fallback).
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-c")
        .arg("50")
        .arg("-i")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .arg("-o")
        .arg(format!("{}/db_c50", dir))
        .assert()
        .success()
        .code(0);

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-c")
        .arg("50")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);

    let dense_db = format!("{}/db_c50.syldb", dir);
    let sample = format!("{}/o157_reads.fastq.gz.sylsp", dir);

    // Convert to a two-stage db (per-record blocks are written individually).
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("db-convert")
        .arg(&dense_db)
        .arg("--screen-c")
        .arg("200")
        .arg("-o")
        .arg(format!("{}/db2", dir))
        .assert()
        .success()
        .code(0);
    let two_stage_db = format!("{}/db2.syl2db", dir);

    // Per-record key = genome_file (col 2) + contig name (last col).
    let detected = |tsv: &str| -> Vec<String> {
        let mut v: Vec<String> = tsv
            .lines()
            .skip(1)
            .filter_map(|l| {
                let cols: Vec<&str> = l.split('\t').collect();
                if cols.len() < 2 {
                    return None;
                }
                Some(format!("{}\t{}", cols[1], cols[cols.len() - 1]))
            })
            .collect();
        v.sort();
        v
    };

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let two = cmd
        .arg("profile")
        .arg("--two-stage")
        .arg(&two_stage_db)
        .arg(&sample)
        .output()
        .expect("Output failed");
    assert!(two.status.success());
    let two = str::from_utf8(&two.stdout).expect("not UTF-8").to_string();
    assert!(two.contains("e.coli-o157.fasta.gz"));

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let single = cmd
        .arg("profile")
        .arg(&dense_db)
        .arg(&sample)
        .output()
        .expect("Output failed");
    let single = str::from_utf8(&single.stdout)
        .expect("not UTF-8")
        .to_string();

    // db-convert + two-stage must reproduce single-stage per-record detections
    // (no collapsing/merging of records sharing a file name).
    assert_eq!(
        detected(&two),
        detected(&single),
        "two-stage .syl2db lost or merged individual records vs single-stage"
    );

    // The densify fallback (raw .syldb --two-stage, no db-convert) cannot handle
    // individual records and must error rather than silently corrupt them.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("profile")
        .arg("--two-stage")
        .arg("--dense-c")
        .arg("50")
        .arg(&dense_db)
        .arg(&sample)
        .assert()
        .failure();
}

fn num_sketched_kmers(inspect_stdout: &str) -> u64 {
    inspect_stdout
        .lines()
        .find_map(|l| l.trim().strip_prefix("num_sketched_kmers:"))
        .and_then(|v| v.trim().parse().ok())
        .expect("inspect output had no num_sketched_kmers")
}

#[serial]
#[test]
fn test_sketch_merge_single_and_paired() {
    let dir = "./tests/results/test_merge_dir";
    if Path::new(dir).exists() {
        let _ = fs::remove_dir_all(dir);
    }
    fs::create_dir_all(dir).unwrap();
    let out = format!("{}/combined", dir);

    // --merge collapses paired + single-end inputs into ONE compressed sketch, with
    // the --compressed-database value used as the single output file path.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--merge")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("--compressed-database")
        .arg(&out)
        .arg("-S")
        .arg("merged_sample")
        .assert()
        .success()
        .code(0);

    let merged_path = format!("{}.sylspc", out);
    assert!(
        Path::new(&merged_path).exists(),
        "merged single output file was not written"
    );

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let inspect = cmd
        .arg("inspect")
        .arg(&merged_path)
        .output()
        .expect("Output failed");
    let inspect = str::from_utf8(&inspect.stdout).expect("Output was not valid UTF-8");
    assert!(inspect.contains("merged_sample"));
    assert!(inspect.contains("paired: true"));
    let merged_kmers = num_sketched_kmers(inspect);

    // Sketching the same inputs individually and then combining them with the `merge`
    // subcommand must yield the identical distinct-k-mer count as the one-pass --merge.
    // Use compressed (*.sylspc) sketches: legacy uncompressed *.sylsp lack read metadata
    // and are rejected as merge inputs.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("--compressed-database")
        .arg(dir)
        .assert()
        .success()
        .code(0);

    let via_sub = format!("{}/via_sub", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("merge")
        .arg(format!("{}/k12_R1.fq.paired.sylspc", dir))
        .arg(format!("{}/o157_reads.fastq.gz.sylspc", dir))
        .arg("-o")
        .arg(&via_sub)
        .assert()
        .success()
        .code(0);

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let inspect = cmd
        .arg("inspect")
        .arg(format!("{}.sylspc", via_sub))
        .output()
        .expect("Output failed");
    let inspect = str::from_utf8(&inspect.stdout).expect("Output was not valid UTF-8");
    assert_eq!(
        merged_kmers,
        num_sketched_kmers(inspect),
        "one-pass --merge disagreed with sketch + merge subcommand"
    );

    // --merge with a default output directory (./) is rejected: it needs an explicit path.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--merge")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .assert()
        .failure();

    // A directory-style value (trailing slash, or an existing directory) is likewise
    // rejected rather than writing a hidden file like `<dir>/.sylspc`.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--merge")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("--compressed-database")
        .arg(format!("{}/", dir))
        .assert()
        .failure();

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--merge")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("--compressed-database")
        .arg(dir)
        .assert()
        .failure();

    // --merge must not write a partial sketch when one input fails to sketch. A valid read
    // combined with a missing file must fail the whole merge rather than silently sketching
    // only the good input.
    let partial_out = format!("{}/partial", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--merge")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-r")
        .arg("test_files/does_not_exist.fq")
        .arg("--compressed-database")
        .arg(&partial_out)
        .assert()
        .failure();
    assert!(
        !Path::new(&format!("{}.sylspc", partial_out)).exists(),
        "a partial merged sketch was written despite a failed input"
    );
}

/// Sketching must croak when read inputs were requested but every stream turned out
/// empty -- otherwise a run where an upstream producer emitted nothing (e.g. through a
/// FIFO) would silently exit 0 having sketched nothing. A single non-empty stream among
/// empty ones is fine; only a zero total is an error.
#[serial]
#[test]
fn test_sketch_empty_reads_croaks() {
    let dir = "./tests/results/test_empty_reads_dir";
    if Path::new(dir).exists() {
        let _ = fs::remove_dir_all(dir);
    }
    fs::create_dir_all(dir).unwrap();
    let empty = format!("{}/empty.fq", dir);
    fs::write(&empty, b"").unwrap();

    // A single empty read stream: nothing was sketched, so exit non-zero rather than
    // silently succeed.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-r")
        .arg(&empty)
        .arg("-d")
        .arg(dir)
        .assert()
        .failure()
        .code(1);

    // --merge over only-empty read streams likewise fails (and writes no output file).
    let merged = format!("{}/merged", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--merge")
        .arg("-r")
        .arg(&empty)
        .arg("--compressed-database")
        .arg(&merged)
        .assert()
        .failure();
    assert!(
        !Path::new(&format!("{}.sylspc", merged)).exists(),
        "an empty merged sketch was written when no reads were sketched"
    );

    // A non-empty stream alongside the empty one is not an error: some reads were sketched.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-r")
        .arg(&empty)
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);
}

/// `--tolerate-empty-inputs` lets a zero-read stream (e.g. a FIFO from an SRA that has no
/// unpaired reads) count as a valid empty sketch instead of aborting a `--merge`. Without
/// the flag such an input fails the merge; with it, the merge succeeds and equals the
/// non-empty input sketched alone.
#[serial]
#[test]
fn test_sketch_merge_tolerate_empty_inputs() {
    let dir = "./tests/results/test_merge_tolerate_empty_dir";
    if Path::new(dir).exists() {
        let _ = fs::remove_dir_all(dir);
    }
    fs::create_dir_all(dir).unwrap();
    let empty = format!("{}/empty.fq", dir);
    fs::write(&empty, b"").unwrap();

    // Reference: the good input merged on its own gives the distinct-k-mer count that an
    // empty input, being zero reads, must not change.
    let good_only = format!("{}/good_only", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--merge")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("--compressed-database")
        .arg(&good_only)
        .arg("-S")
        .arg("s")
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let inspect = cmd
        .arg("inspect")
        .arg(format!("{}.sylspc", good_only))
        .output()
        .expect("Output failed");
    let good_kmers = num_sketched_kmers(str::from_utf8(&inspect.stdout).unwrap());

    // Bug repro: a good input alongside an empty one aborts the merge without the flag,
    // and writes no output.
    let no_flag = format!("{}/no_flag", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--merge")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-r")
        .arg(&empty)
        .arg("--compressed-database")
        .arg(&no_flag)
        .assert()
        .failure();
    assert!(
        !Path::new(&format!("{}.sylspc", no_flag)).exists(),
        "a merged sketch was written despite an untolerated empty input"
    );

    // With --tolerate-empty-inputs the empty stream is a valid zero-read input: the merge
    // succeeds and matches the good-input-only k-mer count.
    let tolerated = format!("{}/tolerated", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--merge")
        .arg("--tolerate-empty-inputs")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-r")
        .arg(&empty)
        .arg("--compressed-database")
        .arg(&tolerated)
        .arg("-S")
        .arg("s")
        .assert()
        .success()
        .code(0);
    assert!(Path::new(&format!("{}.sylspc", tolerated)).exists());
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let inspect = cmd
        .arg("inspect")
        .arg(format!("{}.sylspc", tolerated))
        .output()
        .expect("Output failed");
    assert_eq!(
        good_kmers,
        num_sketched_kmers(str::from_utf8(&inspect.stdout).unwrap()),
        "an empty tolerated input changed the merged k-mer count"
    );

    // The flag also covers empty paired-end mates and an empty interleaved stream mixed in
    // with a good single-end input -- all empties are zero reads, so the k-mer count is
    // still that of the good input alone.
    let empty2 = format!("{}/empty2.fq", dir);
    fs::write(&empty2, b"").unwrap();
    let multi = format!("{}/multi", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--merge")
        .arg("--tolerate-empty-inputs")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-1")
        .arg(&empty)
        .arg("-2")
        .arg(&empty2)
        .arg("--interleaved")
        .arg(&empty)
        .arg("--compressed-database")
        .arg(&multi)
        .arg("-S")
        .arg("s")
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let inspect = cmd
        .arg("inspect")
        .arg(format!("{}.sylspc", multi))
        .output()
        .expect("Output failed");
    assert_eq!(
        good_kmers,
        num_sketched_kmers(str::from_utf8(&inspect.stdout).unwrap()),
        "empty paired/interleaved inputs changed the merged k-mer count"
    );

    // A genuine mismatch -- only one mate of a pair empty -- is a real error, not an empty
    // input, so it must still fail even with the flag.
    let mismatch = format!("{}/mismatch", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--merge")
        .arg("--tolerate-empty-inputs")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("-2")
        .arg(&empty)
        .arg("--compressed-database")
        .arg(&mismatch)
        .assert()
        .failure();

    // Non-merge mode: a run where every input is empty still fails, even with the flag,
    // and must leave NO per-sample sketch files behind -- a tolerated empty input is a
    // zero-read sketch that is not written, so a failed run has no stray artifacts.
    let nm_dir = format!("{}/nonmerge_all_empty", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--tolerate-empty-inputs")
        .arg("-r")
        .arg(&empty)
        .arg("-r")
        .arg(&empty2)
        .arg("-d")
        .arg(&nm_dir)
        .assert()
        .failure();
    let stray: Vec<_> = fs::read_dir(&nm_dir)
        .map(|rd| {
            rd.filter_map(|e| e.ok())
                .map(|e| e.file_name().to_string_lossy().into_owned())
                .filter(|n| n.ends_with(".sylsp") || n.ends_with(".sylspc"))
                .collect()
        })
        .unwrap_or_default();
    assert!(
        stray.is_empty(),
        "a failed all-empty non-merge run left stray sketch files: {:?}",
        stray
    );

    // The all-empty check spans input categories: a run whose single-end AND paired-end
    // inputs are all empty must still fail even with the flag -- the flag only tolerates
    // empties when some other input has reads, never a wholly empty run. Holds for --merge
    // (no output file) and non-merge alike.
    let mixed_merge = format!("{}/mixed_merge", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--merge")
        .arg("--tolerate-empty-inputs")
        .arg("-r")
        .arg(&empty)
        .arg("-1")
        .arg(&empty)
        .arg("-2")
        .arg(&empty2)
        .arg("--compressed-database")
        .arg(&mixed_merge)
        .assert()
        .failure();
    assert!(
        !Path::new(&format!("{}.sylspc", mixed_merge)).exists(),
        "an all-empty single+paired --merge run wrote an output file"
    );

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--tolerate-empty-inputs")
        .arg("-r")
        .arg(&empty)
        .arg("-1")
        .arg(&empty)
        .arg("-2")
        .arg(&empty2)
        .arg("-d")
        .arg(format!("{}/mixed_nm", dir))
        .assert()
        .failure();
}

/// Extract the (single) data row's Sample_file (col 1) and Containment_ind (col 12)
/// from a `profile` TSV, skipping the header line.
fn profile_sample_and_containment(stdout: &str) -> (String, String) {
    let row = stdout
        .lines()
        .find(|l| !l.starts_with("Sample_file") && !l.trim().is_empty())
        .expect("no profile data row");
    let cols: Vec<&str> = row.split('\t').collect();
    (cols[0].to_string(), cols[11].to_string())
}

#[serial]
#[test]
fn test_profile_merge_single_and_paired() {
    let dir = "./tests/results/test_profile_merge_dir";
    if Path::new(dir).exists() {
        let _ = fs::remove_dir_all(dir);
    }
    fs::create_dir_all(dir).unwrap();

    // genome database
    let db = format!("{}/db", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-g")
        .arg("test_files/e.coli-K12.fasta.gz")
        .arg("-o")
        .arg(&db)
        .assert()
        .success()
        .code(0);
    let db_file = format!("{}.syldb", db);

    // profile --merge over paired + single-end reads must emit exactly ONE sample
    // named by -S, not one row-set per input.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let out = cmd
        .arg("profile")
        .arg(&db_file)
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("--merge")
        .arg("-S")
        .arg("combined")
        .output()
        .expect("Output failed");
    let out = str::from_utf8(&out.stdout).expect("Output was not valid UTF-8");
    let data_rows = out
        .lines()
        .filter(|l| !l.starts_with("Sample_file") && !l.trim().is_empty())
        .count();
    assert_eq!(data_rows, 1, "--merge must produce a single merged sample");
    let (name, containment_merge) = profile_sample_and_containment(out);
    assert_eq!(name, "combined");

    // Cross-check: sketch the same inputs, combine them with the `merge` subcommand,
    // and profile the result -- the containment index must be identical. Use compressed
    // (*.sylspc) sketches: legacy *.sylsp lack read metadata and cannot be merged.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("--compressed-database")
        .arg(dir)
        .assert()
        .success()
        .code(0);
    let merged = format!("{}/m", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("merge")
        .arg(format!("{}/k12_R1.fq.paired.sylspc", dir))
        .arg(format!("{}/o157_reads.fastq.gz.sylspc", dir))
        .arg("-o")
        .arg(&merged)
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    let out = cmd
        .arg("profile")
        .arg(&db_file)
        .arg(format!("{}.sylspc", merged))
        .output()
        .expect("Output failed");
    let out = str::from_utf8(&out.stdout).expect("Output was not valid UTF-8");
    let (_, containment_sub) = profile_sample_and_containment(out);
    assert_eq!(
        containment_merge, containment_sub,
        "profile --merge disagreed with merge subcommand + profile"
    );

    // --merge must reject inputs whose sampling rates disagree: summing sketches made
    // at different -c would silently corrupt containment. Here a c=100 pre-sketched
    // sample is mixed with c=200 raw reads against the c=200 database. Use a compressed
    // (*.sylspc) sketch so it carries read metadata and reaches the c-mismatch check.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-c")
        .arg("100")
        .arg("--compressed-database")
        .arg(format!("{}/c100", dir))
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("profile")
        .arg(&db_file)
        .arg(format!("{}/c100/o157_reads.fastq.gz.sylspc", dir))
        .arg("-r")
        .arg("test_files/k12_R1.fq")
        .arg("--merge")
        .assert()
        .failure();

    // --merge must not silently drop an input that fails to load. A pre-sketched sample
    // whose c (300) exceeds the database's (200) is unusable; mixed with a valid raw read
    // the whole merge must fail rather than profile a merged sample that excludes it.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-c")
        .arg("300")
        .arg("--compressed-database")
        .arg(format!("{}/c300", dir))
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("profile")
        .arg(&db_file)
        .arg(format!("{}/c300/o157_reads.fastq.gz.sylspc", dir))
        .arg("-r")
        .arg("test_files/k12_R1.fq")
        .arg("--merge")
        .assert()
        .failure();
}

/// Legacy uncompressed *.sylsp sketches carry no read count, so they cannot be merged
/// (read length is undeterminable). Both the `merge` subcommand and `profile --merge`
/// must reject them rather than silently corrupt the merged read length.
#[serial]
#[test]
fn test_merge_rejects_legacy_sylsp() {
    let dir = "./tests/results/test_merge_legacy_dir";
    if Path::new(dir).exists() {
        let _ = fs::remove_dir_all(dir);
    }
    fs::create_dir_all(dir).unwrap();

    // Sketch two samples to legacy uncompressed *.sylsp (the default -d encoding).
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-r")
        .arg("test_files/k12_R1.fq")
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);

    let legacy_a = format!("{}/o157_reads.fastq.gz.sylsp", dir);
    let legacy_b = format!("{}/k12_R1.fq.sylsp", dir);

    // The `merge` subcommand rejects a legacy *.sylsp input.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("merge")
        .arg(&legacy_a)
        .arg(&legacy_b)
        .arg("-o")
        .arg(format!("{}/merged", dir))
        .assert()
        .failure();

    // profile --merge rejects a legacy *.sylsp input mixed with raw reads.
    let db = format!("{}/db", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-g")
        .arg("test_files/e.coli-K12.fasta.gz")
        .arg("-o")
        .arg(&db)
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("profile")
        .arg(format!("{}.syldb", db))
        .arg(&legacy_a)
        .arg("-r")
        .arg("test_files/k12_R1.fq")
        .arg("--merge")
        .assert()
        .failure();
}

/// `sketch --merge` selects the output encoding from the output flags; an explicit output
/// path whose suffix names a different encoding is unreadable later, so it is rejected.
#[serial]
#[test]
fn test_sketch_merge_rejects_conflicting_suffix() {
    let dir = "./tests/results/test_merge_suffix_dir";
    if Path::new(dir).exists() {
        let _ = fs::remove_dir_all(dir);
    }
    fs::create_dir_all(dir).unwrap();

    // --compressed-database selects the *.sylspc encoding, but the path ends in *.sylsp.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("--merge")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("--compressed-database")
        .arg(format!("{}/out.sylsp", dir))
        .assert()
        .failure();

    // The `merge` subcommand no longer produces legacy *.sylsp, so a *.sylsp output path
    // is rejected rather than silently downgraded.
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("--compressed-database")
        .arg(dir)
        .assert()
        .success()
        .code(0);
    let sample = format!("{}/o157_reads.fastq.gz.sylspc", dir);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("merge")
        .arg(&sample)
        .arg(&sample)
        .arg("-o")
        .arg(format!("{}/merged.sylsp", dir))
        .assert()
        .failure();
}

/// Chop `path` down to `keep` bytes, as a partial write or a bad copy would.
fn truncate_file(path: &str, keep: u64) {
    let f = fs::OpenOptions::new()
        .write(true)
        .open(path)
        .unwrap_or_else(|_| panic!("could not open {} for truncation", path));
    f.set_len(keep).unwrap();
}

fn file_len(path: &str) -> u64 {
    fs::metadata(path).unwrap().len()
}

/// A truncated *.sylspc must be rejected, not silently read as a short sketch.
/// Both a chopped-in-half file and one missing only the zstd frame's trailing
/// checksum are corrupt; the latter still holds every payload byte the reader
/// wants, so it is only caught because the reader drains the frame to its end.
#[serial]
#[test]
fn test_truncated_sylspc_is_rejected() {
    let dir = "./tests/results/truncated_sylspc";
    let _ = fs::remove_dir_all(dir);

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("test_files/e.coli-K12.fasta.gz")
        .arg("test_files/e.coli-o157.fasta.gz")
        .arg("-o")
        .arg(format!("{}/db", dir))
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-r")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("--compressed-database")
        .arg(dir)
        .assert()
        .success()
        .code(0);

    let db = format!("{}/db.syldb", dir);
    let sample = format!("{}/o157_reads.fastq.gz.sylspc", dir);
    let intact_len = file_len(&sample);

    // the intact sketch profiles fine
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("profile")
        .arg(&db)
        .arg(&sample)
        .assert()
        .success()
        .code(0);

    // losing only the trailing checksum is still corruption
    let intact = fs::read(&sample).unwrap();
    truncate_file(&sample, intact_len - 4);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("profile").arg(&db).arg(&sample).assert().failure();

    // ... as is losing half the file
    fs::write(&sample, &intact).unwrap();
    truncate_file(&sample, intact_len / 2);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("profile").arg(&db).arg(&sample).assert().failure();

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("inspect").arg(&sample).assert().failure();
}

/// A truncated *.sylspr must be rejected by every reader of the format: `query`
/// against its reference, and `ref-compress --inspect`.
#[serial]
#[test]
fn test_truncated_sylspr_is_rejected() {
    let dir = "./tests/results/truncated_sylspr";
    let _ = fs::remove_dir_all(dir);

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("test_files/e.coli-K12.fasta.gz")
        .arg("test_files/e.coli-o157.fasta.gz")
        .arg("-o")
        .arg(format!("{}/db", dir))
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("ref-build")
        .arg(format!("{}/db.syldb", dir))
        .arg("-o")
        .arg(format!("{}/ref", dir))
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("ref-compress")
        .arg(format!("{}/o157_reads.fastq.gz.sylsp", dir))
        .arg("-r")
        .arg(format!("{}/ref.sylref", dir))
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);

    let db = format!("{}/db.syldb", dir);
    let refdb = format!("{}/ref.sylref", dir);
    let sample = format!("{}/o157_reads.fastq.gz.sylspr", dir);
    let intact_len = file_len(&sample);

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("query")
        .arg(&db)
        .arg(&sample)
        .arg("--reference")
        .arg(&refdb)
        .assert()
        .success()
        .code(0);

    // losing only the trailing checksum is still corruption
    let intact = fs::read(&sample).unwrap();
    truncate_file(&sample, intact_len - 4);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("query")
        .arg(&db)
        .arg(&sample)
        .arg("--reference")
        .arg(&refdb)
        .assert()
        .failure();

    // ... as is losing half the file
    fs::write(&sample, &intact).unwrap();
    truncate_file(&sample, intact_len / 2);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("query")
        .arg(&db)
        .arg(&sample)
        .arg("--reference")
        .arg(&refdb)
        .assert()
        .failure();

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("ref-compress")
        .arg("--inspect")
        .arg(&sample)
        .assert()
        .failure();
}

/// Flip a bit in the byte at `offset`, as bit rot or a bad disk would.
fn flip_byte(path: &str, offset: u64) {
    let mut bytes = fs::read(path).unwrap();
    bytes[offset as usize] ^= 0x01;
    fs::write(path, &bytes).unwrap();
}

/// A corrupt *.sylref must be caught by `inspect`. The seekable databases are read
/// a block at a time, so nothing validates them end to end in normal use; the
/// whole-file checksum in the header is what `inspect` checks.
#[serial]
#[test]
fn test_corrupt_sylref_is_rejected() {
    let dir = "./tests/results/corrupt_sylref";
    let _ = fs::remove_dir_all(dir);

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("test_files/e.coli-K12.fasta.gz")
        .arg("test_files/e.coli-o157.fasta.gz")
        .arg("-o")
        .arg(format!("{}/db", dir))
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("ref-build")
        .arg(format!("{}/db.syldb", dir))
        .arg("-o")
        .arg(format!("{}/ref", dir))
        .assert()
        .success()
        .code(0);

    let refdb = format!("{}/ref.sylref", dir);
    let intact = fs::read(&refdb).unwrap();
    let intact_len = intact.len() as u64;

    // the intact reference inspects cleanly, and says so
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("inspect")
        .arg(&refdb)
        .assert()
        .success()
        .code(0)
        .stdout(predicate::str::contains("checksum: ok"));

    // a single flipped bit in the body is caught, even though every offset in the
    // header still points somewhere valid
    flip_byte(&refdb, intact_len / 2);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("inspect")
        .arg(&refdb)
        .assert()
        .failure()
        .stderr(predicate::str::contains("corrupt"));

    // ... as is a truncated file
    fs::write(&refdb, &intact).unwrap();
    truncate_file(&refdb, intact_len / 2);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("inspect").arg(&refdb).assert().failure();
}

/// As above for the two-stage genome database: a corrupt *.syl2db must be caught by
/// `inspect`.
#[serial]
#[test]
fn test_corrupt_syl2db_is_rejected() {
    let dir = "./tests/results/corrupt_syl2db";
    let _ = fs::remove_dir_all(dir);

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("sketch")
        .arg("-c")
        .arg("50")
        .arg("test_files/e.coli-K12.fasta.gz")
        .arg("test_files/e.coli-o157.fasta.gz")
        .arg("-o")
        .arg(format!("{}/db", dir))
        .arg("-d")
        .arg(dir)
        .assert()
        .success()
        .code(0);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("db-convert")
        .arg(format!("{}/db.syldb", dir))
        .arg("--screen-c")
        .arg("200")
        .arg("-o")
        .arg(format!("{}/db2", dir))
        .assert()
        .success()
        .code(0);

    let two = format!("{}/db2.syl2db", dir);
    let intact = fs::read(&two).unwrap();
    let intact_len = intact.len() as u64;

    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("inspect")
        .arg(&two)
        .assert()
        .success()
        .code(0)
        .stdout(predicate::str::contains("checksum: ok"));

    flip_byte(&two, intact_len / 2);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("inspect")
        .arg(&two)
        .assert()
        .failure()
        .stderr(predicate::str::contains("corrupt"));

    fs::write(&two, &intact).unwrap();
    truncate_file(&two, intact_len / 2);
    let mut cmd = Command::cargo_bin("weebill").unwrap();
    cmd.arg("inspect").arg(&two).assert().failure();
}
