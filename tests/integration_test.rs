use assert_cmd::prelude::*; // Add methods on commands
use std::str;
use std::fs;
use std::path::Path;
use serial_test::serial;
use std::process::Command; // Run programs

fn fresh(){
    Command::new("rm")
        .arg("-r")
        .args(["./tests/results/test_sketch_dir"])
        .spawn();
}

#[serial]
#[test]
fn test_sketch_commands() {
   let mut cmd = Command::cargo_bin("sylph").unwrap();
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

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("-l")
        .arg("./test_files/list.txt")
        .assert();
    assert.success().code(0);


    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("-i")
        .arg("-m")
        .arg("90")
        .assert();
    assert.success().code(0);

    let mut cmd= Command::cargo_bin("sylph").unwrap();
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
    assert!(Path::new("./tests/results/test_sketch_dir/t1.fq.paired.sylsp").exists(), "Output file was not created");
    fresh();

    let mut cmd= Command::cargo_bin("sylph").unwrap();
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
    assert!(Path::new("./tests/results/test_sketch_dir/t1.fq.paired.sylsp").exists(), "Output file was not created");

    fresh();
    let mut cmd= Command::cargo_bin("sylph").unwrap();
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
    assert!(Path::new("./tests/results/test_sketch_dir/t2.fq.sylsp").exists(), "Output file was not created");
    assert!(Path::new("./tests/results/test_sketch_dir/testdb.syldb").exists(), "Output file was not created");
}

#[serial]
#[test]
fn test_profile_vs_query(){
    fresh();

    let mut output = Command::cargo_bin("sylph").unwrap();
    let output = output
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .output()
        .expect("Output failed");
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    dbg!(stdout.matches('\n').count());
    assert!(stdout.matches('\n').count() == 2);

    let mut output = Command::cargo_bin("sylph").unwrap();
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
    println!("{}",stdout);
    assert!(stdout.matches('\n').count() == 4);
}

#[serial]
#[test]
fn test_sketch_list(){
    fresh();
    let mut cmd = Command::cargo_bin("sylph").unwrap();
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
    assert!(Path::new("./tests/results/test_sketch_dir/e.coli-EC590.fasta.gz.sylsp").exists(), "Output file was not created");
    assert!(Path::new("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp").exists(), "Output file was not created");
    assert!(!Path::new("./tests/results/test_sketch_dir/db.syldb").exists(), "Output file was created");
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
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
    assert!(!Path::new("./tests/results/test_sketch_dir/e.coli-EC590.fasta.gz.sylsp").exists(), "Output file was created");
    assert!(!Path::new("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp").exists(), "Output file was created");
    assert!(Path::new("./tests/results/test_sketch_dir/db.syldb").exists(), "Output file was not created");
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("--gl")
        .arg("test_files/list.txt")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/db")
        .assert();
    assert.success().code(0);
    assert!(Path::new("./tests/results/test_sketch_dir/db.syldb").exists(), "Output file was not created");
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
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
    assert!(!Path::new("./tests/results/test_sketch_dir/db.syldb").exists(), "Output file was not created");
    assert!(Path::new("./tests/results/test_sketch_dir/e.coli-EC590.fasta.gz.sylsp").exists(), "Output file was not created");
    assert!(Path::new("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp").exists(), "Output file was not created");
    fresh();

}
#[serial]
#[test]
fn test_profile_disabling(){
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
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

    let mut output = Command::cargo_bin("sylph").unwrap();
    let assert = output
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .assert();
    assert.failure().code(1);

    let mut output = Command::cargo_bin("sylph").unwrap();
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
fn test_sketch_fasta_fastq_concord(){
    fresh();
    let mut cmd = Command::cargo_bin("sylph").unwrap();
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

    let mut output = Command::cargo_bin("sylph").unwrap();
    let out1 = output
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .output()
        .expect("Fail");

    let mut output = Command::cargo_bin("sylph").unwrap();
    let out2 = output
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .output()
        .expect("Fail");

    let mut output = Command::cargo_bin("sylph").unwrap();
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
fn test_sample_names(){
    fresh();
    let mut cmd = Command::cargo_bin("sylph").unwrap();
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
    assert!(Path::new("./tests/results/test_sketch_dir/SAMPLE_TEST.paired.sylsp").exists(), "Output file was not created");
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
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
    assert!(Path::new("./tests/results/test_sketch_dir/S1.sylsp").exists(), "Output file was not created");
    assert!(Path::new("./tests/results/test_sketch_dir/S2.sylsp").exists(), "Output file was not created");

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd
        .arg("profile")
        .arg("./tests/results/test_sketch_dir/S2.sylsp")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .output().unwrap();
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    dbg!(&stdout);
    assert!(stdout.contains("S2"));
    assert!(!stdout.contains("o157_reads"));

    let mut cmd = Command::cargo_bin("sylph").unwrap();
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
    assert!(Path::new("./tests/results/test_sketch_dir/SAMPLE_TEST_S.paired.sylsp").exists(), "Output file was not created, -S");

    let mut cmd = Command::cargo_bin("sylph").unwrap();
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
    assert!(Path::new("./tests/results/test_sketch_dir/SAMPLE_TEST_S1.paired.sylsp").exists(), "Output file was not created, -S");

    fresh();
}
#[serial]
#[test]
fn test_fpr(){
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .arg("-d ")
        .arg("./tests/results/test_sketch_dir")
        .arg("0")
        .assert();
    assert.success().code(0);
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
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
    let mut cmd = Command::cargo_bin("sylph").unwrap();
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
fn test_raw_inputs_profile_simple(){
    let mut cmd = Command::cargo_bin("sylph").unwrap();
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

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("-1")
        .arg("test_files/t1.fq")
        .assert();
    assert.failure().code(1);
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
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
fn test_estimate_read_counts(){
    let mut cmd = Command::cargo_bin("sylph").unwrap();
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
    let split : Vec<&str> = output.split('\t').collect();
    assert!(split[3].parse::<f64>().unwrap() > 1000.0);

    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
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
    let split : Vec<&str> = output.split('\t').collect();
    assert!(split[3].parse::<f64>().unwrap() < 101.00);

    fresh();

}

#[serial]
#[test]
fn test_raw_inputs_profile_with_sketch(){
    
    let mut output = Command::cargo_bin("sylph").unwrap();
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

    let mut cmd = Command::cargo_bin("sylph").unwrap();
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

    let mut output = Command::cargo_bin("sylph").unwrap();
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
fn test_inspect(){
   let mut cmd = Command::cargo_bin("sylph").unwrap();
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
    let mut cmd = Command::cargo_bin("sylph").unwrap();
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

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd
        .arg("inspect")
        .arg("./tests/results/test_sketch_dir/k12_R1.fq.paired.sylsp")
        .output()
        .expect("Output failed");

    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    assert!(stdout.contains("k12_R1.fq"));

    let mut cmd = Command::cargo_bin("sylph").unwrap();
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
fn test_refdelta_query_with_reference(){
    fresh();
    let dir = "./tests/results/test_sketch_dir";

    // sketch a database and a sample
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("test_files/e.coli-K12.fasta.gz")
        .arg("test_files/e.coli-o157.fasta.gz")
        .arg("test_files/e.coli-EC590.fasta.gz")
        .arg("-o").arg(format!("{}/db", dir))
        .arg("-d").arg(dir)
        .assert().success().code(0);
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-d").arg(dir)
        .assert().success().code(0);

    // build a reference and compress the sample sketch against it
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("ref-build")
        .arg(format!("{}/db.syldb", dir))
        .arg("-o").arg(format!("{}/ref", dir))
        .assert().success().code(0);
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("ref-compress")
        .arg(format!("{}/o157_reads.fastq.gz.sylsp", dir))
        .arg("-r").arg(format!("{}/ref.sylref", dir))
        .arg("-d").arg(dir)
        .assert().success().code(0);
    assert!(Path::new(&format!("{}/o157_reads.fastq.gz.sylspr", dir)).exists(),
        "ref-compress did not produce a .sylspr");

    // querying the .sylspr via --reference must match querying the original .sylsp
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let orig = cmd.arg("query")
        .arg(format!("{}/db.syldb", dir))
        .arg(format!("{}/o157_reads.fastq.gz.sylsp", dir))
        .output().expect("Output failed");
    let orig = str::from_utf8(&orig.stdout).expect("Output was not valid UTF-8");

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let from_ref = cmd.arg("query")
        .arg(format!("{}/db.syldb", dir))
        .arg(format!("{}/o157_reads.fastq.gz.sylspr", dir))
        .arg("--reference").arg(format!("{}/ref.sylref", dir))
        .output().expect("Output failed");
    let from_ref = str::from_utf8(&from_ref.stdout).expect("Output was not valid UTF-8");

    assert!(orig.contains("e.coli-o157.fasta.gz"));
    assert_eq!(orig, from_ref, "query of .sylspr via --reference differs from query of original .sylsp");

    // a .sylspr without --reference must fail
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("query")
        .arg(format!("{}/db.syldb", dir))
        .arg(format!("{}/o157_reads.fastq.gz.sylspr", dir))
        .assert().failure();
}

#[serial]
#[test]
fn test_two_stage_profile(){
    fresh();
    let dir = "./tests/results/two_stage";
    let _ = fs::remove_dir_all(dir);

    // Sparse (-c 200) database that retains the source fasta paths.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("200")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .arg("-o").arg(format!("{}/db_c200", dir))
        .assert().success().code(0);

    // Dense (-c 50) read sample.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-d").arg(dir)
        .assert().success().code(0);

    let db = format!("{}/db_c200.syldb", dir);
    let sample = format!("{}/o157_reads.fastq.gz.sylsp", dir);
    let cache = format!("{}/cache", dir);

    // Two-stage profile: screen at c=200, densely profile the survivors at c=50
    // by re-sketching their source fastas, caching the dense sketches.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd.arg("profile").arg("--two-stage")
        .arg("--dense-c").arg("50")
        .arg("--dense-cache").arg(&cache)
        .arg(&db).arg(&sample)
        .output().expect("Output failed");
    assert!(output.status.success());
    let two_stage = str::from_utf8(&output.stdout).expect("not UTF-8").to_string();
    // The reads are E. coli O157 -> that genome must be profiled.
    assert!(two_stage.contains("e.coli-o157.fasta.gz"));
    // The dense stage caches a per-genome sketch for each screened survivor.
    assert!(Path::new(&cache).exists());
    assert!(fs::read_dir(&cache).unwrap().count() >= 1, "dense cache was not populated");

    // Genomes detected by two-stage must equal those of a plain single-stage
    // profile of the same dense reads against a dense (-c 50) database.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .arg("-o").arg(format!("{}/db_c50", dir))
        .assert().success().code(0);
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd.arg("profile")
        .arg(format!("{}/db_c50.syldb", dir)).arg(&sample)
        .output().expect("Output failed");
    let single = str::from_utf8(&output.stdout).expect("not UTF-8").to_string();

    let detected = |tsv: &str| -> Vec<String> {
        let mut v: Vec<String> = tsv.lines().skip(1)
            .filter_map(|l| l.split('\t').nth(1).map(|s| s.to_string()))
            .collect();
        v.sort();
        v
    };
    assert_eq!(detected(&two_stage), detected(&single),
        "two-stage and single-stage detected different genome sets");
}

#[serial]
#[test]
fn test_two_stage_db_convert_and_profile(){
    fresh();
    let dir = "./tests/results/two_stage_db";
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();

    // Dense (-c 50) database carrying profiling k-mers.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .arg("-o").arg(format!("{}/db_c50", dir))
        .assert().success().code(0);

    // Dense (-c 50) read sample.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-d").arg(dir)
        .assert().success().code(0);

    let dense_db = format!("{}/db_c50.syldb", dir);
    let sample = format!("{}/o157_reads.fastq.gz.sylsp", dir);

    // Convert the dense db into a two-stage seekable database: dense blocks at
    // c=50, sparse stage-1 screen index at c=200.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("db-convert")
        .arg(&dense_db)
        .arg("--screen-c").arg("200")
        .arg("-o").arg(format!("{}/db2", dir))
        .assert().success().code(0);
    let two_stage_db = format!("{}/db2.syl2db", dir);
    assert!(Path::new(&two_stage_db).exists(), "db-convert did not produce a .syl2db");
    // The two-stage db should be no larger than the dense .syldb it came from
    // (dense blocks are Golomb-Rice compressed; only the small sparse index adds).
    let dense_sz = fs::metadata(&dense_db).unwrap().len();
    let two_sz = fs::metadata(&two_stage_db).unwrap().len();
    assert!(two_sz < dense_sz, "compressed two-stage db ({} B) not smaller than dense db ({} B)", two_sz, dense_sz);

    // Profile --two-stage directly against the .syl2db: stage 1 screens via the
    // sparse index, stage 2 decodes only the screened genomes' dense blocks.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd.arg("profile").arg("--two-stage")
        .arg(&two_stage_db).arg(&sample)
        .output().expect("Output failed");
    assert!(output.status.success());
    let from_db2 = str::from_utf8(&output.stdout).expect("not UTF-8").to_string();
    assert!(from_db2.contains("e.coli-o157.fasta.gz"));

    // The detected genome set must equal a plain single-stage dense profile.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd.arg("profile")
        .arg(&dense_db).arg(&sample)
        .output().expect("Output failed");
    let single = str::from_utf8(&output.stdout).expect("not UTF-8").to_string();

    let detected = |tsv: &str| -> Vec<String> {
        let mut v: Vec<String> = tsv.lines().skip(1)
            .filter_map(|l| l.split('\t').nth(1).map(|s| s.to_string()))
            .collect();
        v.sort();
        v
    };
    assert_eq!(detected(&from_db2), detected(&single),
        "two-stage .syl2db and single-stage detected different genome sets");

    // `query` must refuse a .syl2db (it is profile-only).
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("query").arg(&two_stage_db).arg(&sample).assert().failure();
}

#[serial]
#[test]
fn test_two_stage_individual_records(){
    fresh();
    let dir = "./tests/results/two_stage_indiv";
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();

    // Dense (-c 50) database built with --individual-records: e.coli-o157 has two
    // records, so multiple database entries share one file name -- the case that
    // must be preserved per record by db-convert (and rejected by the densify
    // fallback).
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50").arg("-i")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .arg("-o").arg(format!("{}/db_c50", dir))
        .assert().success().code(0);

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-d").arg(dir)
        .assert().success().code(0);

    let dense_db = format!("{}/db_c50.syldb", dir);
    let sample = format!("{}/o157_reads.fastq.gz.sylsp", dir);

    // Convert to a two-stage db (per-record blocks are written individually).
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("db-convert").arg(&dense_db)
        .arg("--screen-c").arg("200")
        .arg("-o").arg(format!("{}/db2", dir))
        .assert().success().code(0);
    let two_stage_db = format!("{}/db2.syl2db", dir);

    // Per-record key = genome_file (col 2) + contig name (last col).
    let detected = |tsv: &str| -> Vec<String> {
        let mut v: Vec<String> = tsv.lines().skip(1)
            .filter_map(|l| {
                let cols: Vec<&str> = l.split('\t').collect();
                if cols.len() < 2 { return None; }
                Some(format!("{}\t{}", cols[1], cols[cols.len() - 1]))
            })
            .collect();
        v.sort();
        v
    };

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let two = cmd.arg("profile").arg("--two-stage").arg(&two_stage_db).arg(&sample)
        .output().expect("Output failed");
    assert!(two.status.success());
    let two = str::from_utf8(&two.stdout).expect("not UTF-8").to_string();
    assert!(two.contains("e.coli-o157.fasta.gz"));

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let single = cmd.arg("profile").arg(&dense_db).arg(&sample)
        .output().expect("Output failed");
    let single = str::from_utf8(&single.stdout).expect("not UTF-8").to_string();

    // db-convert + two-stage must reproduce single-stage per-record detections
    // (no collapsing/merging of records sharing a file name).
    assert_eq!(detected(&two), detected(&single),
        "two-stage .syl2db lost or merged individual records vs single-stage");

    // The densify fallback (raw .syldb --two-stage, no db-convert) cannot handle
    // individual records and must error rather than silently corrupt them.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("profile").arg("--two-stage").arg("--dense-c").arg("50")
        .arg(&dense_db).arg(&sample)
        .assert().failure();
}
