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
fn test_merge_basic() {
    fresh();

    // Sketch paired-end reads
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("-1").arg("test_files/k12_R1.fq")
        .arg("-2").arg("test_files/k12_R2.fq")
        .arg("-d").arg("./tests/results/test_sketch_dir")
        .assert()
        .success();

    // Sketch single-end reads
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("-r").arg("test_files/o157_reads.fastq.gz")
        .arg("-d").arg("./tests/results/test_sketch_dir")
        .assert()
        .success();

    let paired_sketch = "./tests/results/test_sketch_dir/k12_R1.fq.paired.sylsp";
    let single_sketch = "./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp";
    let merged_output = "./tests/results/test_sketch_dir/merged.sylsp";

    assert!(Path::new(paired_sketch).exists());
    assert!(Path::new(single_sketch).exists());

    // Merge the two sketches
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("merge")
        .arg(paired_sketch)
        .arg(single_sketch)
        .arg("-o").arg(merged_output)
        .assert()
        .success();

    assert!(Path::new(merged_output).exists());

    // Verify merged sketch works with profile
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd
        .arg("profile")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg(merged_output)
        .output()
        .expect("Output failed");
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    // Should produce output (header + at least one result line)
    assert!(stdout.matches('\n').count() >= 2);

    fresh();
}

#[serial]
#[test]
fn test_merge_mismatched_params() {
    fresh();

    // Sketch with c=200 (default)
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("-r").arg("test_files/o157_reads.fastq.gz")
        .arg("-d").arg("./tests/results/test_sketch_dir")
        .assert()
        .success();

    // Sketch with c=100
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("-1").arg("test_files/t1.fq")
        .arg("-2").arg("test_files/t2.fq")
        .arg("-c").arg("100")
        .arg("-d").arg("./tests/results/test_sketch_dir/c100")
        .assert()
        .success();

    // Merge should fail due to different c values
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("merge")
        .arg("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp")
        .arg("./tests/results/test_sketch_dir/c100/t1.fq.paired.sylsp")
        .arg("-o").arg("./tests/results/test_sketch_dir/bad_merge.sylsp")
        .assert()
        .failure()
        .code(1);

    fresh();
}

#[serial]
#[test]
fn test_merge_sample_name() {
    fresh();

    // Sketch two single-end files
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("-r").arg("test_files/o157_reads.fastq.gz")
        .arg("-d").arg("./tests/results/test_sketch_dir")
        .assert()
        .success();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("-1").arg("test_files/t1.fq")
        .arg("-2").arg("test_files/t2.fq")
        .arg("-d").arg("./tests/results/test_sketch_dir")
        .assert()
        .success();

    let merged_output = "./tests/results/test_sketch_dir/named_merge.sylsp";

    // Merge with custom sample name
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("merge")
        .arg("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp")
        .arg("./tests/results/test_sketch_dir/t1.fq.paired.sylsp")
        .arg("-o").arg(merged_output)
        .arg("-S").arg("MY_MERGED_SAMPLE")
        .assert()
        .success();

    // Profile with merged sketch and check sample name appears
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd
        .arg("profile")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg(merged_output)
        .output()
        .expect("Output failed");
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    assert!(stdout.contains("MY_MERGED_SAMPLE"));

    fresh();
}

#[serial]
#[test]
fn test_merge_inspect_num_reads() {
    fresh();

    // Sketch paired-end reads
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("-1").arg("test_files/k12_R1.fq")
        .arg("-2").arg("test_files/k12_R2.fq")
        .arg("-d").arg("./tests/results/test_sketch_dir")
        .assert()
        .success();

    // Inspect should show num_reads
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd
        .arg("inspect")
        .arg("./tests/results/test_sketch_dir/k12_R1.fq.paired.sylsp")
        .output()
        .expect("Output failed");
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    assert!(stdout.contains("num_reads"));

    fresh();
}

#[serial]
#[test]
fn test_merge_produces_valid_sketch() {
    fresh();

    // Sketch single-end reads into two separate sketches
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("-1").arg("test_files/k12_R1.fq")
        .arg("-2").arg("test_files/k12_R2.fq")
        .arg("-d").arg("./tests/results/test_sketch_dir")
        .assert()
        .success();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("-1").arg("test_files/t1.fq")
        .arg("-2").arg("test_files/t2.fq")
        .arg("-d").arg("./tests/results/test_sketch_dir")
        .assert()
        .success();

    let merged_output = "./tests/results/test_sketch_dir/merged_valid.sylsp";

    // Merge
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("merge")
        .arg("./tests/results/test_sketch_dir/k12_R1.fq.paired.sylsp")
        .arg("./tests/results/test_sketch_dir/t1.fq.paired.sylsp")
        .arg("-o").arg(merged_output)
        .assert()
        .success();

    // Inspect merged sketch to verify it's valid
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd
        .arg("inspect")
        .arg(merged_output)
        .output()
        .expect("Output failed");
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    assert!(stdout.contains("num_reads"));
    assert!(stdout.contains("paired: true"));

    // Query with merged sketch should work
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("query")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg(merged_output)
        .assert()
        .success();

    fresh();
}

#[serial]
#[test]
fn test_merge_too_few_files() {
    fresh();

    // Sketch one file
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("-r").arg("test_files/o157_reads.fastq.gz")
        .arg("-d").arg("./tests/results/test_sketch_dir")
        .assert()
        .success();

    // Merge with only one file should fail
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("merge")
        .arg("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp")
        .arg("-o").arg("./tests/results/test_sketch_dir/single_merge.sylsp")
        .assert()
        .failure()
        .code(1);

    fresh();
}

#[serial]
#[test]
fn test_interleaved_sketch_matches_paired() {
    fresh();

    // Sketch paired-end reads
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("-1").arg("test_files/k12_R1.fq")
        .arg("-2").arg("test_files/k12_R2.fq")
        .arg("-d").arg("./tests/results/test_sketch_dir")
        .assert()
        .success();

    // Sketch interleaved reads
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("--interleaved").arg("test_files/k12_interleaved.fq")
        .arg("-d").arg("./tests/results/test_sketch_dir/interleaved")
        .assert()
        .success();

    assert!(Path::new("./tests/results/test_sketch_dir/interleaved/k12_interleaved.fq.paired.sylsp").exists(), "Interleaved output not created");

    // Profile both against a genome and compare output
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output_paired = cmd
        .arg("profile")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./tests/results/test_sketch_dir/k12_R1.fq.paired.sylsp")
        .output()
        .expect("Output failed");
    let stdout_paired = str::from_utf8(&output_paired.stdout).expect("Output was not valid UTF-8");

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output_interleaved = cmd
        .arg("profile")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./tests/results/test_sketch_dir/interleaved/k12_interleaved.fq.paired.sylsp")
        .output()
        .expect("Output failed");
    let stdout_interleaved = str::from_utf8(&output_interleaved.stdout).expect("Output was not valid UTF-8");

    // Both should produce output with same number of lines
    assert_eq!(stdout_paired.matches('\n').count(), stdout_interleaved.matches('\n').count());
    // Both should have at least a header + result
    assert!(stdout_paired.matches('\n').count() >= 2);

    fresh();
}

#[serial]
#[test]
fn test_interleaved_sketch_t_files() {
    fresh();

    // Sketch interleaved t files (these have same read names without /1 /2 suffix)
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("--interleaved").arg("test_files/t_interleaved.fq")
        .arg("-d").arg("./tests/results/test_sketch_dir")
        .assert()
        .success();

    assert!(Path::new("./tests/results/test_sketch_dir/t_interleaved.fq.paired.sylsp").exists());

    // Verify it works with profile
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("profile")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./tests/results/test_sketch_dir/t_interleaved.fq.paired.sylsp")
        .assert()
        .success();

    fresh();
}

#[serial]
#[test]
fn test_interleaved_bad_names_fails() {
    fresh();

    // bad_interleaved.fq has two R1 reads with different names — should fail
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch")
        .arg("--interleaved").arg("test_files/bad_interleaved.fq")
        .arg("-d").arg("./tests/results/test_sketch_dir")
        .assert()
        .failure()
        .code(1);

    fresh();
}

#[serial]
#[test]
fn test_interleaved_profile_direct() {
    fresh();

    // Use --interleaved directly with profile (no pre-sketching)
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd
        .arg("profile")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("--interleaved").arg("test_files/k12_interleaved.fq")
        .output()
        .expect("Output failed");
    assert!(output.status.success());
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    assert!(stdout.matches('\n').count() >= 2);

    fresh();
}

#[serial]
#[test]
fn test_interleaved_query_direct() {
    fresh();

    // Use --interleaved directly with query
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd
        .arg("query")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("--interleaved").arg("test_files/k12_interleaved.fq")
        .output()
        .expect("Output failed");
    assert!(output.status.success());
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    assert!(stdout.matches('\n').count() >= 2);

    fresh();
}
