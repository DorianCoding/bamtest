extern crate rand;
extern crate bam;

use std::fs::File;
use std::process::Command;
use std::time::Instant;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rand::Rng;
use glob::glob;

use bam::{RecordReader, RecordWriter};

fn compare_sam_files<P: AsRef<Path>, T: AsRef<Path>>(filename1: P, filename2: T) {
    let mut file1 = BufReader::new(File::open(filename1).unwrap());
    let mut file2 = BufReader::new(File::open(filename2).unwrap());
    let mut line1 = String::new();
    let mut line2 = String::new();

    for i in 1_usize.. {
        line1.clear();
        line2.clear();
        match (file1.read_line(&mut line1), file2.read_line(&mut line2)) {
            (Ok(x), Ok(y)) => {
                if x == 0 && y != 0 {
                    println!("Samtools output: {}", line2.trim());
                    panic!("Samtools output is longer");
                } else if x != 0 && y == 0 {
                    println!("Crate output:    {}", line1.trim());
                    panic!("Crate output is longer");
                } else if x == 0 && y == 0 {
                    break;
                }
            },
            (Ok(_), Err(e)) => panic!("Could not read samtools output: {:?}", e),
            (Err(e), Ok(_)) => panic!("Could not read crate output: {:?}", e),
            (Err(e1), Err(e2)) => panic!("Could not read both outputs: {:?}, {:?}", e1, e2),
        }
        if line1 != line2 {
            println!("Crate output:    {}", line1.trim());
            println!("Samtools output: {}", line2.trim());
            panic!("Outputs do not match on line {}", i);
        }
    }
}

fn test_indexed_reader(path: &str, additional_threads: u16) {
    let mut reader = bam::IndexedReader::build()
        .additional_threads(additional_threads)
        .from_path(path).unwrap();
    let header = reader.header().clone();

    const ITERATIONS: usize = 100;
    let mut rng = rand::thread_rng();
    let mut record = bam::Record::new();

    let output1 = format!("tests/data/tmp/bamcrate.ind_reader_t{}.sam", additional_threads);
    let output2 = format!("tests/data/tmp/samtools.ind_reader_t{}.sam", additional_threads);
    for i in 0..ITERATIONS {
        let ref_id = rng.gen_range(0, header.n_references());
        let length = header.reference_len(ref_id as u32).unwrap();
        let start = rng.gen_range(0, length);
        let end = rng.gen_range(start + 1, length + 1);

        let mut count = 0;
        let ref_name = header.reference_name(ref_id as u32).unwrap();
        println!("    Iteration {}", i);
        println!("    Fetching {}:{}-{}", ref_name, start + 1, end);

        let timer = Instant::now();
        let mut sam_writer = bam::SamWriter::from_path(&output1, header.clone())
            .unwrap();
        let mut viewer = reader.fetch(&bam::Region::new(ref_id as u32, start, end)).unwrap();
        loop {
            match viewer.read_into(&mut record) {
                Ok(true) => {},
                Ok(false) => break,
                Err(e) => panic!("{}", e),
            }
            count += 1;
            sam_writer.write(&record).unwrap();
        }
        sam_writer.finish().unwrap();
        println!("        bam::IndexedReader: {:?}", timer.elapsed());

        let timer = Instant::now();
        let mut child = Command::new("samtools")
            .args(&["view", "-h"])
            .arg(path)
            .arg(format!("{}:{}-{}", ref_name, start + 1, end))
            .args(&["-o", &output2])
            .spawn()
            .expect("Failed to run samtools view");
        let ecode = child.wait().expect("Failed to wait on samtools view");
        assert!(ecode.success());
        println!("        samtools view:      {:?}", timer.elapsed());
        println!("        total {} records", count);
        compare_sam_files(&output1, &output2);
    }
}

fn test_bam_reader(path: &str, additional_threads: u16) {
    let mut reader = bam::BamReader::from_path(path, additional_threads).unwrap();

    let mut record = bam::Record::new();
    let mut count = 0;

    let output1 = format!("tests/data/tmp/bamcrate.bam_reader_t{}.sam", additional_threads);
    let output2 = format!("tests/data/tmp/samtools.bam_reader_t{}.sam", additional_threads);
    let timer = Instant::now();
    let mut sam_writer = bam::SamWriter::from_path(&output1, reader.header().clone()).unwrap();
    loop {
        match reader.read_into(&mut record) {
            Ok(true) => {},
            Ok(false) => break,
            Err(e) => panic!("{}", e),
        }
        count += 1;
        sam_writer.write(&record).unwrap();
    }
    sam_writer.finish().unwrap();
    println!("        bam::BamReader: {:?}", timer.elapsed());

    let timer = Instant::now();
    let mut child = Command::new("samtools")
        .args(&["view", "-h"])
        .arg(path)
        .args(&["-o", &output2])
        .spawn()
        .expect("Failed to run samtools view");
    let ecode = child.wait().expect("Failed to wait on samtools view");
    assert!(ecode.success());
    println!("        samtools view:  {:?}", timer.elapsed());
    println!("        total {} records", count);
    compare_sam_files(&output1, &output2);
}

fn test_bam_to_bam(path: &str, additional_threads: u16) {
    let mut reader = bam::BamReader::from_path(path, additional_threads).unwrap();
    let mut record = bam::Record::new();
    let mut count = 0;

    let bam_output = format!("tests/data/tmp/bamcrate.bam_to_bam_t{}.bam", additional_threads);
    let output1 = format!("tests/data/tmp/bamcrate.bam_to_bam_t{}.sam", additional_threads);
    let output2 = format!("tests/data/tmp/samtools.bam_to_bam_t{}.sam", additional_threads);
    let timer = Instant::now();
    let mut bam_writer = bam::BamWriter::build()
        .additional_threads(additional_threads)
        .from_path(&bam_output, reader.header().clone())
        .unwrap();
    loop {
        match reader.read_into(&mut record) {
            Ok(true) => {},
            Ok(false) => break,
            Err(e) => panic!("{}", e),
        }
        count += 1;
        bam_writer.write(&record).unwrap();
    }
    bam_writer.finish().unwrap();
    println!("        bam::BamWriter: {:?}", timer.elapsed());

    let timer = Instant::now();
    let mut child = Command::new("samtools")
        .args(&["view", "-h"])
        .arg(&path)
        .args(&["-o", &output2])
        .spawn()
        .expect("Failed to run samtools view");
    let ecode = child.wait().expect("Failed to wait on samtools view");
    assert!(ecode.success());
    println!("        samtools view:  {:?}", timer.elapsed());
    println!("        total {} records", count);

    let mut child = Command::new("samtools")
        .args(&["view", "-h"])
        .arg(&bam_output)
        .args(&["-o", &output1])
        .spawn()
        .expect("Failed to run samtools view");
    let ecode = child.wait().expect("Failed to wait on samtools view");
    assert!(ecode.success());

    compare_sam_files(&output1, &output2);
}

fn test_sam_to_bam(path: &str) {
    let mut reader = bam::SamReader::from_path(path).unwrap();
    let mut record = bam::Record::new();
    let mut count = 0;

    let bam_output = format!("tests/data/tmp/bamcrate.sam_to_bam.bam");
    let output1 = format!("tests/data/tmp/bamcrate.sam_to_bam.sam");
    let timer = Instant::now();
    let mut bam_writer = bam::BamWriter::from_path(&bam_output, reader.header().clone())
        .unwrap();
    loop {
        match reader.read_into(&mut record) {
            Ok(true) => {},
            Ok(false) => break,
            Err(e) => panic!("{}", e),
        }
        count += 1;
        bam_writer.write(&record).unwrap();
    }
    bam_writer.finish().unwrap();
    println!("        bam::BamWriter: {:?}", timer.elapsed());

    let timer = Instant::now();
    let mut child = Command::new("samtools")
        .args(&["view", "-h"])
        .arg(&bam_output)
        .args(&["-o", &output1])
        .spawn()
        .expect("Failed to run samtools view");
    let ecode = child.wait().expect("Failed to wait on samtools view");
    assert!(ecode.success());
    println!("        samtools view:  {:?}", timer.elapsed());
    println!("        total {} records", count);

    compare_sam_files(&output1, &path);
}

#[test]
fn indexed_reader_singlethread() {
    for entry in glob("tests/data/iread*.bam").unwrap() {
        let entry = entry.unwrap();
        println!("Analyzing {}", entry.display());
        let entry_str = entry.as_os_str().to_str().unwrap();
        test_indexed_reader(entry_str, 0);
    }
}

#[test]
fn bam_reader_singlethread() {
    for entry in glob("tests/data/read*.bam").unwrap() {
        let entry = entry.unwrap();
        println!("Analyzing {}", entry.display());
        let entry_str = entry.as_os_str().to_str().unwrap();
        test_bam_reader(entry_str, 0);
    }
}

#[test]
fn bam_to_bam_singlethread() {
    for entry in glob("tests/data/read*.bam").unwrap() {
        let entry = entry.unwrap();
        println!("Analyzing {}", entry.display());
        let entry_str = entry.as_os_str().to_str().unwrap();
        test_bam_to_bam(entry_str, 0);
    }
}

#[test]
fn sam_to_bam() {
    for entry in glob("tests/data/read*.sam").unwrap() {
        let entry = entry.unwrap();
        println!("Analyzing {}", entry.display());
        let entry_str = entry.as_os_str().to_str().unwrap();
        test_sam_to_bam(entry_str);
    }
}

#[test]
fn indexed_reader_multithread() {
    for entry in glob("tests/data/iread*.bam").unwrap() {
        let entry = entry.unwrap();
        println!("Analyzing {}", entry.display());
        let entry_str = entry.as_os_str().to_str().unwrap();
        test_indexed_reader(entry_str, 2);
    }
}

#[test]
fn bam_reader_multithread() {
    for entry in glob("tests/data/read*.bam").unwrap() {
        let entry = entry.unwrap();
        println!("Analyzing {}", entry.display());
        let entry_str = entry.as_os_str().to_str().unwrap();
        test_bam_reader(entry_str, 2);
    }
}

#[test]
fn bam_to_bam_multithread() {
    for entry in glob("tests/data/read*.bam").unwrap() {
        let entry = entry.unwrap();
        println!("Analyzing {}", entry.display());
        let entry_str = entry.as_os_str().to_str().unwrap();
        test_bam_to_bam(entry_str, 2);
    }
}
