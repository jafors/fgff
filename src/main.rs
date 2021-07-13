#[macro_use]
extern crate log;
extern crate env_logger;
extern crate fern;
extern crate csv;
extern crate serde;
#[macro_use]
extern crate clap;
extern crate bio;
extern crate itertools;
extern crate rust_htslib;
extern crate vec_map;

use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader};
use std::process;

use clap::{App, ArgMatches, SubCommand};

use bio::io::{fasta, gff};


pub mod common;
pub mod build;


pub fn run() -> Result<(), Box<dyn Error>> {
    let yaml = load_yaml!("cli.yaml");
    let build_yaml = load_yaml!("build_cli.yaml");
    let matches = App::from_yaml(yaml)
                    .version(env!("CARGO_PKG_VERSION"))
                    .subcommand(SubCommand::from_yaml(build_yaml))
                    .get_matches();

    match matches.subcommand() {
        ("build", Some(m)) => run_build(m),
        _ => Ok(())
    }
}

pub fn run_build(matches: &ArgMatches) -> Result<(), Box<dyn Error>> {
    fern::Dispatch::new()
        .format(|out, message, _| out.finish(format_args!("{}", message)))
        .level(if matches.is_present("verbose") {
            log::LevelFilter::Debug
        } else {
            log::LevelFilter::Info
        })
        .chain(std::io::stderr())
        .apply()
        .unwrap();
    let mut gff_reader = gff::Reader::new(io::stdin(), gff::GffType::GFF3);
    let mut gff_writer = gff::Writer::new(io::stdout(), gff::GffType::GFF3);
    let mut reclist_buffer = std::io::BufReader::new(File::open(matches.value_of("recordlist").unwrap()).unwrap());
    let mut tsv_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&matches.value_of("tsv").unwrap())?;
    let operation = matches.value_of("operation").unwrap();
    let mut biotypes: Vec<&str> = Vec::new();
    let tsl = matches.value_of("transcript_support_level").unwrap();
    let tpm = matches.value_of("tpm_threshold").unwrap().parse::<f64>().unwrap();

    if matches.is_present("ig") {
        biotypes.append(&mut vec!["IG_C_gene",
            "IG_D_gene",
            "IG_J_gene",
            "IG_LV_gene",
            "IG_V_gene",
            "TR_C_gene",
            "TR_J_gene",
            "TR_V_gene",
            "TR_D_gene",
            "IG_C_pseudogene",
            "IG_pseudogene",
            "IG_J_pseudogene",
            "IG_V_pseudogene",
            "TR_J_pseudogene",
            "TR_V_pseudogene",
            ]
        )
    }
    if matches.is_present("ncRNA") {
        biotypes.append(&mut vec!["Mt_rRNA",
            "Mt_tRNA",
            "miRNA",
            "misc_RNA",
            "rRNA",
            "scRNA",
            "snRNA",
            "snoRNA",
            "ribozyme",
            "sRNA",
            "scaRNA",
            "Mt_tRNA_pseudogene",
            "tRNA_pseudogene",
            "snoRNA_pseudogene",
            "snRNA_pseudogene",
            "scRNA_pseudogene",
            "rRNA_pseudogene",
            "misc_RNA_pseudogene",
            "miRNA_pseudogene",
            "lncRNA"
            ]
        )
    }
    if matches.is_present("pseudogenes") {
        biotypes.append(&mut vec!["pseudogene",
        "processed_transcript",
        "processed_pseudogene",
        "unitary_pseudogene",
        "unprocessed_pseudogene",
        "polymorphic_pseudogene",
        "transcribed_processed_pseudogene",
        "transcribed_unitary_pseudogene",
        "transcribed_unprocessed_pseudogene",
        "translated_processed_pseudogene",
        "translated_unprocessed_pseudogene",
            ]
        )
    }
    if matches.is_present("protein_coding") {
        biotypes.append(&mut vec!["protein_coding",
        "nonsense_mediated_decay",
        "non_stop_decay",
        "processed_transcript",
        "retained_intron"])
    }
    if matches.is_present("strict") {
        //println!("Strict Mode");
        biotypes = vec!["protein_coding"];
    }

    debug!("GTF");
    build::phase(
        &mut gff_reader,
        &mut gff_writer,
        &mut reclist_buffer,
        &mut tsv_reader,
        operation,
        biotypes,
        tsl,
        tpm,
    )
}

pub fn main() {
    if let Err(e) = run() {
        error!("{}", e);
        process::exit(1);
    }
}
