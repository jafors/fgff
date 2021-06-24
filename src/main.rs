#[macro_use]
extern crate log;
extern crate env_logger;
extern crate fern;
#[macro_use]
extern crate clap;
extern crate bio;
extern crate itertools;
extern crate rust_htslib;
extern crate vec_map;

use std::error::Error;
use std::fs::File;
use std::io;
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

    debug!("GTF");
    build::phase(
        &mut gff_reader,
        &mut gff_writer,
    )
}

pub fn main() {
    if let Err(e) = run() {
        error!("{}", e);
        process::exit(1);
    }
}
