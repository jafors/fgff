use std::collections::{BTreeMap, VecDeque};
use std::error::Error;
use std::fs;
use std::io;

use itertools::Itertools;

use bio::io::fasta;
use bio::io::gff;

use rust_htslib::bcf;

use bio_types::strand::Strand;

use crate::common::{Gene, Exon, Interval, Transcript};


pub fn phase<G: io::Read, O: io::Write>(
    gtf_reader: &mut gff::Reader<G>,
    gff_writer: &mut gff::Writer<O>,
) -> Result<(), Box<dyn Error>> {
    let mut genes = BTreeMap::new();
    //let mut gene: Option<Gene> = None;
    let mut gene = Gene::new(&gff::Record::new(), "", Interval::new(0,0,"."), "");

    for record in gtf_reader.records() {
        debug!("New Record!");
        let record = record?;
        println!("Record: {:?}", record);
        continue;
        match record.feature_type() {
            "gene" | "ncRNA_gene" | "pseudogene" => {
                // register new Gene
                debug!("Gene found");
                println!("Gene_id {}, ID={}", record.attributes().get("gene_id").unwrap(), record.attributes().get("ID").unwrap());
                // store last gene
                if gene.biotype != "" {
                    genes.entry(gene.rec.attributes().get("gene_id").expect("missing gene_id in GFF").to_string())
                        .or_insert(gene);
                }
                gene = Gene::new(
                    &record,
                    record.seqname(),
                    Interval::new(*record.start() as u32 - 1, *record.end() as u32, record.frame()),
                    record
                        .attributes()
                        .get("biotype")
                        .expect("missing gene_biotype in GTF"),
                );
            }
            "unconfirmed_transcript" | "pseudogenic_transcript" | "mRNA" | "lnc_RNA" | "snRNA" | "rRNA" | "ncRNA" | "miRNA" | "scRNA" | "snoRNA" | "C_gene_segment" |  "J_gene_segment" | "V_gene_segment"=> {
                // register new transcript
                debug!("Transcript found");
                println!("ID: {}, Parent: {}", record.attributes().get("ID").unwrap(), record.attributes().get("Parent").unwrap());
                gene.transcripts
                    .entry(record
                        .attributes()
                        .get("transcript_id")
                        .expect("missing transcript_id attribute in GTF").to_string())
                    .or_insert(
                        Transcript::new(
                            &record,
                            Interval::new(*record.start() as u32 - 1, *record.end() as u32, record.frame()),
                            record
                                .attributes()
                                .get("biotype")
                                .expect("missing transcript_biotype in GTF"),
                    ));
                
            }
            "exon" => {
                debug!("Exon found");
                // register exon
                println!("ID: {}, Parent: {}", record.attributes().get("ID").unwrap(), record.attributes().get("Parent").unwrap());
                let key = str::replace(record.attributes().get("Parent").expect("No parent transcript for exon"), "transcript:", "");
                gene.transcripts
                    .get_mut(&key)
                    .expect("no transcript record before exon in GTF")
                    .exons
                    .entry(
                        record
                            .attributes()
                            .get("exon_id")
                            .expect("missing exon_id in GFF").to_string()
                    ).or_insert(
                        Exon::new(
                            &record,
                            Interval::new(*record.start() as u32 - 1, *record.end() as u32, record.frame()),
                    ));
            }
            "CDS" | "three_prime_UTR" | "five_prime_UTR" => {
                debug!("exonic region found");
                // register exon
                println!("ID: {}, Parent: {}", record.attributes().get("ID").unwrap(), record.attributes().get("Parent").unwrap());
                let key = str::replace(record.attributes().get("Parent").expect("No parent transcript for exon"), "transcript:", "");
                gene.transcripts
                    .get_mut(&key)
                    .expect("no transcript record before exon in GTF")
                    .exons
                    .entry(
                        format!("{}-{}", key, record.feature_type())
                    ).or_insert(
                        Exon::new(
                            &record,
                            Interval::new(*record.start() as u32 - 1, *record.end() as u32, record.frame()),
                    ));
            }
            _ => continue,
        }
    }
    // genes.entry(gene.rec.attributes().get("gene_id").expect("missing gene_id in GFF").to_string())
    // .or_insert(gene);

    // // let mut tcount = 0;
    // for (k, v) in &genes {
    //     //println!("key: {}", k);
    //     //println!("{}", genes.get(k).unwrap().biotype);
    //     let gene = genes.get(k).unwrap();
    //     if gene.biotype != "protein_coding" {
    //         continue;
    //     }
    //     gff_writer.write(&gene.rec);
    //     for (tk, tv) in &gene.transcripts {
    //         //let transcript = gene.transcripts(tk);
    //         gff_writer.write(&tv.rec);
    //         for (ek, ev) in &tv.exons {
    //             gff_writer.write(&ev.rec);
    //         }
    //     }
    // }
    // // println!("Len of all transcripts: {}", tcount);

    Ok(())
}


