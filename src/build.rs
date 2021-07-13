use std::collections::{BTreeMap, VecDeque};
use std::error::Error;
use std::fs;
use std::io::{self, BufRead, BufReader};
use std::iter::FromIterator;

use itertools::Itertools;

use bio::io::fasta;
use bio::io::gff;

use rust_htslib::bcf;

use bio_types::strand::Strand;

use crate::common::{Gene, Exon, Interval, Transcript, CDS, Kallisto};


pub fn kill(kids: Vec<String>, mut genes: BTreeMap<String, Gene>) -> Result<BTreeMap<String, Gene>, Box<dyn Error>> {
    for kid in kids {
        if genes.contains_key(&kid) {
            genes.remove(&kid);
        }
        for (_gid, g) in &mut genes {

            if g.transcripts.contains_key(&kid) {
                g.transcripts.remove(&kid);
                break;
            }
            let mut keep_transcripts = BTreeMap::new();
            for (tid, t) in &mut g.transcripts {
                if t.exons.contains_key(&kid) {
                    t.exons.remove(&kid);
                }
                if !t.exons.is_empty() {
                    keep_transcripts.insert(tid.clone(), t.clone());
                }
            }
            g.transcripts = keep_transcripts;
        }
    }
    Ok(genes)
}

pub fn keep(kids: Vec<String>, mut genes: BTreeMap<String, Gene>) -> Result<BTreeMap<String, Gene>, Box<dyn Error>> {
    let _out: BTreeMap<String, Gene> = BTreeMap::new();
    if kids[0].contains("ENSG") {
        genes.retain(|k, _| kids.contains(k));
        Ok(genes)
    }
    else {
        for (gid, g) in &mut genes {
            if kids[0].contains("ENST") {
                g.transcripts.retain(|k, _| kids.contains(k));
            }
            let mut empty_transcripts = Vec::new();
            let mut keep_transcripts = BTreeMap::new();
            for (tid, t) in &mut g.transcripts {
                if kids[0].contains("ENSE") {
                    t.exons.retain(|k, _| kids.contains(k));
                }
                if t.exons.is_empty() {
                    empty_transcripts.push(tid);
                }
                else {
                    keep_transcripts.insert(tid.clone(), t.clone());
                }
            }
            g.transcripts = keep_transcripts;
        }

        Ok(genes)
    }
}


pub fn filter<F: io::Read + io::Seek,>(
    tsv_reader: &mut csv::Reader<F>,
    threshold: f64,
    mut genes: BTreeMap<String, Gene>,
    mut parents: BTreeMap<String, String>,
) -> Result<BTreeMap<String, Gene>, Box<dyn Error>> {
    let mut out_genes = BTreeMap::new();
    for record in tsv_reader.records() {
        let record = record?;
        let row: Kallisto = record.deserialize(None)?;
        let v: Vec<&str> = row.transcript_id.split('.').collect();
        let tid = v[0];
        if row.tpm >= threshold {
            if !parents.contains_key(tid) {
                continue
            };
            let gid = parents.get(tid).unwrap();
            let g = genes.get_mut(gid).unwrap();
            let t = g.transcripts.get(tid).unwrap();
            let mut new_gene = g.clone();
            new_gene.transcripts.clear();
            out_genes.entry(gid.clone()).or_insert_with(|| new_gene).transcripts.insert(tid.to_string(), t.clone());
        }
    }
    Ok(out_genes)
}

pub fn parse_tsl (tsl: &str) -> u32 {
    let mut tsl_threshold = match tsl.parse::<u32>() {
        Ok(t) => t,
        Err(_e) => 6,
    };
    if tsl_threshold > 6 {
        tsl_threshold = 6;
    }
    tsl_threshold
}

pub fn phase<G: io::Read, O: io::Write, F: io::Read + io::Seek,>(
    gtf_reader: &mut gff::Reader<G>,
    gff_writer: &mut gff::Writer<O>,
    reclist_buffer: &mut BufReader<F>,
    tsv_reader: &mut csv::Reader<F>,
    operation: &str,
    biotypes: Vec<&str>,
    tsl: &str,
    tpm: f64,
) -> Result<(), Box<dyn Error>> {
    let mut genes = BTreeMap::new();
    let mut gene = Gene::new(&gff::Record::new(), "", Interval::new(0,0,"."), "");
    let mut parents = BTreeMap::new();
    // get tsl threshold - 1-5 or NA (represented as 6)
    let tsl_threshold = parse_tsl(tsl);
    let reclist = reclist_buffer.lines().map(|l| l.unwrap());
    for record in gtf_reader.records() {
        debug!("New Record!");
        let record = record?;
        //debug!("Record: {:?}", record);
        match record.feature_type() {
            "gene" | "ncRNA_gene" | "pseudogene" => {
                // register new Gene
                debug!("Gene found");
                debug!("Gene_id {}, ID={}", record.attributes().get("gene_id").unwrap(), record.attributes().get("ID").unwrap());

                let genetype = match record.attributes().contains_key("biotype") {
                    true => "biotype",
                    false => "gene_type"
                };
                if !biotypes.is_empty() && !biotypes.contains(&record.attributes().get(genetype).unwrap().as_str()) {
                    continue;
                }
                gene = Gene::new(
                    &record,
                    record.seqname(),
                    Interval::new(*record.start() as u64 - 1, *record.end() as u64, record.frame()),
                    record
                        .attributes()
                        .get(genetype)
                        .expect("missing gene_biotype in GTF"),
                );
                

                genes.entry(gene.rec.attributes().get("gene_id").expect("missing gene_id in GFF").to_string())
                    .or_insert(gene);


            }
            "unconfirmed_transcript" | "pseudogenic_transcript" | "mRNA" | "lnc_RNA" | "snRNA" | "rRNA" | "ncRNA" | "miRNA" | "scRNA" | "snoRNA" | "tRNA" | "D_gene_segment" | "C_gene_segment" |  "J_gene_segment" | "V_gene_segment" => {
                // register new transcript
                debug!("Transcript found");
                let transcript_type = match record.attributes().contains_key("biotype") {
                    true => "biotype",
                    false => "transcript_type"
                };
                
                //debug!("ID: {}, Parent: {}", record.attributes().get("ID").unwrap(), record.attributes().get("Parent").unwrap());
                let transcript_id = record
                    .attributes()
                    .get("transcript_id")
                    .expect("missing transcript_id attribute in GTF").to_string();
                
                let gene_id = str::replace(record.attributes().get("Parent").expect("No parent gene for transcript"), "gene:", "");
                if !biotypes.is_empty() && (!biotypes.contains(&record.attributes().get(transcript_type).unwrap().as_str()) || !genes.contains_key(&gene_id)) {
                    continue;
                }
                if record.attributes().contains_key("transcript_support_level") {
                    if parse_tsl(record.attributes().get("transcript_support_level").expect(format!("No TSL for this transcript {}", transcript_id).as_str())) > tsl_threshold {
                        continue;
                    }
                }
                let gene = genes.get_mut(&gene_id).expect(&format!("no gene with this ID {}", gene_id));
                parents.insert(transcript_id, gene_id);
                gene.transcripts
                    .entry(record
                        .attributes()
                        .get("transcript_id")
                        .expect("missing transcript_id attribute in GTF").to_string())
                    .or_insert(
                        Transcript::new(
                            &record,
                            Interval::new(*record.start() as u64 - 1, *record.end() as u64, record.frame()),
                            record
                                .attributes()
                                .get(transcript_type)
                                .expect("missing transcript_biotype in GTF"),
                    ));
                
            }
            "exon" => {
                debug!("Exon found");
                debug!("{:?}", parents);
                // register exon
                //debug!("ID: {}, Parent: {}", record.attributes().get("exon_id").unwrap(), record.attributes().get("Parent").unwrap());
                let key = str::replace(record.attributes().get("Parent").expect("No parent transcript for exon"), "transcript:", "");
                if (!biotypes.is_empty() || tsl_threshold < 6) && !parents.contains_key(&key) {
                    continue;
                }
                let gene_id = parents.get(&key).expect(&format!("no transcript with this ID {}", key));
                let gene = genes.get_mut(gene_id).expect(&format!("no gene with this ID {}", gene_id));
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
                            Interval::new(*record.start() as u64 - 1, *record.end() as u64, record.frame()),
                    ));
            }
            "CDS" => {
                debug!("exonic region found");
                // register CDS
                //debug!("ID: {}, Parent: {}", record.attributes().get("ID").unwrap(), record.attributes().get("Parent").unwrap());
                let key = str::replace(record.attributes().get("Parent").expect("No parent transcript for exon"), "transcript:", "");
                if (!biotypes.is_empty() || tsl_threshold < 6) && !parents.contains_key(&key) {
                    continue;
                }
                let gene_id = parents.get(&key).expect(&format!("no transcript with this ID {}", key));
                let gene = genes.get_mut(gene_id).expect(&format!("no gene with this ID {}", gene_id));
                gene.transcripts
                    .get_mut(&key)
                    .expect("no transcript record before exon in GTF")
                    .cds
                    .entry(
                        (*record.start() as u64 - 1, *record.end() as u64)
                    ).or_insert(
                        CDS::new(
                            &record,
                            Interval::new(*record.start() as u64 - 1, *record.end() as u64, record.frame()),
                            String::from(&(*record.feature_type())),
                    ));
            }
            "three_prime_UTR" | "five_prime_UTR" => {
                debug!("UTR region found");
                // register UTR
                //debug!("ID: {}, Parent: {}", record.attributes().get("ID").unwrap(), record.attributes().get("Parent").unwrap());
                let key = str::replace(record.attributes().get("Parent").expect("No parent transcript for exon"), "transcript:", "");
                if (!biotypes.is_empty() || tsl_threshold < 6) && !parents.contains_key(&key) {
                    continue;
                }
                let gene_id = parents.get(&key).expect(&format!("no transcript with this ID {}", key));
                let gene = genes.get_mut(gene_id).expect(&format!("no gene with this ID {}", gene_id));
                gene.transcripts
                    .get_mut(&key)
                    .expect("no transcript record before exon in GTF")
                    .cds
                    .entry(
                        (*record.start() as u64 - 1, *record.end() as u64)
                    ).or_insert(
                        CDS::new(
                            &record,
                            Interval::new(*record.start() as u64 - 1, *record.end() as u64, record.frame()),
                            String::from(&(*record.feature_type())),
                    ));
            }
            _ => continue,
        }
    }

    let tpm_threshold = tpm;

    let out_genes = filter(tsv_reader, tpm_threshold, genes, parents).unwrap();

    let mut out_vec = Vec::from_iter(out_genes);

    out_vec.sort_by(|(_, a), (_, b)| a.chrom.cmp(&b.chrom).then(a.interval.cmp(&b.interval)));

    for (_k, v) in &mut out_vec {
        //debug!("key: {}", k);
        //debug!("{}", genes.get(k).unwrap().biotype);
        if v.transcripts.is_empty() {
            continue;
        }
        // if gene.biotype != "protein_coding" {
        //     continue;
        // }
        // if gene.transcripts.is_empty() {
        //     out_genes.remove(&k);
        // }
        gff_writer.write(&v.rec)?;
        let mut out_t_vec = Vec::from_iter(v.transcripts.clone());
        out_t_vec.sort_by(|(_, a), (_, b)| a.interval.cmp(&b.interval));
        for (_tk, tv) in &mut out_t_vec {
            //let transcript = gene.transcripts(tk);
            if tv.exons.is_empty() {
                continue;
            }
            let mut out_e_vec  = Vec::from_iter(tv.exons.clone());
            out_e_vec.sort_by(|(_, a), (_, b)| a.interval.cmp(&b.interval));
            gff_writer.write(&tv.rec)?;
            for (_ek, ev) in &out_e_vec {
                gff_writer.write(&ev.rec)?;
            }
            for (_ck, cv) in &tv.cds {
                gff_writer.write(&cv.rec)?;
            }
        }
    }
    Ok(())
}