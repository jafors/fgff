use std::cmp::Ordering;
use std::collections::{BTreeMap};
use bio_types::strand::Strand;

use std::str;
use std::str::FromStr;

use bio::io::gff;
use serde::Deserialize;

use std::borrow::ToOwned;

#[derive(Debug, Clone)]
pub struct Gene {
    pub rec: gff::Record,
    pub transcripts: BTreeMap<String, Transcript>,
    pub interval: Interval,
    pub chrom: String,
    pub biotype: String,
}

impl Gene {
    pub fn new(rec: &gff::Record, chrom: &str, interval: Interval, biotype: &str) -> Self {
        Gene {
            rec: rec.to_owned(),
            transcripts: BTreeMap::new(),
            interval: interval,
            chrom: chrom.to_owned(),
            biotype: biotype.to_owned(),
        }
    }

    pub fn start(&self) -> u64 {
        self.interval.start
    }

    pub fn end(&self) -> u64 {
        self.interval.end
    }
}

#[derive(Debug, PartialEq)]
pub enum PhasingStrand {
    Forward,
    Reverse,
}

impl From<Strand> for PhasingStrand {
    fn from(strand: Strand) -> Self {
        let s = match strand {
            Strand::Forward => PhasingStrand::Forward,
            Strand::Reverse => PhasingStrand::Reverse,
            _ => panic!("Unsupported Strand orientation! Only Forward (+) and Reverse(-) allowed"),
        };
        s
    }
}

#[derive(Debug, Clone)]
pub struct Transcript {
    pub rec: gff::Record,
    pub exons: BTreeMap<String, Exon>,
    pub cds: BTreeMap<(u64, u64), CDS>,
    pub interval: Interval,
    pub biotype: String,
}

impl Transcript {
    pub fn new(rec: &gff::Record, interval: Interval, biotype: &str) -> Self {
        Transcript {
            rec: rec.to_owned(),
            exons: BTreeMap::new(),
            cds: BTreeMap::new(),
            interval: interval,
            biotype: biotype.to_owned(),
        }
    }
    pub fn is_coding(&self) -> bool {
        return !(self.exons.is_empty());
    }
}


#[derive(Debug, Clone)]
pub struct Exon {
    pub rec: gff::Record,
    pub interval: Interval,
}

impl Exon {
    pub fn new(rec: &gff::Record, interval: Interval) -> Self {
        Exon {
            rec: rec.to_owned(),
            interval: interval
        }
    }
}

#[derive(Debug, Clone)]
pub struct CDS {
    pub rec: gff::Record,
    pub interval: Interval,
    pub feature_type: String,
}

impl CDS {
    pub fn new(rec: &gff::Record, interval: Interval, feature_type: String) -> Self {
        CDS {
            rec: rec.to_owned(),
            interval,
            feature_type,
        }
    }
}



#[derive(Debug)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
    pub frame: u64,
}

impl Ord for Interval {
    fn cmp(&self, other: &Interval) -> Ordering {
        self.start.cmp(&other.start).then(self.end.cmp(&other.end))
    }
}

impl PartialOrd for Interval {
    fn partial_cmp(&self, other: &Interval) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Interval {
    fn eq(&self, other: &Interval) -> bool {
        self.start == other.start
    }
}

impl Eq for Interval {}

impl Clone for Interval {
    fn clone(&self) -> Interval {
        *self
    }
}

impl Copy for Interval {}

impl Interval {
    pub fn new(start: u64, end: u64, frame: &str) -> Self {
        Interval {
            start: start,
            end: end,
            frame: match frame {
                "." => 0 as u64,
                _ => u64::from_str(frame).unwrap(),
            }
        }
    }
}

#[derive(Debug, Deserialize)]
pub struct Kallisto {
    pub transcript_id: String,
    pub length: u64,
    pub	eff_length: f64,
    pub	est_counts: f64,
    pub	tpm: f64,
}



