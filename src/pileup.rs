//! BAM file pileup.

use std::rc::Rc;
use std::io;
use std::cmp::min;

use super::{Record, RecordReader, Region};

/// Record sequence and qualities matching a single reference position. It can be
/// - `Nucleotide(nt: u8, qual: u8)`, where `nt` is a letter, like `b'A'`, and quality does not have +33 added.
/// if the qualities are not present, `qual` is equal to `255`.
/// - `Insertion(seq: Vec<u8>, qual: Vec<u8>)`, `seq` and `qual` have the same format as in `Nucleotide(nt, qual)`.
/// First letter represents the nucleotide matching to the reference position, and the next letters represent the
/// insertion.
/// - `Deletion` - this reference position is absent from the record.
#[derive(Clone)]
pub enum PosSequence {
    Nucleotide(u8, u8),
    Insertion(Vec<u8>, Vec<u8>),
    Deletion,
}

impl PosSequence {
    /// Returns 1 for `Nucleotide`, 0 for `Deletion`, and 1 + length of the insertion for `Insertion`.
    pub fn len(&self) -> u32 {
        match self {
            PosSequence::Nucleotide(_, _) => 1,
            PosSequence::Insertion(seq, _) => 1 + seq.len() as u32,
            PosSequence::Deletion => 0,
        }
    }

    /// Returns `true` for `Nucleotide` variant, and `false` otherwise.
    pub fn is_single_nt(&self) -> bool {
        match self {
            PosSequence::Nucleotide(_, _) => true,
            _ => false,
        }
    }

    /// Returns `true` for `Insertion` variant, and `false` otherwise.
    pub fn is_insertion(&self) -> bool {
        match self {
            PosSequence::Insertion(_, _) => true,
            _ => false,
        }
    }

    /// Returns `true` for `Deletion` variant, and `false` otherwise.
    pub fn is_deletion(&self) -> bool {
        match self {
            PosSequence::Deletion => true,
            _ => false,
        }
    }

    /// Returns sequence and qualities iterator. It has [len()](#method.len) number of elements.
    pub fn seq_qual<'a>(&'a self) -> PosSequenceIterator<'a> {
        PosSequenceIterator {
            parent: self,
            i: 0,
            j: self.len(),
        }
    }
}

/// Iterator over pairs (nt, qual) for a single record and a single reference [position](struct.PosSequence.html).
#[derive(Clone)]
pub struct PosSequenceIterator<'a> {
    parent: &'a PosSequence,
    i: u32,
    j: u32,
}

impl<'a> Iterator for PosSequenceIterator<'a> {
    type Item = (u8, u8);

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < self.j {
            self.i += 1;
            match self.parent {
                PosSequence::Insertion(seq, qual) => {
                    let index = self.i as usize - 1;
                    Some((seq[index], qual[index]))
                },
                PosSequence::Nucleotide(nt, qual) => {
                    Some((*nt, *qual))
                },
                PosSequence::Deletion => unreachable!(),
            }
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = (self.j - self.i) as usize;
        (len, Some(len))
    }
}

impl<'a> DoubleEndedIterator for PosSequenceIterator<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.i < self.j {
            self.j -= 1;
            match self.parent {
                PosSequence::Insertion(seq, qual) => {
                    Some((seq[self.j as usize], qual[self.j as usize]))
                },
                PosSequence::Nucleotide(nt, qual) => {
                    Some((*nt, *qual))
                },
                PosSequence::Deletion => unreachable!(),
            }
        } else {
            None
        }
    }
}

impl<'a> ExactSizeIterator for PosSequenceIterator<'a> {}

impl<'a> std::iter::FusedIterator for PosSequenceIterator<'a> {}

/// Type of the record sequence, matching a single reference position:
/// * `Match` - single base-pair match or mismatch,
/// * `Insertion(len)` - single base-pair match followed by the insertion of length `len`,
/// * `Deletion` - this position is not present in the record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EntryType {
    Match,
    Insertion(u32),
    Deletion,
}

#[derive(Clone)]
pub struct PileupEntry {
    record: Rc<Record>,
    query_start: u32,
    query_end: u32,
}

struct PileupRecord {
    record: Rc<Record>,
    query_pos: u32,
    aln_query_end: u32,

    ref_pos: u32,
    cigar_index: usize,
    cigar_remaining: u32,
}

impl PileupRecord {
    /// Creates a new PileupRecord from a mapped `Record`.
    fn new(record: Rc<Record>) -> Self {
        let ref_pos = record.start();
        assert!(ref_pos >= 0, "Pileup Record cannot be unmapped");

        let mut cigar_index = 0;
        let mut query_pos = 0;
        let cigar_remaining = loop {
            let (len, op) = record.cigar().at(cigar_index);
            if op.consumes_ref() {
                break len;
            }
            if op.consumes_query() {
                query_pos += len;
            }
            cigar_index += 1;
        };
        assert!(cigar_index < record.cigar().len(), "CIGAR cannot contain only insertions");
        let aln_query_end = record.aligned_query_end();

        PileupRecord {
            record,

            query_pos,
            aln_query_end,

            ref_pos: ref_pos as u32,
            cigar_index,
            cigar_remaining,
        }
    }

    fn to_entry(&self) -> PileupEntry {
        let (len, op) = self.record.cigar().at(self.cigar_index);
        let query_end = if !op.consumes_query() {
            self.query_pos
        } else if self.cigar_remaining == 1 {
            let mut query_end = self.query_pos + 1;
            let mut i = self.cigar_index + 1;
            while i < self.record.cigar().len() && query_end < self.aln_query_end {
                let (len, op) = self.record.cigar().at(self.cigar_index + 1);
                if op.consumes_ref() {
                    break;
                } else if op.consumes_query() {
                    query_end += len;
                }
                i += 1;
            }
            min(query_end, self.aln_query_end)
        } else {
            self.query_pos + 1
        };

        PileupEntry {
            record: Rc::clone(&self.record),
            query_start: self.query_pos,
            query_end,
        }
    }

    /// Moves the record to the next reference position. Returns false if reached the end of the alignment.
    fn move_forward(&mut self) -> bool {
        let (len, op) = self.record.cigar().at(self.cigar_index);
        self.cigar_remaining -= 1;
        if op.consumes_ref() {
            self.ref_pos += 1;
        }
        if op.consumes_query() {
            self.query_pos += 1;
        }

        while self.cigar_remaining == 0 {
            self.cigar_index += 1;
            if self.cigar_index == self.record.cigar().len() || self.query_pos >= self.aln_query_end {
                return false;
            }
            let (len, op) = self.record.cigar().at(self.cigar_index);
            if op.consumes_ref() {
                self.cigar_remaining = len;
            } else if op.consumes_query() {
                self.query_pos += len;
            }
        }
        assert!(self.query_pos < self.aln_query_end);
        true
    }
}

/// [Pileup](struct.Pileup.html) builder.
pub struct PileupBuilder<'r, R: RecordReader> {
    reader: &'r mut R,
    read_filter: Option<Box<dyn Fn(&Record) -> bool>>,
}

impl<'r, R: RecordReader> PileupBuilder<'r, R> {
    pub(crate) fn new(reader: &'r mut R) -> Self {
        Self {
            reader,
            read_filter: None,
        }
    }

    /// Sets the read filter. The function panics if the filter was previously set.
    pub fn read_filter(&mut self, filter: impl 'static + Fn(&Record) -> bool) -> &mut Self {
        assert!(self.read_filter.is_none());
        self.read_filter = Some(Box::new(filter));
        self
    }
}

pub struct Pileup<'r, R: RecordReader> {
    reader: &'r mut R,
    read_filter: Option<Box<dyn Fn(&Record) -> bool>>,
    records: Vec<PileupRecord>,
    error: Option<io::Error>,

    ref_id: u32,
    ref_pos: u32,
    last_ref_id: u32,
    last_ref_pos: u32,
}

impl<'r, R: RecordReader> Pileup<'r, R> {
    fn new(reader: &'r mut R, read_filter: Option<Box<dyn Fn(&Record) -> bool>>) -> Self {
        let mut res = Pileup {
            reader,
            read_filter,
            records: Vec::new(),
            error: None,

            ref_id: 0,
            ref_pos: 0,
            last_ref_id: 0,
            last_ref_pos: 0,
        };
        res.read_next();
        res.ref_id = res.last_ref_id;
        res.ref_pos = res.last_ref_pos;
        res
    }

    fn record_passes(&self, record: &Record) -> bool {
        if !record.flag().is_mapped() {
            return false;
        }
        assert!(record.ref_id() >= 0 && record.start() >= 0);
        if let Some(filter) = &self.read_filter {
            filter(record)
        } else {
            true
        }
    }

    fn read_next(&mut self) {
        if self.last_ref_id == std::u32::MAX || self.error.is_some() {
            return;
        }
        loop {
            match self.reader.next() {
                None => self.last_ref_id = std::u32::MAX,
                Some(Ok(record)) => {
                    if !self.record_passes(&record) {
                        continue;
                    }
                    let rec_ref_id = record.ref_id() as u32;
                    let rec_start = record.start() as u32;
                    if rec_ref_id < self.last_ref_id
                            || (rec_ref_id == self.last_ref_id && rec_start < self.last_ref_pos) {
                        self.error = Some(io::Error::new(io::ErrorKind::InvalidData, "Input file is unsorted"));
                        self.last_ref_id = std::u32::MAX;
                    }
                    self.last_ref_id = rec_ref_id;
                    self.last_ref_pos = rec_start;
                    self.records.push(PileupRecord::new(Rc::new(record)));
                },
                Some(Err(e)) => {
                    self.error = Some(e);
                    self.last_ref_id = std::u32::MAX;
                },
            }
            return;
        }
    }
}

impl<'r, R: RecordReader> Iterator for Pileup<'r, R> {
    type Item = io::Result<PileupColumn>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(e) = &self.error {
            // return Some(Err(e));
            // TODO FIX
            return None;
        }

        let mut new_ref_id = std::u32::MAX;
        let mut new_ref_pos = std::u32::MAX;
        while new_ref_id == std::u32::MAX && (!self.records.is_empty() || self.last_ref_id < std::u32::MAX) {
            for pil_rec in self.records.iter() {
                let rec_ref_id = pil_rec.record.ref_id() as u32;
                if rec_ref_id < new_ref_id {
                    new_ref_id = rec_ref_id;
                    new_ref_pos = pil_rec.ref_pos;
                } else if rec_ref_id == new_ref_id {
                    new_ref_pos = min(new_ref_pos, pil_rec.ref_pos);
                }
            }

            while self.last_ref_id < std::u32::MAX
                    // new_ref_id == last_ref_id or new_ref_id == u32::MAX, same with pos.
                    && self.last_ref_id <= new_ref_id && self.last_ref_pos <= new_ref_pos {
                self.read_next();
            }
            if let Some(e) = &self.error {
                //return Some(Err(e.clone()));
                // TODO FIX
                return None;
            }
        }

        let mut entries = Vec::new();
        let mut i = self.records.len();
        for i in (0..self.records.len()).rev() {
            let mut pil_rec = &mut self.records[i];
            let rec_ref_id = pil_rec.record.ref_id() as u32;
            if rec_ref_id == new_ref_id && pil_rec.ref_pos == new_ref_pos {
                entries.push(pil_rec.to_entry());
                if !pil_rec.move_forward() {
                    std::mem::drop(pil_rec);
                    self.records.swap_remove(i);
                }
            } else {
                assert!(rec_ref_id > new_ref_id || pil_rec.ref_pos > new_ref_pos,
                    "Record is to the left of the new pileup position");
            }
        }

        if entries.is_empty() {
            None
        } else {
            Some(Ok(PileupColumn {
                entries,
                ref_id: new_ref_id,
                ref_pos: new_ref_pos,
            }))
        }
    }
}

pub struct PileupColumn {
    entries: Vec<PileupEntry>,
    ref_id: u32,
    ref_pos: u32,
}
