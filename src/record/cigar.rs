//! Cigar and operations on it.

use std::io::{self, Write};
use std::fmt::{self, Display, Formatter};
use std::slice::Iter;

use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

/// Cigar operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operation {
    AlnMatch = 0,
    Insertion = 1,
    Deletion = 2,
    Skip = 3,
    Soft = 4,
    Hard = 5,
    Padding = 6,
    SeqMatch = 7,
    SeqMismatch = 8,
}

impl Operation {
    /// Convert `u8` symbol (for example `b'M'`) into [Operation](enum.Operation.html).
    ///
    /// To convert a number `0-8` into [Operation](enum.Operation.html) use `Operation::from(number)`.
    pub fn from_symbol(symbol: u8) -> Result<Operation, String> {
        use Operation::*;
        match symbol {
            b'M' => Ok(AlnMatch),
            b'I' => Ok(Insertion),
            b'D' => Ok(Deletion),
            b'N' => Ok(Skip),
            b'S' => Ok(Soft),
            b'H' => Ok(Hard),
            b'P' => Ok(Padding),
            b'=' => Ok(SeqMatch),
            b'X' => Ok(SeqMismatch),
            _ => Err(format!("Unexpected cigar operation: {}", symbol as char)),
        }
    }

    /// Convert [Operation](enum.Operation.html) into `u8` symbol (for example `b'M'`).
    ///
    /// To convert [Operation](enum.Operation.html) into a number `0-8` use `operation as u8`.
    pub fn to_byte(self) -> u8 {
        b"MIDNSHP=X"[self as usize]
    }

    /// Checks if the operation consumes query. For example, `M` consumes query, while `D` does not.
    pub fn consumes_query(self) -> bool {
        match self {
            Operation::AlnMatch
            | Operation::Insertion
            | Operation::Soft
            | Operation::SeqMatch
            | Operation::SeqMismatch => true,
            _ => false
        }
    }

    /// Checks if the operation consumes reference.
    /// For example, `M` consumes reference, while `I` does not.
    pub fn consumes_ref(self) -> bool {
        match self {
            Operation::AlnMatch
            | Operation::Deletion
            | Operation::Skip
            | Operation::SeqMatch
            | Operation::SeqMismatch => true,
            _ => false
        }
    }

    /// Returns `true` if the operation consumes both query and reference (M, = or X).
    pub fn is_match(self) -> bool {
        match self {
            Operation::AlnMatch | Operation::SeqMatch | Operation::SeqMismatch => true,
            _ => false,
        }
    }

    /// Returns `true` if the operation consumes only reference (I or S).
    pub fn is_insertion(self) -> bool {
        match self {
            Operation::Insertion | Operation::Soft => true,
            _ => false,
        }
    }

    /// Returns `true` if the operation consumes only query (D or N).
    pub fn is_deletion(self) -> bool {
        match self {
            Operation::Deletion | Operation::Skip => true,
            _ => false,
        }
    }

    /// Returns `true` if the operation does not consume query nor reference (H or P).
    pub fn is_hard_clipping(self) -> bool {
        match self {
            Operation::Hard | Operation::Padding => true,
            _ => false,
        }
    }
}

impl Display for Operation {
    fn fmt(&self, f: &mut Formatter) -> Result<(), fmt::Error> {
        write!(f, "{}", self.to_byte() as char)
    }
}

impl From<u32> for Operation {
    fn from(value: u32) -> Operation {
        use Operation::*;
        match value {
            0 => AlnMatch,
            1 => Insertion,
            2 => Deletion,
            3 => Skip,
            4 => Soft,
            5 => Hard,
            6 => Padding,
            7 => SeqMatch,
            8 => SeqMismatch,
            _ => panic!("Unexpected cigar operation: {}", value),
        }
    }
}

/// A wrapper around raw Cigar.
pub struct Cigar(pub(crate) Vec<u32>);

impl Cigar {
    pub(crate) fn new() -> Self {
        Cigar(Vec::new())
    }

    /// Clears the contents but does not touch capacity.
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Fills Cigar from raw data.
    pub fn fill_from_raw<I: IntoIterator<Item = u32>>(&mut self, iter: I) {
        self.0.clear();
        self.0.extend(iter);
    }

    /// Fills Cigar from text. If the error occured, the cigar may be filled partly.
    pub fn fill_from_text<I: IntoIterator<Item = u8>>(&mut self, text: I) -> Result<(), String> {
        self.0.clear();
        let mut op_len = 0_u32;
        for symb in text {
            if symb >= b'0' && symb <= b'9' {
                op_len = 10 * op_len + (symb - b'0') as u32;
            } else {
                let op = Operation::from_symbol(symb)?;
                if op_len == 0 {
                    return Err("Operation length cannot be zero".to_string());
                }
                self.0.push(op_len << 4 | op as u32);
                op_len = 0;
            }
        }
        Ok(())
    }

    pub(crate) fn fill_from<R: ReadBytesExt>(&mut self, stream: &mut R, len: usize)
            -> io::Result<()> {
        unsafe {
            super::resize(&mut self.0, len);
        }
        stream.read_u32_into::<LittleEndian>(&mut self.0)?;
        Ok(())
    }

    /// Returns a pair `(length, operation)` by its index.
    pub fn at(&self, index: usize) -> (u32, Operation) {
        let v = self.0[index];
        (v >> 4, Operation::from(v & 0xf))
    }

    /// Returns an iterator over typles `(length, operation)`.
    pub fn iter(&self) -> impl Iterator<Item = (u32, Operation)> + '_ {
        (0..self.0.len()).map(move |i| self.at(i))
    }

    /// Cigar length.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns raw Cigar. Each `u32` value represents `length << 4 | operation`, where
    /// operations are encoded from 0 to 8.
    pub fn raw(&self) -> &[u32] {
        &self.0
    }

    /// Calculates reference alignment length. Consider using
    /// [Record::calculate_end](../record/struct.Record.html#method.calculate_end), as the
    /// record alignment end is stored once calculated.
    pub fn calculate_ref_len(&self) -> u32 {
        self.0.iter().map(|value| match value & 0xf {
            0 | 2 | 3 | 7 | 8 => value >> 4,
            1 | 4 | 5 | 6 => 0,
            _ => panic!("Unexpected Cigar operation: {}", value & 0xf),
        }).sum::<u32>()
    }

    /// Calculates query length.
    pub fn calculate_query_len(&self) -> u32 {
        self.0.iter().map(|value| match value & 0xf {
            0 | 1 | 4 | 7 | 8 => value >> 4,
            2 | 3 | 5 | 6 => 0,
            _ => panic!("Unexpected Cigar operation: {}", value & 0xf),
        }).sum::<u32>()
    }

    /// Shrink inner vector
    pub fn shrink_to_fit(&mut self) {
        self.0.shrink_to_fit();
    }

    /// Writes to `f` in a human readable format. Write `*` if empty.
    pub fn write_readable<W: Write>(&self, f: &mut W) -> io::Result<()> {
        if self.len() == 0 {
            return f.write_u8(b'*');
        }
        for (len, op) in self.iter() {
            write!(f, "{}", len)?;
            f.write_u8(op.to_byte())?;
        }
        Ok(())
    }

    pub(crate) fn aligned_pairs(&self, r_pos: u32) -> AlignedPairs {
        AlignedPairs {
            raw_iter: self.0.iter(),
            q_pos: 0,
            r_pos,
            remaining_len: 0,
            operation: Operation::AlnMatch,
        }
    }

    pub(crate) fn matching_pairs(&self, r_pos: u32) -> MatchingPairs {
        MatchingPairs {
            raw_iter: self.0.iter(),
            q_pos: 0,
            r_pos,
            remaining_len: 0,
        }
    }
}

/// Iterator over pairs `(Option<u32>, Option<u32>)`.
/// The first element represents a sequence index, and the second element represents a
/// reference index. If the current operation is an insertion or a deletion, the respective
/// element will be `None.`
pub struct AlignedPairs<'a> {
    raw_iter: Iter<'a, u32>,
    q_pos: u32,
    r_pos: u32,
    remaining_len: u32,
    operation: Operation,
}

impl<'a> Iterator for AlignedPairs<'a> {
    type Item = (Option<u32>, Option<u32>);

    fn next(&mut self) -> Option<Self::Item> {
        while self.remaining_len == 0 {
            let v = self.raw_iter.next()?;
            self.operation = Operation::from(v & 0xf);
            if !self.operation.is_hard_clipping() {
                self.remaining_len = v >> 4;
                break;
            }
        }
        self.remaining_len -= 1;
        let q_pos = if self.operation.consumes_query() {
            self.q_pos += 1;
            Some(self.q_pos - 1)
        } else {
            None
        };
        let r_pos = if self.operation.consumes_ref() {
            self.r_pos += 1;
            Some(self.r_pos - 1)
        } else {
            None
        };
        Some((q_pos, r_pos))
    }
}

/// Iterator over pairs `(u32, u32)`.
/// The first element represents a sequence index, and the second element represents a
/// reference index. This iterator skips insertions and deletions.
pub struct MatchingPairs<'a> {
    raw_iter: Iter<'a, u32>,
    q_pos: u32,
    r_pos: u32,
    remaining_len: u32,
}

impl<'a> Iterator for MatchingPairs<'a> {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        use Operation::*;

        while self.remaining_len == 0 {
            let v = self.raw_iter.next()?;
            let operation = Operation::from(v & 0xf);
            let op_len = v >> 4;
            match operation {
                AlnMatch | SeqMatch | SeqMismatch => self.remaining_len = v >> 4,
                Insertion | Soft => self.q_pos += op_len,
                Deletion | Skip => self.r_pos += op_len,
                _ => {},
            }
        }

        self.remaining_len -= 1;
        self.q_pos += 1;
        self.r_pos += 1;
        Some((self.q_pos - 1, self.r_pos - 1))
    }
}

impl Display for Cigar {
    fn fmt(&self, f: &mut Formatter) -> Result<(), fmt::Error> {
        if self.len() == 0 {
            return write!(f, "*");
        }
        for (len, op) in self.iter() {
            write!(f, "{}{}", len, op)?;
        }
        Ok(())
    }
}