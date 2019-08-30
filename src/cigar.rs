use std::io::{self, Write};
use std::fmt::{self, Display, Formatter};

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
    pub fn from_symbol(symbol: u8) -> Operation {
        use Operation::*;
        match symbol {
            b'M' => AlnMatch,
            b'I' => Insertion,
            b'D' => Deletion,
            b'N' => Skip,
            b'S' => Soft,
            b'H' => Hard,
            b'P' => Padding,
            b'=' => SeqMatch,
            b'X' => SeqMismatch,
            _ => panic!("Unexpected cigar operation: {}", symbol as char),
        }
    }

    /// Convert [Operation](enum.Operation.html) into `u8` symbol (for example `b'M'`).
    ///
    /// To convert [Operation](enum.Operation.html) into a number `0-8` use `operation as u8`.
    pub fn to_byte(self) -> u8 {
        b"MIDNSHP=X"[self as usize]
    }

    /// Checks if the operation consumes query. For example, `M` consumes query, while `D` does not.
    pub fn consumes_query(&self) -> bool {
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
    pub fn consumes_ref(&self) -> bool {
        match self {
            Operation::AlnMatch
            | Operation::Deletion
            | Operation::Skip
            | Operation::SeqMatch
            | Operation::SeqMismatch => true,
            _ => false
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

    pub(crate) fn fill_from<R: ReadBytesExt>(&mut self, stream: &mut R, len: usize)
            -> io::Result<()> {
        unsafe {
            super::record::resize(&mut self.0, len);
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
    pub fn calculate_aligned_len(&self) -> u32 {
        self.0.iter().map(|value| match value & 0xf {
            0 | 2 | 3 | 7 | 8 => value >> 4,
            1 | 4 | 5 | 6 => 0,
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
