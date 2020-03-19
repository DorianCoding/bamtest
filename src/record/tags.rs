//! Record tags and operations on them.

use std::io::{self, Read, Write};
use std::io::ErrorKind::InvalidData;
use std::mem;

use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

use super::resize;

/// Enum that represents tag type for the cases when a tag contains integer.
///
/// Possible values are `I8` (`c`), `U8` (`C`), `I16` (`s`), `U16` (`S`), `I32` (`i`) and `U32` (`I`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum IntegerType {
    I8,
    U8,
    I16,
    U16,
    I32,
    U32,
}

impl IntegerType {
    /// Returns a letter that represents the integer type. For example, `i8` corresponds to `c`.
    pub fn letter(self) -> u8 {
        use IntegerType::*;
        match self {
            I8 => b'c',
            U8 => b'C',
            I16 => b's',
            U16 => b'S',
            I32 => b'i',
            U32 => b'I',
        }
    }

    /// Returns IntegerType from a letter, such as `c`.
    pub fn from_letter(ty: u8) -> Option<Self> {
        use IntegerType::*;
        match ty {
            b'c' => Some(I8),
            b'C' => Some(U8),
            b's' => Some(I16),
            b'S' => Some(U16),
            b'i' => Some(I32),
            b'I' => Some(U32),
            _ => None,
        }
    }

    pub fn size_of(self) -> usize {
        use IntegerType::*;
        match self {
            I8 | U8 => 1,
            I16 | U16 => 2,
            I32 | U32 => 4,
        }
    }

    fn parse_raw(self, mut raw_tag: &[u8]) -> i64 {
        use IntegerType::*;
        match self {
            I8 => unsafe { mem::transmute::<u8, i8>(raw_tag[0]) as i64 },
            U8 => raw_tag[0] as i64,
            I16 => raw_tag.read_i16::<LittleEndian>()
                .expect("Failed to extract tag from a raw representation.") as i64,
            U16 => raw_tag.read_u16::<LittleEndian>()
                .expect("Failed to extract tag from a raw representation.") as i64,
            I32 => raw_tag.read_i32::<LittleEndian>()
                .expect("Failed to extract tag from a raw representation.") as i64,
            U32 => raw_tag.read_u32::<LittleEndian>()
                .expect("Failed to extract tag from a raw representation.") as i64,
        }
    }
}

fn parse_float(mut raw_tag: &[u8]) -> f32 {
    raw_tag.read_f32::<LittleEndian>()
        .expect("Failed to extract tag from a raw representation.")
}

/// Enum that represents tag type for `String` and `Hex` types (`Z` and `H`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum StringType {
    String,
    Hex,
}

impl StringType {
    /// Returns a letter that represents the string type.
    pub fn letter(&self) -> u8 {
        match self {
            StringType::String => b'Z',
            StringType::Hex => b'H',
        }
    }

    /// Returns StringType from letters `Z` and `H`.
    pub fn from_letter(ty: u8) -> Option<Self> {
        match ty {
            b'Z' => Some(StringType::String),
            b'H' => Some(StringType::Hex),
            _ => None,
        }
    }
}

/// Wrapper around raw integer array stored in a tag.
pub struct IntArrayView<'a> {
    raw: &'a [u8],
    int_type: IntegerType,
}

impl<'a> IntArrayView<'a> {
    /// Get the type of the inner array.
    pub fn int_type(&self) -> IntegerType {
        self.int_type
    }

    pub fn len(&self) -> usize {
        self.raw.len() / self.int_type.size_of()
    }

    /// Returns an element at the `index`. Returns `i64` to include both `i32` and `u32`.
    pub fn at(&self, index: usize) -> i64 {
        let start = self.int_type.size_of() * index;
        let end = start + self.int_type.size_of();
        assert!(end <= self.raw.len(), "Index out of bounds: index {}, len {}", index, self.len());
        self.int_type.parse_raw(&self.raw[start..end])
    }

    /// Returns iterator over values (converted into `i64`).
    pub fn iter<'b: 'a>(&'b self) -> IntArrayViewIter<'b> {
        IntArrayViewIter {
            chunks: self.raw.chunks(self.int_type.size_of()),
            int_type: self.int_type,
        }
    }

    /// Returns raw array.
    pub fn raw(&self) -> &[u8] {
        self.raw
    }
}

/// Double-ended iterator over [IntArrayView](struct.IntArrayView.html)
/// values converted into `i64`.
#[derive(Clone)]
pub struct IntArrayViewIter<'a> {
    chunks: std::slice::Chunks<'a, u8>,
    int_type: IntegerType,
}

impl<'a> Iterator for IntArrayViewIter<'a> {
    type Item = i64;

    fn next(&mut self) -> Option<Self::Item> {
        self.chunks.next().map(|chunk| self.int_type.parse_raw(chunk))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.chunks.size_hint()
    }
}

impl<'a> DoubleEndedIterator for IntArrayViewIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.chunks.next_back().map(|chunk| self.int_type.parse_raw(chunk))
    }
}

impl<'a> ExactSizeIterator for IntArrayViewIter<'a> {}
impl<'a> std::iter::FusedIterator for IntArrayViewIter<'a> {}

/// Wrapper around raw float array stored in a tag.
#[derive(Clone)]
pub struct FloatArrayView<'a> {
    raw: &'a [u8],
}

impl<'a> FloatArrayView<'a> {
    pub fn len(&self) -> usize {
        self.raw.len() / 4
    }

    /// Returns an element at the `index`.
    pub fn at(&self, index: usize) -> f32 {
        let start = 4 * index;
        let end = start + 4;
        assert!(end <= self.raw.len(), "Index out of bounds: index {}, len {}", index, self.len());
        parse_float(&self.raw[start..end])
    }

    pub fn iter<'b: 'a>(&'b self) -> FloatArrayViewIter<'b> {
        FloatArrayViewIter {
            chunks: self.raw.chunks(4),
        }
    }

    /// Returns raw array.
    pub fn raw(&self) -> &[u8] {
        self.raw
    }
}

/// Double-ended iterator over [FloatArrayView](struct.FloatArrayView.html) values.
pub struct FloatArrayViewIter<'a> {
    chunks: std::slice::Chunks<'a, u8>,
}

impl<'a> Iterator for FloatArrayViewIter<'a> {
    type Item = f32;

    fn next(&mut self) -> Option<Self::Item> {
        self.chunks.next().map(|chunk| parse_float(chunk))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.chunks.size_hint()
    }
}

impl<'a> DoubleEndedIterator for FloatArrayViewIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.chunks.next_back().map(|chunk| parse_float(chunk))
    }
}

impl<'a> ExactSizeIterator for FloatArrayViewIter<'a> {}
impl<'a> std::iter::FusedIterator for FloatArrayViewIter<'a> {}


/// Enum with all possible tag values.
///
/// # Variants
/// * `Char` - contains a one-byte character,
/// * `Int(i64, IntegerType)` - contains an integer in `i64` format to be able to store both
/// `i32` and `u32`. Enum [IntegerType](enum.IntegerType.html) specifies initial integer size.
/// * `Float` - contains a float,
/// * `String(&[u8], StringType)` - contains a string as bytes and a [type](enum.StringType.html)
/// of the string - `String` or `Hex`.
/// * `IntArray` - contains a [view](struct.IntArrayView.html) over an integer array, which
/// allows to get an element at a specific index and iterate over all values,
/// * `FloatArray` - contains a [view](struct.FloatArrayView.html) over a float array.
pub enum TagValue<'a> {
    Char(u8),
    Int(i64, IntegerType),
    Float(f32),
    String(&'a [u8], StringType),
    IntArray(IntArrayView<'a>),
    FloatArray(FloatArrayView<'a>),
}

impl<'a> TagValue<'a> {
    /// Get [TagValue](enum.TagValue.html) from a raw representation.
    /// This function expects the raw tag to have a correct length.
    fn from_raw(ty: u8, raw_tag: &'a [u8]) -> TagValue<'a> {
        use TagValue::*;
        if let Some(int_type) = IntegerType::from_letter(ty) {
            return Int(int_type.parse_raw(raw_tag), int_type);
        }
        if let Some(str_type) = StringType::from_letter(ty) {
            return String(&raw_tag[..raw_tag.len() - 1], str_type);
        }

        match ty {
            b'A' => Char(raw_tag[0]),
            b'f' => Float(parse_float(raw_tag)),
            b'B' => {
                let arr_ty = raw_tag[0];
                if arr_ty == b'f' {
                    return FloatArray(FloatArrayView {
                        raw: &raw_tag[5..],
                    });
                }
                if let Some(int_type) = IntegerType::from_letter(arr_ty) {
                    return IntArray(IntArrayView {
                        raw: &raw_tag[5..],
                        int_type,
                    });
                }
                panic!("Unexpected tag array type: {}", arr_ty as char);
            }
            _ => panic!("Unexpected tag type: {}", ty as char),
        }
    }

    /// Write the tag value in a sam format.
    pub fn write_sam<W: Write>(&self, f: &mut W) -> io::Result<()> {
        use TagValue::*;
        match self {
            Char(value) => f.write_all(&[b'A', b':', *value]),
            Int(value, _) => write!(f, "i:{}", *value),
            Float(value) => write!(f, "f:{}", *value),
            String(u8_slice, str_type) => {
                f.write_all(&[str_type.letter(), b':'])?;
                f.write_all(u8_slice)
            },
            IntArray(arr_view) => {
                f.write_all(b"B:i")?;
                for value in arr_view.iter() {
                    write!(f, ",{}", value)?;
                }
                Ok(())
            },
            FloatArray(arr_view) => {
                f.write_all(b"B:f")?;
                for value in arr_view.iter() {
                    write!(f, ",{}", value)?;
                }
                Ok(())
            },
        }
    }
}

/// Wrapper around `&u8`, used to specify tag type as _Hex_, not `String` and not `&[u8]`.
/// Should have an even number of characters.
pub struct Hex<'a>(pub &'a [u8]);

mod private {
    /// Sealed trait, designed to forbid new implementations of `WriteValue`.
    pub trait Sealed {}
    impl Sealed for char {}
    impl Sealed for i8 {}
    impl Sealed for u8 {}
    impl Sealed for i16 {}
    impl Sealed for u16 {}
    impl Sealed for i32 {}
    impl Sealed for u32 {}
    impl Sealed for f32 {}
    
    impl Sealed for &str {}
    impl<'a> Sealed for super::Hex<'a> {}
    
    impl Sealed for &[i8] {}
    impl Sealed for &[u8] {}
    impl Sealed for &[i16] {}
    impl Sealed for &[u16] {}
    impl Sealed for &[i32] {}
    impl Sealed for &[u32] {}
    impl Sealed for &[f32] {}
}

/// A trait for writing tag values. Implemented for:
/// * `char`, but will return error, if char takes more than one byte,
/// * numeric values `i8`, `u8`, `i16`, `u16`, `i32`, `u32`, `f32`,
/// * numeric slices `&[i8]`, `&[u8]`, `&[i16]`, ..., `&[f32]`,
/// * string slice `&str`. You can use `std::str::from_utf8_unchecked` to convert `&[u8]` to `&str`.
/// We use `&str` here to distinguish between string and int array types.
/// * [hex wrapper](struct.Hex.html) over `&[u8]`,
///
/// String and Hex values cannot contain null symbols, even at the end.
///
/// The trait cannot be implemented for new types.
pub trait WriteValue: private::Sealed {
    /// Returns the number of written bytes (does not include name).
    fn write<W: Write>(&self, f: &mut W) -> io::Result<usize>;
}

impl WriteValue for char {
    fn write<W: Write>(&self, f: &mut W) -> io::Result<usize> {
        if self.len_utf8() != 1 {
            return Err(io::Error::new(InvalidData,
                "Cannot write tag value: char value takes more than one byte"));
        }
        f.write_all(&[b'A', *self as u8])?;
        Ok(2)
    }
}

impl WriteValue for i8 {
    fn write<W: Write>(&self, f: &mut W) -> io::Result<usize> {
        f.write_all(&[b'c', unsafe { std::mem::transmute::<i8, u8>(*self) }])?;
        Ok(2)
    }
}

impl WriteValue for u8 {
    fn write<W: Write>(&self, f: &mut W) -> io::Result<usize> {
        f.write_all(&[b'C', *self])?;
        Ok(2)
    }
}

macro_rules! write_value_prim {
    ($name:ty, $letter:expr, $fun:ident) => {
        impl WriteValue for $name {
            fn write<W: Write>(&self, f: &mut W) -> io::Result<usize> {
                f.write_u8($letter)?;
                f.$fun::<LittleEndian>(*self)?;
                Ok(1 + mem::size_of::<$name>())
            }
        }
    }
}

write_value_prim!(i16, b's', write_i16);
write_value_prim!(u16, b'S', write_u16);
write_value_prim!(i32, b'i', write_i32);
write_value_prim!(u32, b'I', write_u32);
write_value_prim!(f32, b'f', write_f32);

impl WriteValue for &str {
    fn write<W: Write>(&self, f: &mut W) -> io::Result<usize> {
        if self.bytes().any(|ch| ch == 0) {
            return Err(io::Error::new(InvalidData, "Cannot write tag value: String contains null"));
        }
        f.write_u8(b'Z')?;
        f.write_all(self.as_bytes())?;
        f.write_u8(0)?;
        Ok(2 + self.len())
    }
}

impl<'a> WriteValue for Hex<'a> {
    fn write<W: Write>(&self, f: &mut W) -> io::Result<usize> {
        if self.0.len() % 2 != 0 {
            return Err(io::Error::new(InvalidData,
                "Cannot write tag value: Hex string contains an odd number of bytes"));
        }
        if self.0.iter().any(|&ch| ch == 0) {
            return Err(io::Error::new(InvalidData, "Cannot write tag value: Hex contains null"));
        }
        f.write_u8(b'H')?;
        f.write_all(self.0)?;
        f.write_u8(0)?;
        Ok(2 + self.0.len())
    }
}

impl WriteValue for &[i8] {
    fn write<W: Write>(&self, f: &mut W) -> io::Result<usize> {
        f.write_all(b"Bc")?;
        f.write_i32::<LittleEndian>(self.len() as i32)?;
        unsafe {
            f.write_all(std::slice::from_raw_parts(self.as_ptr() as *const u8, self.len()))?;
        }
        Ok(6 + self.len())
    }
}

impl WriteValue for &[u8] {
    fn write<W: Write>(&self, f: &mut W) -> io::Result<usize> {
        f.write_all(b"BC")?;
        f.write_i32::<LittleEndian>(self.len() as i32)?;
        f.write_all(self)?;
        Ok(6 + self.len())
    }
}

macro_rules! write_value_array {
    ($name:ty, $letter:expr, $fun:ident) => {
        impl WriteValue for &[$name] {
            fn write<W: Write>(&self, f: &mut W) -> io::Result<usize> {
                f.write_all(&[b'B', $letter])?;
                f.write_i32::<LittleEndian>(self.len() as i32)?;
                for v in self.iter() {
                    f.$fun::<LittleEndian>(*v)?;
                }
                Ok(6 + self.len() * mem::size_of::<$name>())
            }
        }
    }
}

write_value_array!(i16, b's', write_i16);
write_value_array!(u16, b'S', write_u16);
write_value_array!(i32, b'i', write_i32);
write_value_array!(u32, b'I', write_u32);
write_value_array!(f32, b'f', write_f32);

pub trait PushNum: Copy + Sized {
    #[doc(hidden)]
    fn push_individually(self, raw: &mut Vec<u8>) -> usize;

    #[doc(hidden)]
    fn push_array(arr: &[Self], raw: &mut Vec<u8>) -> usize;
}

macro_rules! push_num {
    (
        $push_fun:ident($type:ty) { $write_fun:ident($letter:expr) }
    ) => {
        #[inline]
        fn $push_fun(value: $type, raw: &mut Vec<u8>) -> usize {
            raw.push($letter);
            raw.$write_fun::<LittleEndian>(value).unwrap();
            1 + mem::size_of::<$type>()
        }
    }
}

#[inline]
fn push_i8(value: i8, raw: &mut Vec<u8>) -> usize {
    raw.push(b'c');
    raw.push(unsafe { std::mem::transmute::<i8, u8>(value) });
    2
}

#[inline]
fn push_u8(value: u8, raw: &mut Vec<u8>) -> usize {
    raw.push(b'C');
    raw.push(value);
    2
}

push_num!{
    push_i16(i16) { write_i16(b's') }
}

push_num!{
    push_u16(u16) { write_u16(b'S') }
}

push_num!{
    push_i32(i32) { write_i32(b'i') }
}

push_num!{
    push_u32(u32) { write_u32(b'I') }
}

impl PushNum for i8 {
    fn push_individually(self, raw: &mut Vec<u8>) -> usize {
        push_i8(self, raw)
    }

    fn push_array(arr: &[Self], raw: &mut Vec<u8>) -> usize {
        raw.push(b'B');
        raw.push(b'c');
        raw.write_i32::<LittleEndian>(arr.len() as i32).unwrap();
        unsafe {
            raw.write_all(std::slice::from_raw_parts(arr.as_ptr() as *const u8, arr.len())).unwrap();
        }
        6 + arr.len()
    }
}

impl PushNum for u8 {
    fn push_individually(self, raw: &mut Vec<u8>) -> usize {
        push_u8(self, raw)
    }

    fn push_array(arr: &[Self], raw: &mut Vec<u8>) -> usize {
        raw.push(b'B');
        raw.push(b'C');
        raw.write_i32::<LittleEndian>(arr.len() as i32).unwrap();
        raw.write_all(arr).unwrap();
        6 + arr.len()
    }
}

macro_rules! push_num_array {
    ($fun:ident, $letter:expr) => {
        fn push_array(arr: &[Self], raw: &mut Vec<u8>) -> usize {
            raw.push(b'B');
            raw.push($letter);
            raw.write_i32::<LittleEndian>(arr.len() as i32).unwrap();
            for v in arr.iter() {
                raw.$fun::<LittleEndian>(*v).unwrap();
            }
            6 + arr.len() * mem::size_of::<Self>()
        }
    }
}

impl PushNum for i16 {
    fn push_individually(self, raw: &mut Vec<u8>) -> usize {
        if self >= 0 && self < 0x100 {
            push_u8(self as u8, raw)
        } else if self < 0 && self >= -0x80 {
            push_i8(self as i8, raw)
        } else {
            push_i16(self, raw)
        }
    }

    push_num_array!(write_i16, b's');
}

impl PushNum for u16 {
    fn push_individually(self, raw: &mut Vec<u8>) -> usize {
        if self < 0x100 {
            push_u8(self as u8, raw)
        } else {
            push_u16(self, raw)
        }
    }

    push_num_array!(write_u16, b'S');
}

impl PushNum for i32 {
    fn push_individually(self, raw: &mut Vec<u8>) -> usize {
        if self >= 0 {
            if self < 0x100 {
                push_u8(self as u8, raw)
            } else if self < 0x10000 {
                push_u16(self as u16, raw)
            } else {
                push_i32(self, raw)
            }
        } else {
            if self >= -0x80 {
                push_i8(self as i8, raw)
            } else if self >= -0x8000 {
                push_i16(self as i16, raw)
            } else {
                push_i32(self, raw)
            }
        }
    }

    push_num_array!(write_i32, b'i');
}

impl PushNum for u32 {
    fn push_individually(self, raw: &mut Vec<u8>) -> usize {
        if self < 0x100 {
            push_u8(self as u8, raw)
        } else if self < 0x10000 {
            push_u16(self as u16, raw)
        } else {
            push_u32(self, raw)
        }
    }

    push_num_array!(write_u32, b'I');
}

impl PushNum for f32 {
    fn push_individually(self, raw: &mut Vec<u8>) -> usize {
        raw.push(b'f');
        raw.write_f32::<LittleEndian>(self).unwrap();
        1 + mem::size_of::<f32>()
    }

    push_num_array!(write_f32, b'f');
}

/// Wrapper around raw tags.
///
/// Allows to get and modify record tags. Method [get](#method.get) returns a
/// [TagValue](enum.TagValue.html), which is a viewer of raw data (it does not copy the data,
/// unless it has a numeric/char type).
///
/// ```rust
/// // Read tag ZZ, which can be either u32 array or a char:
/// match record.tags().get(b"ZZ") {
///     Some(TagValue::IntArray(array_view)) => {
///         assert!(array_view.int_type() == IntegerType::U32);
///         println!("Array[2] = {}", array_view.at(2));
///     },
///     Some(TagValue::Char(value)) => println!("Char = {}", value),
///     _ => panic!("Unexpected type"),
/// }
/// ```
///
/// Methods [push](#method.push) and [insert](#method.insert) add a tag (or modify a tag),
/// and they provide a convinient way to specify tag type. For numeric and string types, you
/// can just write
/// ```rust
/// // Add a new tag with name `AA`, type `i32` and value `10`.
/// record.tags_mut().push(b"AA", 10_i32);
/// // Add a new tag with name `ZZ`, type `string` and value `"abcd"`.
/// record.tags_mut().push(b"ZZ", "abcd");
/// ```
/// To add a Hex value you need to wrap `&[u8]` in a [Hex](struct.Hex.html) wrapper:
/// ```rust
/// // Add a new tag with name `HH`, type `hex` and value `"FF00"`.
/// record.tags_mut().push(b"HH", Hex(b"FF00"));
/// ```
/// Finally, to specify a numeric array, you may need to coerce the array to a splice:
/// ```rust
/// // Add a new tag with name `BB`, type `i16 array` and value `[3, 4, 5, 6]`.
/// record.tags_mut().push(b"BB", &[3_i16, 4, 5, 6] as &[i16]);
/// ```
#[derive(Clone)]
pub struct TagViewer {
    raw: Vec<u8>,
    lengths: Vec<u32>,
}

/// Get a size of type from letter (c -> 1), (i -> 4). Returns Error for non int/float types
fn tag_type_size(ty: u8) -> io::Result<u32> {
    match ty {
        b'c' | b'C' | b'A' => Ok(1),
        b's' | b'S' => Ok(2),
        b'i' | b'I' | b'f' => Ok(4),
        _ => Err(io::Error::new(InvalidData, format!("Corrupted record: Unexpected tag type: {}", ty as char))),
    }
}

/// Get the length of the first tag (including name) in a raw tags array.
///
/// For example, the function would return 7 for the raw representation of `"AA:i:10    BB:i:20"`.
fn get_length(raw_tags: &[u8]) -> io::Result<u32> {
    if raw_tags.len() < 4 {
        return Err(io::Error::new(InvalidData, "Corrupted record: Truncated tags"));
    }
    let ty = raw_tags[2];
    match ty {
        b'Z' | b'H' => {
            for i in 3..raw_tags.len() {
                if raw_tags[i] == 0 {
                    if ty == b'H' && i % 2 != 0 {
                        return Err(io::Error::new(InvalidData, "Corrupted record: Hex tag has an odd number of bytes"));
                    }
                    return Ok(1 + i as u32);
                }
            }
            Err(io::Error::new(InvalidData, "Corrupted record: Truncated tags"))
        },
        b'B' => {
            if raw_tags.len() < 8 {
                return Err(io::Error::new(InvalidData, "Corrupted record: Truncated tags"));
            }
            let arr_len = (&raw_tags[4..8]).read_i32::<LittleEndian>()? as u32;
            Ok(8 + tag_type_size(raw_tags[3])? * arr_len)
        },
        _ => Ok(3 + tag_type_size(raw_tags[2])?),
    }
}

/// Alias for a tag name.
pub type TagName = [u8; 2];

impl TagViewer {
    /// Create a new tag viewer
    pub(crate) fn new() -> Self {
        Self {
            raw: Vec::new(),
            lengths: Vec::new(),
        }
    }

    pub(crate) fn shrink_to_fit(&mut self) {
        self.raw.shrink_to_fit();
        self.lengths.shrink_to_fit();
    }

    pub(crate) fn fill_from<R: Read>(&mut self, stream: &mut R, length: usize) -> io::Result<()> {
        unsafe {
            resize(&mut self.raw, length);
        }
        stream.read_exact(&mut self.raw)?;
        self.lengths.clear();
        let mut sum_len = 0;
        while sum_len < self.raw.len() {
            let tag_len = get_length(&self.raw[sum_len..])?;
            self.lengths.push(tag_len);
            sum_len += tag_len as usize;
        }
        if sum_len > self.raw.len() {
            Err(io::Error::new(InvalidData, "Corrupted record: Truncated tags"))
        } else {
            Ok(())
        }
    }

    /// Returns raw data.
    pub fn raw(&self) -> &[u8] {
        &self.raw
    }

    /// Clears the contents but does not touch capacity.
    pub fn clear(&mut self) {
        self.lengths.clear();
        self.raw.clear();
    }

    /// Returns a value of a tag with `name`. Value integer/float types, returns copied value, for
    /// array and string types returns a wrapper over reference. Takes `O(n_tags)`.
    pub fn get<'a>(&'a self, name: &TagName) -> Option<TagValue<'a>> {
        let mut start = 0;
        for &tag_len in self.lengths.iter() {
            let tag_len = tag_len as usize;
            if name == &self.raw[start..start + 2] {
                return Some(TagValue::from_raw(self.raw[start + 2],
                    &self.raw[start + 3..start + tag_len]));
            }
            start += tag_len;
        }
        None
    }

    /// Iterate over tuples `(name, tag_value)`, where `name: [u8; 2]` and `tag_value: TagValue`.
    pub fn iter<'a>(&'a self) -> TagIter<'a> {
        TagIter {
            pos: 0,
            lengths: self.lengths.iter(),
            raw: &self.raw,
        }
    }

    /// Appends a new tag. Trait [WriteValue](trait.WriteValue.html) is implemented for possible
    /// tag values, so you can add new tags like this:
    /// ```rust
    /// record.tags_mut().push(b"AA", 10)
    /// ```
    /// Due to Rust constraints, you may need to explicitly coerce a numeric array to a slice,
    /// for example
    /// ```rust
    /// record.tags_mut().push(b"BB", &[10_i16, 20, 30] as &[i16]).unwrap();
    /// ```
    ///
    /// This function does not check if there is already a tag with the same name.
    /// Takes `O(new_tag_len)`.
    ///
    /// See [WriteValue](trait.WriteValue.html) for more information. 
    pub fn push<V: WriteValue>(&mut self, name: &TagName, value: V) -> io::Result<()> {
        self.raw.push(name[0]);
        self.raw.push(name[1]);
        self.lengths.push(2 + value.write(&mut self.raw)? as u32);
        Ok(())
    }

    /// Inserts a new tag instead of existing.
    /// If there is no tags with the same name, pushes to
    /// the end and returns `None`. Takes `O(raw_tags_len + new_tag_len)`.
    pub fn insert<'a, V: WriteValue>(&'a mut self, name: &TagName, value: V) -> io::Result<()> {
        let mut start = 0;
        for i in 0..self.lengths.len() {
            let tag_len = self.lengths[i] as usize;
            if name == &self.raw[start..start + 2] {
                let mut new_tag = Vec::new();
                let new_len = value.write(&mut new_tag)? as u32;
                self.raw.splice(start + 2..start + tag_len, new_tag);
                self.lengths[i] = 2 + new_len;
                return Ok(());
            }
            start += tag_len;
        }
        self.push(name, value)?;
        Ok(())
    }

    /// Removes a tag if present. Returns `true` if the tag existed and `false` otherwise.
    /// Takes `O(raw_tags_len)`.
    pub fn remove(&mut self, name: &TagName) -> bool {
        let mut start = 0;
        for i in 0..self.lengths.len() {
            let tag_len = self.lengths[i] as usize;
            if name == &self.raw[start..start + 2] {
                self.raw.drain(start..start + tag_len);
                self.lengths.remove(i);
                return true;
            }
            start += tag_len;
        }
        false
    }

    /// Writes tags in a SAM format.
    pub fn write_sam<W: Write>(&self, f: &mut W) -> io::Result<()> {
        for (name, value) in self.iter() {
            f.write_all(&[b'\t', name[0], name[1], b':'])?;
            value.write_sam(f)?;
        }
        Ok(())
    }

    /// Pushes integer in a smallest possible format.
    pub fn push_int(&mut self, name: &TagName, value: i64) {
        if value >= 0 {
            if value < 0x100 {
                self.push(name, value as u8).expect("Failed to push int tag");
            } else if value < 0x10000 {
                self.push(name, value as u16).expect("Failed to push int tag");
            } else {
                self.push(name, value as u32).expect("Failed to push int tag");
            }
        } else {
            if value >= -0x80 {
                self.push(name, value as i8).expect("Failed to push int tag");
            } else if value >= -0x8000 {
                self.push(name, value as i16).expect("Failed to push int tag");
            } else {
                self.push(name, value as i32).expect("Failed to push int tag");
            }
        }
    }

    /// Returns false if failed to parse.
    fn push_sam_array(&mut self, name: &TagName, value: &str) -> bool {
        let mut split = value.split(',');
        let arr_type = match split.next() {
            Some(v) => v,
            None => return false,
        };

        match arr_type {
            "c" => split.map(|s| s.parse::<i8>()).collect::<Result<Vec<_>, _>>()
                .map(|values| self.push(name, &values as &[i8])).is_ok(),
            "C" => split.map(|s| s.parse::<u8>()).collect::<Result<Vec<_>, _>>()
                .map(|values| self.push(name, &values as &[u8])).is_ok(),
            "s" => split.map(|s| s.parse::<i16>()).collect::<Result<Vec<_>, _>>()
                .map(|values| self.push(name, &values as &[i16])).is_ok(),
            "S" => split.map(|s| s.parse::<u16>()).collect::<Result<Vec<_>, _>>()
                .map(|values| self.push(name, &values as &[u16])).is_ok(),
            "i" => split.map(|s| s.parse::<i32>()).collect::<Result<Vec<_>, _>>()
                .map(|values| self.push(name, &values as &[i32])).is_ok(),
            "I" => split.map(|s| s.parse::<u32>()).collect::<Result<Vec<_>, _>>()
                .map(|values| self.push(name, &values as &[u32])).is_ok(),
            "f" => split.map(|s| s.parse::<f32>()).collect::<Result<Vec<_>, _>>()
                .map(|values| self.push(name, &values as &[f32])).is_ok(),
            _ => false,
        }
    }

    /// Returns false if failed to parse.
    fn inner_push_sam(&mut self, tag: &str) -> bool {
        let tag_bytes = tag.as_bytes();
        // 012345...
        // nn:t:value
        if tag_bytes.len() < 5 || tag_bytes[2] != b':' || tag_bytes[4] != b':' {
            return true;
        }
        let tag_name = &[tag_bytes[0], tag_bytes[1]];
        let tag_type = tag_bytes[3];
        let tag_value = unsafe { std::str::from_utf8_unchecked(&tag_bytes[5..]) };

        match tag_type {
            b'A' => {
                if tag_bytes.len() != 6 {
                    return false;
                }
                self.push(tag_name, tag_bytes[5] as char).is_ok()
            },
            b'i' => {
                let number: i64 = match tag_value.parse() {
                    Ok(value) => value,
                    Err(_) => return false,
                };
                self.push_int(tag_name, number);
                true
            },
            b'f' => {
                let float: f32 = match tag_value.parse() {
                    Ok(value) => value,
                    Err(_) => return false,
                };
                self.push(tag_name, float).is_ok()
            },
            b'Z' => {
                self.push(tag_name, tag_value).is_ok()
            },
            b'H' => {
                self.push(tag_name, Hex(tag_value.as_bytes())).is_ok()
            },
            b'B' => {
                self.push_sam_array(tag_name, tag_value)
            },
            _ => false,
        }
    }

    /// Adds a new tag in SAM format (name:type:value).
    pub fn push_sam(&mut self, tag: &str) -> io::Result<()> {
        if self.inner_push_sam(tag) {
            Ok(())
        } else {
            Err(io::Error::new(InvalidData, format!("Cannot parse tag '{}'", tag)))
        }
    }
}

/// Iterator over tags.
#[derive(Clone)]
pub struct TagIter<'a> {
    pos: usize,
    lengths: std::slice::Iter<'a, u32>,
    raw: &'a [u8],
}

impl<'a> Iterator for TagIter<'a> {
    type Item = (TagName, TagValue<'a>);

    fn next(&mut self) -> Option<Self::Item> {
        match self.lengths.next() {
            Some(tag_len) => {
                let tag_len = *tag_len as usize;
                let start = self.pos;
                self.pos += tag_len;
                let name = [self.raw[start], self.raw[start + 1]];
                Some((name, TagValue::from_raw(self.raw[start + 2],
                        &self.raw[start + 3..start + tag_len])))
            },
            None => None,
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.lengths.size_hint()
    }
}

impl<'a> ExactSizeIterator for TagIter<'a> {}

impl<'a> std::iter::FusedIterator for TagIter<'a> {}
