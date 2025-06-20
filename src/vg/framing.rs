//! Implementation of vg's Protobuf framing format
//!
//! This format is used to store multiple Protobuf messages in a single file.
//! Each message is preceded by its length as a varint64, and messages are
//! grouped with a count of messages per group. The first message in a group
//! can optionally be a type tag string.

use anyhow::{anyhow, Result};
use protobuf::Message;
use std::io::{Read, BufRead};

/// Reads a varint64 from a reader
fn read_varint64<R: Read>(reader: &mut R) -> Result<u64> {
    let mut result: u64 = 0;
    let mut shift: u32 = 0;

    loop {
        let mut buf = [0u8];
        reader.read_exact(&mut buf)?;
        let byte = buf[0];

        result |= ((byte & 0x7f) as u64) << shift;
        if (byte & 0x80) == 0 {
            break;
        }
        shift += 7;
        if shift >= 64 {
            return Err(anyhow!("Malformed varint64"));
        }
    }

    Ok(result)
}

/// Represents a group in a VG file
#[derive(Debug)]
pub struct Group {
    /// Number of messages in the group
    pub count: u64,
    /// Optional type tag (usually present as first message)
    pub type_tag: Option<String>,
    /// Raw message data
    pub messages: Vec<Vec<u8>>,
}

/// Reads a single group from a reader
pub fn read_group<R: BufRead>(reader: &mut R) -> Result<Option<Group>> {
    // Try to read the count, return None if we're at EOF
    let count = match read_varint64(reader) {
        Ok(count) => count,
        Err(e) if e.to_string().contains("UnexpectedEof") => return Ok(None),
        Err(e) => return Err(e),
    };

    let mut messages = Vec::new();
    let mut type_tag = None;

    // Read the first message (potential type tag)
    if count > 0 {
        let size = read_varint64(reader)?;
        let mut buf = vec![0u8; size as usize];
        reader.read_exact(&mut buf)?;

        // Try to interpret as type tag
        if let Ok(tag) = String::from_utf8(buf.clone()) {
            if tag.chars().all(|c| c.is_ascii_alphanumeric()) {
                type_tag = Some(tag);
            } else {
                messages.push(buf);
            }
        } else {
            messages.push(buf);
        }

        // Read remaining messages
        for _ in 1..count {
            let size = read_varint64(reader)?;
            let mut buf = vec![0u8; size as usize];
            reader.read_exact(&mut buf)?;
            messages.push(buf);
        }
    }

    Ok(Some(Group {
        count,
        type_tag,
        messages,
    }))
}

/// Iterator over groups in a VG file
pub struct GroupIterator<R: BufRead> {
    reader: R,
}

impl<R: BufRead> GroupIterator<R> {
    /// Creates a new iterator over groups in a VG file
    pub fn new(reader: R) -> Self {
        Self { reader }
    }
}

impl<R: BufRead> Iterator for GroupIterator<R> {
    type Item = Result<Group>;

    fn next(&mut self) -> Option<Self::Item> {
        match read_group(&mut self.reader) {
            Ok(Some(group)) => Some(Ok(group)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_read_varint64() {
        // Test cases with known values
        let cases = vec![
            (vec![0x00], 0u64),
            (vec![0x01], 1u64),
            (vec![0x7f], 127u64),
            (vec![0x80, 0x01], 128u64),
            (vec![0xff, 0x01], 255u64),
        ];

        for (bytes, expected) in cases {
            let mut reader = Cursor::new(bytes);
            let result = read_varint64(&mut reader).unwrap();
            assert_eq!(result, expected);
        }
    }

    #[test]
    fn test_group_iterator() {
        // Create a test file with two groups
        let mut data = Vec::new();

        // Group 1: GAM type tag + 1 message
        data.extend(&[0x02]); // count = 2
        data.extend(&[0x03]); // size = 3
        data.extend(b"GAM"); // type tag
        data.extend(&[0x02]); // size = 2
        data.extend(&[0x01, 0x02]); // message data

        // Group 2: 1 message without type tag
        data.extend(&[0x01]); // count = 1
        data.extend(&[0x03]); // size = 3
        data.extend(&[0x03, 0x04, 0x05]); // message data

        let mut reader = Cursor::new(data);
        let groups: Vec<Group> = GroupIterator::new(&mut reader)
            .collect::<Result<Vec<_>>>()
            .unwrap();

        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].type_tag, Some("GAM".to_string()));
        assert_eq!(groups[0].messages.len(), 1);
        assert_eq!(groups[1].type_tag, None);
        assert_eq!(groups[1].messages.len(), 1);
    }
}