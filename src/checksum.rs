//! Integrity checksums for the seekable databases (`.sylref`, `.syl2db`).
//!
//! The single-frame formats (`.sylspc`, `.syldbc`, `.sylspr`) get corruption
//! detection for free from zstd's trailing content checksum: reading one means
//! decoding the whole frame, so the checksum is always verified. The two seekable
//! databases are different — a query reads a header, a footer and a handful of
//! blocks, never the whole file, so nothing validates an end-to-end digest as a
//! side effect. They instead store an XXH64 of everything after their fixed-size
//! header, checked on demand ([`crate::refdelta::RefIndex::verify_checksum`],
//! [`crate::twostage_db::TwoStageDb::verify_checksum`]); `weebill inspect` does so
//! for every database it is given.
//!
//! The header itself is not covered: it is patched last (the offsets are only
//! known once the body has been streamed out), and every field in it — magic,
//! version, section offsets — is validated when the file is opened.

use std::hash::Hasher;
use std::io::{self, Read, Write};
use twox_hash::XxHash64;

/// XXH64 of every byte `r` yields.
pub(crate) fn hash_reader<R: Read>(mut r: R) -> io::Result<u64> {
    let mut h = XxHash64::with_seed(0);
    let mut buf = vec![0u8; 1 << 20];
    loop {
        match r.read(&mut buf)? {
            0 => return Ok(h.finish()),
            n => h.write(&buf[..n]),
        }
    }
}

/// The error raised when a database's contents do not hash to the checksum stored
/// in its header.
pub(crate) fn mismatch(what: &str, expected: u64, got: u64) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "{} is corrupt: its contents hash to {:016x}, but {:016x} was recorded when it was written",
            what, got, expected
        ),
    )
}

/// A writer that XXH64s everything written through it, so a database body that is
/// streamed straight to disk (and never held in RAM) can still be checksummed.
pub(crate) struct HashingWriter<W> {
    inner: W,
    hasher: XxHash64,
}

impl<W: Write> HashingWriter<W> {
    pub(crate) fn new(inner: W) -> Self {
        HashingWriter {
            inner,
            hasher: XxHash64::with_seed(0),
        }
    }

    /// The wrapped writer, plus the digest of everything written through it.
    pub(crate) fn finish(self) -> (W, u64) {
        let digest = self.hasher.finish();
        (self.inner, digest)
    }
}

impl<W: Write> Write for HashingWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let n = self.inner.write(buf)?;
        self.hasher.write(&buf[..n]);
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}
