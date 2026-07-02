// AVX-512 k-mer marker extraction.
//
// This is a wider-SIMD counterpart to `avx2_seeding.rs`. The AVX2 path packs
// four k-mers into a `__m256i` (four 64-bit lanes); AVX-512 doubles that to
// eight 64-bit lanes in a `__m512i`, so twice as many k-mers are rolled,
// hashed and thresholded per iteration.
//
// The idea (process many overlapping k-mers in parallel with the widest vector
// the CPU offers, and compact the survivors of the sampling threshold to the
// front of the output with a single native instruction) is taken from Igor
// Martayan's PhD thesis on SIMD k-mer algorithms
// (https://github.com/imartayan/phd, "simd-minimizers" chapter). There the
// compaction is a shuffle + lookup table on AVX2; on AVX-512 it is the single
// `vpcompressq` instruction, exposed here as `_mm512_mask_compressstoreu_epi64`.
//
// Output invariant: this path emits exactly the same *multiset* of hashes as
// `avx2_seeding::extract_markers_avx2` for the same input (only the order in
// which they land in `kmer_vec` differs, which does not matter — sample
// sketches are hash multisets and genome sketches are sorted downstream). In
// particular it reproduces the AVX2 path's trailing-window behaviour: both
// process exactly the first `4 * ((n - k + 1) / 4)` windows of the sequence,
// so switching between the two never changes a sketch's contents. This keeps
// weebill/sylph sketch compatibility intact regardless of which CPU built the
// sketch.

use crate::types::*;
use std::arch::x86_64::*;

/// AVX-512 form of `mm_hash64` / `mm_hash256`, applied to eight 64-bit lanes.
/// Bit-for-bit identical to the scalar `mm_hash64`, including the historical
/// first-step quirk (`!(key + (key << 21))`).
#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn mm_hash512(kmer: __m512i) -> __m512i {
    let ones = _mm512_set1_epi64(-1);
    let mut key = kmer;
    let s1 = _mm512_slli_epi64(key, 21);
    key = _mm512_add_epi64(key, s1);
    // NOT: xor with all-ones (the AVX2 path uses cmpeq(x,x); AVX-512 compares
    // produce a mask register, so use an explicit all-ones vector instead).
    key = _mm512_xor_si512(key, ones);

    key = _mm512_xor_si512(key, _mm512_srli_epi64(key, 24));
    let s2 = _mm512_slli_epi64(key, 3);
    let s3 = _mm512_slli_epi64(key, 8);
    key = _mm512_add_epi64(key, s2);
    key = _mm512_add_epi64(key, s3);
    key = _mm512_xor_si512(key, _mm512_srli_epi64(key, 14));
    let s4 = _mm512_slli_epi64(key, 2);
    let s5 = _mm512_slli_epi64(key, 4);
    key = _mm512_add_epi64(key, s4);
    key = _mm512_add_epi64(key, s5);
    key = _mm512_xor_si512(key, _mm512_srli_epi64(key, 28));
    let s6 = _mm512_slli_epi64(key, 31);
    key = _mm512_add_epi64(key, s6);
    key
}

/// Number of trailing windows shared with the AVX2 path: both process exactly
/// the first `4 * ((n - k + 1) / 4)` windows and drop the (< 4) remainder.
#[inline]
fn avx2_compatible_window_count(len: usize, k: usize) -> usize {
    if len < k {
        return 0;
    }
    let nwin = len - k + 1;
    nwin - (nwin % 4)
}

#[target_feature(enable = "avx512f")]
pub unsafe fn extract_markers_avx512(string: &[u8], kmer_vec: &mut Vec<u64>, c: usize, k: usize) {
    let total_windows = avx2_compatible_window_count(string.len(), k);
    if total_windows == 0 {
        return;
    }

    // Eight contiguous lanes over the window range [0, total_windows). Windows
    // [8*lane_len, total_windows) are handled by a scalar remainder loop so the
    // union of processed windows is exactly [0, total_windows) — matching AVX2.
    let lane_len = total_windows / 8;
    let simd_windows = lane_len * 8;

    if lane_len > 0 {
        let two_k_minus_2 = 2 * (k - 1) as i32;
        let marker_mask = (Kmer::MAX >> (std::mem::size_of::<Kmer>() * 8 - 2 * k)) as i64;
        let rev_marker_mask: i64 = !(3i64 << (2 * k - 2));
        let threshold_marker = u64::MAX / c as u64;

        // Per-lane starting byte offset (lane j covers windows [j*lane_len, ..)).
        let off: [usize; 8] = std::array::from_fn(|j| j * lane_len);

        let mm_marker_mask = _mm512_set1_epi64(marker_mask);
        let mm_rev_marker_mask = _mm512_set1_epi64(rev_marker_mask);
        let mm_threshold = _mm512_set1_epi64(threshold_marker as i64);
        let rev_sub = _mm512_set1_epi64(3);

        let mut rolling_f = _mm512_setzero_si512();
        let mut rolling_r = _mm512_setzero_si512();

        // Prime the first k-1 bases of every lane.
        for i in 0..k - 1 {
            let f_nucs = _mm512_set_epi64(
                BYTE_TO_SEQ[string[off[7] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[6] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[5] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[4] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[3] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[2] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[1] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[0] + i] as usize] as i64,
            );
            let r_nucs = _mm512_sub_epi64(rev_sub, f_nucs);
            rolling_f = _mm512_slli_epi64(rolling_f, 2);
            rolling_f = _mm512_or_si512(rolling_f, f_nucs);
            rolling_r = _mm512_srli_epi64(rolling_r, 2);
            rolling_r = _mm512_or_si512(
                rolling_r,
                _mm512_sllv_epi64(r_nucs, two_k_minus_2_vec(two_k_minus_2)),
            );
        }

        for i in k - 1..lane_len + k - 1 {
            let f_nucs = _mm512_set_epi64(
                BYTE_TO_SEQ[string[off[7] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[6] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[5] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[4] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[3] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[2] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[1] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[0] + i] as usize] as i64,
            );
            let r_nucs = _mm512_sub_epi64(rev_sub, f_nucs);

            rolling_f = _mm512_slli_epi64(rolling_f, 2);
            rolling_f = _mm512_or_si512(rolling_f, f_nucs);
            rolling_f = _mm512_and_si512(rolling_f, mm_marker_mask);

            rolling_r = _mm512_srli_epi64(rolling_r, 2);
            rolling_r = _mm512_and_si512(rolling_r, mm_rev_marker_mask);
            rolling_r = _mm512_or_si512(
                rolling_r,
                _mm512_sllv_epi64(r_nucs, two_k_minus_2_vec(two_k_minus_2)),
            );

            // canonical = min(forward, reverse) as unsigned (k <= 31 so the top
            // bit is always 0; this matches the AVX2 signed-cmp+blend exactly).
            let canonical = _mm512_min_epu64(rolling_f, rolling_r);
            let hash = mm_hash512(canonical);

            // Branchless survivor compaction (vpcompressq): write only lanes
            // whose hash is below the sampling threshold, contiguously.
            let mask = _mm512_cmplt_epu64_mask(hash, mm_threshold);
            if mask != 0 {
                if kmer_vec.capacity() - kmer_vec.len() < 8 {
                    kmer_vec.reserve(8);
                }
                let ptr = kmer_vec.as_mut_ptr().add(kmer_vec.len());
                _mm512_mask_compressstoreu_epi64(ptr as *mut i64, mask, hash);
                kmer_vec.set_len(kmer_vec.len() + mask.count_ones() as usize);
            }
        }
    }

    // Scalar remainder: windows [simd_windows, total_windows).
    if simd_windows < total_windows {
        scalar_window_range(string, kmer_vec, c, k, simd_windows, total_windows);
    }
}

/// Broadcast a shift count for `_mm512_sllv_epi64` (per-lane variable shift).
#[inline]
#[target_feature(enable = "avx512f")]
unsafe fn two_k_minus_2_vec(shift: i32) -> __m512i {
    _mm512_set1_epi64(shift as i64)
}

/// Scalar fallback that processes a half-open range of *window indices*
/// `[start_window, end_window)` of `string`, matching `seeding::fmh_seeds`.
#[inline]
fn scalar_window_range(
    string: &[u8],
    kmer_vec: &mut Vec<u64>,
    c: usize,
    k: usize,
    start_window: usize,
    end_window: usize,
) {
    let marker_mask = u64::MAX >> (64 - 2 * k);
    let marker_rev_mask = !(3u64 << (2 * k - 2));
    let reverse_shift_dist = 2 * (k - 1);
    let threshold_marker = u64::MAX / (c as u64);

    let mut f: u64 = 0;
    let mut r: u64 = 0;
    // Prime k-1 bases ending just before `start_window`'s k-mer start.
    let base = start_window;
    for i in 0..k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[base + i] as usize] as u64;
        let nuc_r = 3 - nuc_f;
        f <<= 2;
        f |= nuc_f;
        r >>= 2;
        r |= nuc_r << reverse_shift_dist;
    }
    for w in start_window..end_window {
        let i = w + k - 1;
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        let nuc_r = 3 - nuc_f;
        f <<= 2;
        f |= nuc_f;
        f &= marker_mask;
        r >>= 2;
        r &= marker_rev_mask;
        r |= nuc_r << reverse_shift_dist;
        let canonical = if f < r { f } else { r };
        let hash = crate::seeding::mm_hash64(canonical);
        if hash < threshold_marker {
            kmer_vec.push(hash);
        }
    }
}

#[target_feature(enable = "avx512f")]
pub unsafe fn extract_markers_avx512_positions(
    string: &[u8],
    kmer_vec: &mut Vec<(usize, usize, u64)>,
    c: usize,
    k: usize,
    contig_number: usize,
) {
    let total_windows = avx2_compatible_window_count(string.len(), k);
    if total_windows == 0 {
        return;
    }
    let lane_len = total_windows / 8;
    let simd_windows = lane_len * 8;

    if lane_len > 0 {
        let two_k_minus_2 = 2 * (k - 1) as i32;
        let marker_mask = (Kmer::MAX >> (std::mem::size_of::<Kmer>() * 8 - 2 * k)) as i64;
        let rev_marker_mask: i64 = !(3i64 << (2 * k - 2));
        let threshold_marker = u64::MAX / c as u64;

        let off: [usize; 8] = std::array::from_fn(|j| j * lane_len);
        let mm_marker_mask = _mm512_set1_epi64(marker_mask);
        let mm_rev_marker_mask = _mm512_set1_epi64(rev_marker_mask);
        let rev_sub = _mm512_set1_epi64(3);

        let mut rolling_f = _mm512_setzero_si512();
        let mut rolling_r = _mm512_setzero_si512();

        for i in 0..k - 1 {
            let f_nucs = _mm512_set_epi64(
                BYTE_TO_SEQ[string[off[7] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[6] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[5] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[4] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[3] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[2] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[1] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[0] + i] as usize] as i64,
            );
            let r_nucs = _mm512_sub_epi64(rev_sub, f_nucs);
            rolling_f = _mm512_slli_epi64(rolling_f, 2);
            rolling_f = _mm512_or_si512(rolling_f, f_nucs);
            rolling_r = _mm512_srli_epi64(rolling_r, 2);
            rolling_r = _mm512_or_si512(
                rolling_r,
                _mm512_sllv_epi64(r_nucs, two_k_minus_2_vec(two_k_minus_2)),
            );
        }

        let mut hashes = [0u64; 8];
        for i in k - 1..lane_len + k - 1 {
            let f_nucs = _mm512_set_epi64(
                BYTE_TO_SEQ[string[off[7] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[6] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[5] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[4] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[3] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[2] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[1] + i] as usize] as i64,
                BYTE_TO_SEQ[string[off[0] + i] as usize] as i64,
            );
            let r_nucs = _mm512_sub_epi64(rev_sub, f_nucs);

            rolling_f = _mm512_slli_epi64(rolling_f, 2);
            rolling_f = _mm512_or_si512(rolling_f, f_nucs);
            rolling_f = _mm512_and_si512(rolling_f, mm_marker_mask);

            rolling_r = _mm512_srli_epi64(rolling_r, 2);
            rolling_r = _mm512_and_si512(rolling_r, mm_rev_marker_mask);
            rolling_r = _mm512_or_si512(
                rolling_r,
                _mm512_sllv_epi64(r_nucs, two_k_minus_2_vec(two_k_minus_2)),
            );

            let canonical = _mm512_min_epu64(rolling_f, rolling_r);
            let hash = mm_hash512(canonical);
            _mm512_storeu_si512(hashes.as_mut_ptr() as *mut __m512i, hash);

            let w = i - (k - 1);
            for (j, &h) in hashes.iter().enumerate() {
                if h < threshold_marker {
                    kmer_vec.push((contig_number, off[j] + w, h));
                }
            }
        }
    }

    if simd_windows < total_windows {
        for w in simd_windows..total_windows {
            let mut f: u64 = 0;
            let mut r: u64 = 0;
            let marker_mask = u64::MAX >> (64 - 2 * k);
            let marker_rev_mask = !(3u64 << (2 * k - 2));
            let reverse_shift_dist = 2 * (k - 1);
            let threshold_marker = u64::MAX / (c as u64);
            for x in 0..k {
                let nuc_f = BYTE_TO_SEQ[string[w + x] as usize] as u64;
                let nuc_r = 3 - nuc_f;
                f <<= 2;
                f |= nuc_f;
                f &= marker_mask;
                r >>= 2;
                r &= marker_rev_mask;
                r |= nuc_r << reverse_shift_dist;
            }
            let canonical = if f < r { f } else { r };
            let hash = crate::seeding::mm_hash64(canonical);
            if hash < threshold_marker {
                kmer_vec.push((contig_number, w, hash));
            }
        }
    }
}
