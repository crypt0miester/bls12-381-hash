//! Witness-assisted RFC 9380 hash-to-curve for BLS12-381, for Solana SBF.
//!
//! Suites are feature gated. The DST is a runtime parameter. Host-side witness
//! generation ships in the same crate under `witness`, gated on the target.

#![no_std]
#![allow(dead_code)]
#![allow(unexpected_cfgs)]

extern crate alloc;

mod consts_g1;
mod consts_g2;
mod fp;
mod macros;
mod fp2;
mod g1;
mod g2;

use solana_program_error::ProgramError;

/// Verification failure: a witness did not satisfy its check, or an input was
/// malformed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Error {
    BadWitness,
    BadInput,
}

impl From<Error> for ProgramError {
    fn from(_: Error) -> Self {
        ProgramError::InvalidInstructionData
    }
}

/// Standard CFRG / RFC 9380 domain separation tags for the RO and NU suites.
pub mod dst {
    pub const G1_RO: &[u8] = b"BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_POP_";
    pub const G2_RO: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_";
    pub const G1_NU: &[u8] = b"BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_NU_POP_";
    pub const G2_NU: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_NU_POP_";
}

#[cfg(feature = "g1-ro")]
pub use crate::g1::hash_to_g1;
#[cfg(feature = "g2-ro")]
pub use crate::g2::hash_to_g2;
#[cfg(feature = "g2-ro")]
pub use crate::g2::hash_to_g2_compact;
#[cfg(feature = "g2-ro")]
pub use crate::g2::hash_to_g2_compact_xgcd;
#[doc(hidden)]
#[cfg(feature = "g1-ro")]
pub use crate::g1::hash_to_g1_prefix;
#[doc(hidden)]
#[cfg(feature = "g2-ro")]
pub use crate::g2::hash_to_g2_prefix;
#[cfg(feature = "g1-nu")]
pub use crate::g1::encode_to_g1;
#[cfg(feature = "g2-nu")]
pub use crate::g2::encode_to_g2;
#[cfg(feature = "modexp")]
pub use crate::g1::run as hash_to_g1_modexp;

/// Per function CU probes for the bench harness; not part of the API
#[doc(hidden)]
pub mod probe;

/// Host-side witness generation, one module per suite.
#[cfg(not(target_os = "solana"))]
pub mod witness {
    pub use crate::g1::witness as g1;
    pub use crate::g2::witness as g2;
}
