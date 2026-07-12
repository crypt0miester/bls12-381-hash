//! Pins the subgroup-check contract of the bls12-381 syscalls, which the
//! pipeline design leans on: cofactor clearing feeds pre-subgroup points
//! through the add syscall (must be accepted), the final g2_validate is a
//! real subgroup check (not just on-curve), and the pairing syscall
//! re-validates its inputs (which is what makes the final validate pure
//! defense-in-depth for a pairing-bound consumer like the min-pk verify).

use bls12_381::{G1Affine, G2Affine};
use mollusk_svm::{program::loader_keys::LOADER_V3, Mollusk};
use solana_instruction::Instruction;
use solana_pubkey::Pubkey;

const ID: Pubkey = Pubkey::new_from_array([7u8; 32]);
const MESSAGE: &[u8] = b"tapedrive vote payload: epoch 42, slot 1337, snapshot root cafebabe";

#[test]
fn syscall_subgroup_contract() {
    let elf = std::fs::read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../target/deploy/bls381_bench.so"
    ))
    .unwrap();
    let mut mollusk = Mollusk::default();
    mollusk.add_program_with_loader_and_elf(&ID, &LOADER_V3, &elf);
    mollusk.compute_budget.compute_unit_limit = 500_000_000;

    // stage 2 of the witnessed pipeline returns the pre-clearing iso output:
    // on the curve, outside the subgroup
    let witness = bls381_hash::witness::g2::generate(MESSAGE);
    let mut payload = vec![47u8, 2u8];
    payload.extend_from_slice(&witness);
    payload.extend_from_slice(MESSAGE);
    let r = mollusk.process_instruction(&Instruction::new_with_bytes(ID, &payload, vec![]), &[]);
    assert!(!r.program_result.is_err());
    let uncleared: [u8; 192] = r.return_data.as_slice().try_into().unwrap();

    let point = G2Affine::from_uncompressed_unchecked(&uncleared).unwrap();
    assert!(bool::from(point.is_on_curve()));
    assert!(!bool::from(point.is_torsion_free()), "probe point must sit outside the subgroup");

    // validate must reject it: the pipeline-final check is a subgroup check
    let mut data = vec![10u8];
    data.extend_from_slice(&uncleared);
    let v = mollusk.process_instruction(&Instruction::new_with_bytes(ID, &data, vec![]), &[]);
    assert_ne!(v.return_data[0], 0, "g2 validate accepted a non-subgroup point");

    // the pairing must reject it as an input: a pairing-bound consumer gets
    // the subgroup check again for free
    let mut data = vec![14u8];
    data.extend_from_slice(&G1Affine::generator().to_uncompressed());
    data.extend_from_slice(&uncleared);
    let p = mollusk.process_instruction(&Instruction::new_with_bytes(ID, &data, vec![]), &[]);
    assert_ne!(p.return_data[0], 0, "pairing accepted a non-subgroup g2 input");

    // the add group op must accept it: cofactor clearing runs pre-subgroup
    // points through this syscall
    let mut data = vec![11u8];
    data.extend_from_slice(&uncleared);
    data.extend_from_slice(&uncleared);
    let a = mollusk.process_instruction(&Instruction::new_with_bytes(ID, &data, vec![]), &[]);
    assert_eq!(a.return_data[0], 0, "g2 add rejected a pre-clearing point");
}
