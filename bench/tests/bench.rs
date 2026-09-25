use bls12_381::{
    hash_to_curve::{ExpandMsgXmd, HashToCurve},
    G1Affine, G1Projective, G2Affine, G2Projective, Scalar,
};
use mollusk_svm::{program::loader_keys::LOADER_V3, result::InstructionResult, Mollusk};
use solana_instruction::Instruction;
use solana_pubkey::Pubkey;

const ID: Pubkey = Pubkey::new_from_array([7u8; 32]);

const DST_G2: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_";
const DST_G1: &[u8] = b"BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_POP_";

const MESSAGE: &[u8] = b"tapedrive vote payload: epoch 42, slot 1337, snapshot root cafebabe";

fn mollusk() -> Mollusk {
    let elf = std::fs::read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../target/deploy/bls381_bench.so"
    ))
    .expect("build the program first: cd program && cargo build-sbf");
    let mut mollusk = Mollusk::default();
    mollusk.add_program_with_loader_and_elf(&ID, &LOADER_V3, &elf);
    mollusk.compute_budget.compute_unit_limit = 500_000_000;
    mollusk
}

fn run(mollusk: &Mollusk, tag: u8, payload: &[u8]) -> InstructionResult {
    let mut data = vec![tag];
    data.extend_from_slice(payload);
    let instruction = Instruction::new_with_bytes(ID, &data, vec![]);
    mollusk.process_instruction(&instruction, &[])
}

fn cu(mollusk: &Mollusk, tag: u8, payload: &[u8], label: &str) -> u64 {
    let result = run(mollusk, tag, payload);
    assert!(
        !result.program_result.is_err(),
        "{label} failed: {:?}",
        result.program_result
    );
    println!("{label}: {} CU", result.compute_units_consumed);
    result.compute_units_consumed
}

#[test]
#[ignore = "naive zkcrypto rows miscompile on current platform-tools, see BENCHMARKS.md"]
fn bench_hash_to_curve_pipeline() {
    let mollusk = mollusk();

    let full_g2 = {
        let result = run(&mollusk, 0, MESSAGE);
        assert!(!result.program_result.is_err(), "hash_to_g2 failed: {:?}", result.program_result);
        let expected = G2Affine::from(
            <G2Projective as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::hash_to_curve(
                MESSAGE, DST_G2,
            ),
        )
        .to_compressed();
        assert_eq!(result.return_data, expected.to_vec(), "hash_to_g2 output mismatch");
        println!("hash_to_g2 full (with compress): {} CU", result.compute_units_consumed);
        result.compute_units_consumed
    };

    let full_g1 = {
        let result = run(&mollusk, 1, MESSAGE);
        assert!(!result.program_result.is_err(), "hash_to_g1 failed: {:?}", result.program_result);
        let expected = G1Affine::from(
            <G1Projective as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::hash_to_curve(
                MESSAGE, DST_G1,
            ),
        )
        .to_compressed();
        assert_eq!(result.return_data, expected.to_vec(), "hash_to_g1 output mismatch");
        println!("hash_to_g1 full (with compress): {} CU", result.compute_units_consumed);
        result.compute_units_consumed
    };

    let field = cu(&mollusk, 2, MESSAGE, "hash_to_field only");
    let mapped = cu(&mollusk, 3, MESSAGE, "field + 2x map_to_curve + add");
    let cleared = cu(&mollusk, 4, MESSAGE, "field + map + add + clear_h");

    println!();
    println!("phase breakdown (G2):");
    println!("  hash_to_field:      {field} CU");
    println!("  2x map_to_curve:    {} CU", mapped.saturating_sub(field));
    println!("  clear_cofactor:     {} CU", cleared.saturating_sub(mapped));
    println!("  affine + compress:  {} CU", full_g2.saturating_sub(cleared));
    println!("  total G2:           {full_g2} CU");
    println!("  total G1:           {full_g1} CU");
}

#[test]
#[ignore = "naive zkcrypto rows miscompile on current platform-tools, see BENCHMARKS.md"]
fn bench_hash_to_g2_matches_blst() {
    let mollusk = mollusk();
    let result = run(&mollusk, 0, MESSAGE);
    assert!(!result.program_result.is_err());

    let mut point = blst::blst_p2::default();
    let mut compressed = [0u8; 96];
    unsafe {
        blst::blst_hash_to_g2(
            &mut point,
            MESSAGE.as_ptr(),
            MESSAGE.len(),
            DST_G2.as_ptr(),
            DST_G2.len(),
            std::ptr::null(),
            0,
        );
        blst::blst_p2_compress(compressed.as_mut_ptr(), &point);
    }
    assert_eq!(result.return_data, compressed.to_vec(), "SBF output differs from blst");
}

#[test]
fn bench_syscalls() {
    let mollusk = mollusk();

    let g2_gen = G2Affine::generator().to_uncompressed();
    let g1_gen = G1Affine::generator().to_uncompressed();

    // Validate.
    let validate_g2 = run(&mollusk, 10, &g2_gen);
    println!(
        "g2 validate: rc={:?} {} CU",
        validate_g2.return_data.first(),
        validate_g2.compute_units_consumed
    );
    let validate_g1 = run(&mollusk, 18, &g1_gen);
    println!(
        "g1 validate: rc={:?} {} CU",
        validate_g1.return_data.first(),
        validate_g1.compute_units_consumed
    );

    // Add: gen + gen == 2*gen.
    let mut payload = g2_gen.to_vec();
    payload.extend_from_slice(&g2_gen);
    let add_g2 = run(&mollusk, 11, &payload);
    let expected = G2Affine::from(G2Projective::generator().double()).to_uncompressed();
    println!(
        "g2 add: rc={:?} match={} {} CU",
        add_g2.return_data.first(),
        add_g2.return_data.get(1..) == Some(expected.as_ref()),
        add_g2.compute_units_consumed
    );

    let mut payload = g1_gen.to_vec();
    payload.extend_from_slice(&g1_gen);
    let add_g1 = run(&mollusk, 15, &payload);
    let expected = G1Affine::from(G1Projective::generator().double()).to_uncompressed();
    println!(
        "g1 add: rc={:?} match={} {} CU",
        add_g1.return_data.first(),
        add_g1.return_data.get(1..) == Some(expected.as_ref()),
        add_g1.compute_units_consumed
    );

    // Mul: 7 * gen, scalar big-endian.
    let scalar_be = {
        let mut b = [0u8; 32];
        b[31] = 7;
        b
    };
    let mut payload = scalar_be.to_vec();
    payload.extend_from_slice(&g2_gen);
    let mul_g2 = run(&mollusk, 12, &payload);
    let expected =
        G2Affine::from(G2Projective::generator() * Scalar::from(7u64)).to_uncompressed();
    println!(
        "g2 mul: rc={:?} match={} {} CU",
        mul_g2.return_data.first(),
        mul_g2.return_data.get(1..) == Some(expected.as_ref()),
        mul_g2.compute_units_consumed
    );

    let mut payload = scalar_be.to_vec();
    payload.extend_from_slice(&g1_gen);
    let mul_g1 = run(&mollusk, 16, &payload);
    let expected =
        G1Affine::from(G1Projective::generator() * Scalar::from(7u64)).to_uncompressed();
    println!(
        "g1 mul: rc={:?} match={} {} CU",
        mul_g1.return_data.first(),
        mul_g1.return_data.get(1..) == Some(expected.as_ref()),
        mul_g1.compute_units_consumed
    );

    // Decompress.
    let decompress_g2 = run(&mollusk, 13, &G2Affine::generator().to_compressed());
    println!(
        "g2 decompress: rc={:?} match={} {} CU",
        decompress_g2.return_data.first(),
        decompress_g2.return_data.get(1..) == Some(g2_gen.as_ref()),
        decompress_g2.compute_units_consumed
    );
    let decompress_g1 = run(&mollusk, 17, &G1Affine::generator().to_compressed());
    println!(
        "g1 decompress: rc={:?} match={} {} CU",
        decompress_g1.return_data.first(),
        decompress_g1.return_data.get(1..) == Some(g1_gen.as_ref()),
        decompress_g1.compute_units_consumed
    );

    // Pairing, 1 and 2 pairs.
    let mut payload = g1_gen.to_vec();
    payload.extend_from_slice(&g2_gen);
    let pair_one = run(&mollusk, 14, &payload);
    println!(
        "pairing 1 pair: rc={:?} {} CU",
        pair_one.return_data.first(),
        pair_one.compute_units_consumed
    );

    let mut payload = g1_gen.to_vec();
    payload.extend_from_slice(&g1_gen);
    payload.extend_from_slice(&g2_gen);
    payload.extend_from_slice(&g2_gen);
    let pair_two = run(&mollusk, 14, &payload);
    println!(
        "pairing 2 pairs: rc={:?} {} CU",
        pair_two.return_data.first(),
        pair_two.compute_units_consumed
    );

    // big_mod_exp with 48-byte operands: base^exp mod p.
    let p_be = hex::decode(
        "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",
    )
    .unwrap();
    let mut payload = vec![0u8; 48];
    payload[47] = 5;
    let mut exp = vec![0u8; 48];
    exp[47] = 3;
    payload.extend_from_slice(&exp);
    payload.extend_from_slice(&p_be);
    let modexp = run(&mollusk, 20, &payload);
    println!(
        "big_mod_exp 48B: rc={:?} out_tail={:?} {} CU",
        modexp.return_data.first(),
        modexp.return_data.last(),
        modexp.compute_units_consumed
    );
}

#[test]
fn bench_u128_mac_loop() {
    let mollusk = mollusk();
    let base = run(&mollusk, 21, &0u64.to_le_bytes());
    let loop_100k = run(&mollusk, 21, &100_000u64.to_le_bytes());
    assert!(!base.program_result.is_err());
    assert!(!loop_100k.program_result.is_err());
    let per_op = (loop_100k.compute_units_consumed - base.compute_units_consumed) as f64 / 100_000.0;
    println!(
        "u128 mac: base={} loop={} per_op={:.2} CU",
        base.compute_units_consumed, loop_100k.compute_units_consumed, per_op
    );
}

// Per function CU costs, svm-unit-test style: each probe id runs one
// function in a loop-carried loop, so the number is that function alone.
#[test]
fn bench_function_costs() {
    let mollusk = mollusk();
    let probes: &[(u8, &str, u64)] = &[
        (0, "mont_mul", 20_000),
        (1, "mont_sqr", 20_000),
        (2, "from_mont (redc)", 20_000),
        (3, "to_mont", 20_000),
        (4, "add_mod", 100_000),
        (5, "sub_mod", 100_000),
        (6, "neg_mod", 100_000),
        (7, "mul2", 10_000),
        (8, "sq2", 10_000),
        (9, "add2", 50_000),
        (10, "sub2", 50_000),
        (11, "mul_by_xi2", 50_000),
        (12, "mul_by_a2i", 50_000),
        (13, "from_mont2", 10_000),
        (14, "to_mont2", 10_000),
        (15, "iso11_adapted", 1_000),
        (16, "iso3_adapted", 1_000),
        (17, "gx_at", 5_000),
        (18, "gx2_at", 2_000),
        (19, "mul_by_xi", 50_000),
        (20, "wit48 parse", 50_000),
        (21, "wit96 parse", 20_000),
        (22, "be/le round trip", 50_000),
        (23, "expand_message_xmd g1", 100),
        (24, "expand_message_xmd g2", 100),
        (25, "inv_divsteps", 200),
        (26, "fixed mul (a3)", 20_000),
        (27, "quotient check", 20_000),
        (28, "iso3_velu body", 1_000),
        (29, "pin (branch 1)", 2_000),
        (30, "pin (branch 2)", 2_000),
    ];
    println!("per function CU:");
    for &(id, name, n) in probes {
        let mut payload = vec![id];
        payload.extend_from_slice(&0u64.to_le_bytes());
        let base = run(&mollusk, 24, &payload);
        let mut payload = vec![id];
        payload.extend_from_slice(&n.to_le_bytes());
        let looped = run(&mollusk, 24, &payload);
        assert!(!looped.program_result.is_err(), "{name} probe failed");
        let per_op =
            (looped.compute_units_consumed - base.compute_units_consumed) as f64 / n as f64;
        println!("  {name:<22} {per_op:>10.1}");
    }
}

#[test]
fn bench_field_primitives() {
    let mollusk = mollusk();
    for (tag, name, n) in [(22u8, "mont_mul", 20_000u64), (23, "mul2", 10_000)] {
        let base = run(&mollusk, tag, &0u64.to_le_bytes());
        let looped = run(&mollusk, tag, &n.to_le_bytes());
        assert!(!base.program_result.is_err());
        assert!(!looped.program_result.is_err());
        let per_op =
            (looped.compute_units_consumed - base.compute_units_consumed) as f64 / n as f64;
        println!("{name}: per_op={per_op:.1} CU");
    }
}

// Stage by stage CU for the witnessed min-sig and min-pk pipelines.
#[test]
fn bench_witness_stage_breakdown() {
    let mollusk = mollusk();
    let names = ["hash_to_field", "sswu maps", "add + iso", "clear + validate"];

    for (tag, label, witness) in [
        (46u8, "min-sig G1", bls381_hash::witness::g1::generate(MESSAGE)),
        (47, "min-pk G2", bls381_hash::witness::g2::generate(MESSAGE)),
    ] {
        let mut cumulative = [0u64; 4];
        for stage in 0..4u8 {
            let mut payload = vec![stage];
            payload.extend_from_slice(&witness);
            payload.extend_from_slice(MESSAGE);
            let r = run(&mollusk, tag, &payload);
            assert!(!r.program_result.is_err(), "{label} stage {stage} failed");
            cumulative[stage as usize] = r.compute_units_consumed;
        }
        println!("witnessed {label} stages:");
        for (i, name) in names.iter().enumerate() {
            let delta = if i == 0 {
                cumulative[0]
            } else {
                cumulative[i] - cumulative[i - 1]
            };
            println!("  {name:<18} {delta:>7} CU (cumulative {})", cumulative[i]);
        }
    }
}

#[test]
fn bench_syscall_assisted_hash_to_g1() {
    use bls12_381::hash_to_curve::{HashToField, MapToCurve};

    let mollusk = mollusk();

    let field = run(&mollusk, 30, MESSAGE);
    assert!(!field.program_result.is_err(), "stage field: {:?}", field.program_result);

    let mapped = run(&mollusk, 31, MESSAGE);
    assert!(!mapped.program_result.is_err(), "stage maps: {:?}", mapped.program_result);

    let iso = run(&mollusk, 32, MESSAGE);
    assert!(!iso.program_result.is_err(), "stage iso: {:?}", iso.program_result);

    let full = run(&mollusk, 33, MESSAGE);
    assert!(!full.program_result.is_err(), "stage full: {:?}", full.program_result);

    // Reference: sum of the two mapped points before cofactor clearing.
    type F = <G1Projective as MapToCurve>::Field;
    let mut u = [F::default(); 2];
    F::hash_to_field::<ExpandMsgXmd<sha2::Sha256>>(MESSAGE, DST_G1, &mut u);
    let sum = G1Projective::map_to_curve(&u[0]) + G1Projective::map_to_curve(&u[1]);
    let expected_uncleared = G1Affine::from(sum).to_uncompressed();
    assert_eq!(
        iso.return_data,
        expected_uncleared.to_vec(),
        "pre-clearing point differs from zkcrypto"
    );

    // Reference: full hash_to_curve, zkcrypto and blst.
    let expected_full = G1Affine::from(
        <G1Projective as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::hash_to_curve(MESSAGE, DST_G1),
    )
    .to_uncompressed();
    assert_eq!(full.return_data, expected_full.to_vec(), "final point differs from zkcrypto");

    let mut point = blst::blst_p1::default();
    let mut serialized = [0u8; 96];
    unsafe {
        blst::blst_hash_to_g1(
            &mut point,
            MESSAGE.as_ptr(),
            MESSAGE.len(),
            DST_G1.as_ptr(),
            DST_G1.len(),
            std::ptr::null(),
            0,
        );
        blst::blst_p1_serialize(serialized.as_mut_ptr(), &point);
    }
    assert_eq!(full.return_data, serialized.to_vec(), "final point differs from blst");

    let f = field.compute_units_consumed;
    let m = mapped.compute_units_consumed;
    let i = iso.compute_units_consumed;
    let t = full.compute_units_consumed;
    println!();
    println!("syscall-assisted hash_to_G1 (min-sig):");
    println!("  hash_to_field:          {f} CU");
    println!("  2x sswu map:            {} CU", m.saturating_sub(f));
    println!("  E' add + iso-11:        {} CU", i.saturating_sub(m));
    println!("  clear_cofactor + check: {} CU", t.saturating_sub(i));
    println!("  TOTAL:                  {t} CU");
}

#[test]
fn bench_witness_hash_to_g1() {
    let mollusk = mollusk();

    let witnesses = bls381_hash::witness::g1::generate(MESSAGE);
    let mut payload = witnesses.clone();
    payload.extend_from_slice(MESSAGE);

    let result = run(&mollusk, 40, &payload);
    assert!(
        !result.program_result.is_err(),
        "witnessed hash_to_g1 failed: {:?}",
        result.program_result
    );

    let mut point = blst::blst_p1::default();
    let mut serialized = [0u8; 96];
    unsafe {
        blst::blst_hash_to_g1(
            &mut point,
            MESSAGE.as_ptr(),
            MESSAGE.len(),
            DST_G1.as_ptr(),
            DST_G1.len(),
            std::ptr::null(),
            0,
        );
        blst::blst_p1_serialize(serialized.as_mut_ptr(), &point);
    }
    assert_eq!(result.return_data, serialized.to_vec(), "differs from blst");

    println!(
        "witness-assisted hash_to_G1: {} CU ({} witness bytes)",
        result.compute_units_consumed,
        witnesses.len()
    );

    // corrupted witness must abort, not produce a different point
    let mut bad = payload.clone();
    bad[60] ^= 1;
    let rejected = run(&mollusk, 40, &bad);
    assert!(rejected.program_result.is_err(), "corrupt witness was accepted");
}

#[test]
fn bench_witness_hash_to_g2() {
    let mollusk = mollusk();

    let witnesses = bls381_hash::witness::g2::generate(MESSAGE);
    let mut payload = witnesses.clone();
    payload.extend_from_slice(MESSAGE);

    let result = run(&mollusk, 41, &payload);
    assert!(
        !result.program_result.is_err(),
        "witnessed hash_to_g2 failed: {:?}",
        result.program_result
    );

    let mut point = blst::blst_p2::default();
    let mut serialized = [0u8; 192];
    unsafe {
        blst::blst_hash_to_g2(
            &mut point,
            MESSAGE.as_ptr(),
            MESSAGE.len(),
            DST_G2.as_ptr(),
            DST_G2.len(),
            std::ptr::null(),
            0,
        );
        blst::blst_p2_serialize(serialized.as_mut_ptr(), &point);
    }
    assert_eq!(result.return_data, serialized.to_vec(), "differs from blst");

    println!(
        "witness-assisted hash_to_G2: {} CU ({} witness bytes)",
        result.compute_units_consumed,
        witnesses.len()
    );

    let mut bad = payload.clone();
    bad[120] ^= 1;
    let rejected = run(&mollusk, 41, &bad);
    assert!(rejected.program_result.is_err(), "corrupt witness was accepted");
}

fn blst_hash_g2_serialized(msg: &[u8]) -> [u8; 192] {
    let mut point = blst::blst_p2::default();
    let mut serialized = [0u8; 192];
    unsafe {
        blst::blst_hash_to_g2(
            &mut point,
            msg.as_ptr(),
            msg.len(),
            DST_G2.as_ptr(),
            DST_G2.len(),
            std::ptr::null(),
            0,
        );
        blst::blst_p2_serialize(serialized.as_mut_ptr(), &point);
    }
    serialized
}

#[test]
fn bench_witness_hash_to_g2_compact() {
    let mollusk = mollusk();

    let witnesses = bls381_hash::witness::g2::generate_compact(MESSAGE);
    let mut payload = witnesses.clone();
    payload.extend_from_slice(MESSAGE);

    let result = run(&mollusk, 48, &payload);
    assert!(
        !result.program_result.is_err(),
        "compact hash_to_g2 failed: {:?}",
        result.program_result
    );
    assert_eq!(result.return_data, blst_hash_g2_serialized(MESSAGE).to_vec(), "differs from blst");

    // and byte-identical to the default witnessed path
    let fat = bls381_hash::witness::g2::generate(MESSAGE);
    let mut fat_payload = fat.clone();
    fat_payload.extend_from_slice(MESSAGE);
    let fat_result = run(&mollusk, 41, &fat_payload);
    assert_eq!(result.return_data, fat_result.return_data, "compact and default paths disagree");

    println!(
        "compact witness-assisted hash_to_G2: {} CU ({} witness bytes)",
        result.compute_units_consumed,
        witnesses.len()
    );

    let mut bad = payload.clone();
    bad[100] ^= 1;
    let rejected = run(&mollusk, 48, &bad);
    assert!(rejected.program_result.is_err(), "corrupt compact witness was accepted");
}

#[test]
fn bench_witness_hash_to_g2_compact_xgcd() {
    let mollusk = mollusk();

    let witnesses = bls381_hash::witness::g2::generate_compact_xgcd(MESSAGE);
    assert_eq!(witnesses.len(), 97);
    let mut payload = witnesses.clone();
    payload.extend_from_slice(MESSAGE);

    let result = run(&mollusk, 49, &payload);
    assert!(
        !result.program_result.is_err(),
        "xgcd hash_to_g2 failed: {:?}",
        result.program_result
    );
    assert_eq!(result.return_data, blst_hash_g2_serialized(MESSAGE).to_vec(), "differs from blst");

    println!(
        "xgcd witness-assisted hash_to_G2: {} CU ({} witness bytes)",
        result.compute_units_consumed,
        witnesses.len()
    );

    let mut bad = payload.clone();
    bad[60] ^= 1;
    let rejected = run(&mollusk, 49, &bad);
    assert!(rejected.program_result.is_err(), "corrupt xgcd witness was accepted");
}

#[test]
fn bench_witness_hash_to_g2_compact_parity() {
    let mollusk = mollusk();

    let witnesses = bls381_hash::witness::g2::generate_compact_parity(MESSAGE);
    assert_eq!(witnesses.len(), 96);
    let mut payload = witnesses.clone();
    payload.extend_from_slice(MESSAGE);

    let result = run(&mollusk, 50, &payload);
    assert!(
        !result.program_result.is_err(),
        "parity hash_to_g2 failed: {:?}",
        result.program_result
    );
    assert_eq!(result.return_data, blst_hash_g2_serialized(MESSAGE).to_vec(), "differs from blst");

    println!(
        "parity witness-assisted hash_to_G2: {} CU ({} witness bytes)",
        result.compute_units_consumed,
        witnesses.len()
    );

    let mut bad = payload.clone();
    bad[60] ^= 1;
    let rejected = run(&mollusk, 50, &bad);
    assert!(rejected.program_result.is_err(), "corrupt parity witness was accepted");
}

// Zero-witness hash_to_G2 through big_mod_exp (SIMD-0529): the payload is
// the message alone, every root, character and inverse recomputed through
// the syscall. Pins blst byte-equality and the min-pk e2e verify.
#[test]
fn bench_modexp_hash_to_g2() {
    use blst::min_pk::{AggregatePublicKey, PublicKey, SecretKey, Signature};

    let mollusk = mollusk();

    let result = run(&mollusk, 57, MESSAGE);
    assert!(
        !result.program_result.is_err(),
        "modexp hash_to_g2 failed: {:?}",
        result.program_result
    );
    assert_eq!(result.return_data, blst_hash_g2_serialized(MESSAGE).to_vec(), "differs from blst");
    println!(
        "modexp zero-witness hash_to_G2: {} CU (0 witness bytes)",
        result.compute_units_consumed
    );

    let keys: Vec<SecretKey> = (0..20u8)
        .map(|i| {
            let ikm = [i + 1; 32];
            SecretKey::key_gen(&ikm, &[]).unwrap()
        })
        .collect();
    let pks: Vec<PublicKey> = keys.iter().map(|s| s.sk_to_pk()).collect();
    let all_refs: Vec<&PublicKey> = pks.iter().collect();
    let agg_all = AggregatePublicKey::aggregate(&all_refs, false)
        .unwrap()
        .to_public_key();
    let sigs: Vec<Signature> = keys
        .iter()
        .map(|s| s.sign(MESSAGE, DST_G2, &[]))
        .collect();
    let sig_refs: Vec<&Signature> = sigs.iter().collect();
    let agg_sig = blst::min_pk::AggregateSignature::aggregate(&sig_refs, false)
        .unwrap()
        .to_signature();

    let mut payload = vec![0u8];
    payload.extend_from_slice(&agg_all.serialize());
    payload.extend_from_slice(&agg_sig.compress());
    payload.extend_from_slice(MESSAGE);

    let result = run(&mollusk, 58, &payload);
    assert!(
        !result.program_result.is_err(),
        "modexp min-pk verify failed: {:?}",
        result.program_result
    );
    println!(
        "modexp min-pk end-to-end verify k=20: {} CU (0 witness bytes)",
        result.compute_units_consumed
    );

    // tampered signature must fail
    let mut bad = payload.clone();
    bad[1 + 96 + 10] ^= 1;
    let rejected = run(&mollusk, 58, &bad);
    assert!(rejected.program_result.is_err(), "tampered modexp verify was accepted");
}

// CU spread across messages: the divsteps batch count varies with the
// input (typical convergence ~27-28 of the 37-batch cap), so the parity
// hash cost moves a little per message. Also pins blst byte-equality
// across inputs, not just the fixture message.
#[test]
fn bench_parity_message_spread() {
    let mollusk = mollusk();
    let (mut lo, mut hi, mut sum) = (u64::MAX, 0u64, 0u64);
    let messages: Vec<Vec<u8>> = (0..8u8)
        .map(|i| format!("tapedrive vote payload: epoch {i}, slot {}", 1337 + i as u32).into_bytes())
        .collect();
    for msg in &messages {
        let mut payload = bls381_hash::witness::g2::generate_compact_parity(msg);
        payload.extend_from_slice(msg);
        let r = run(&mollusk, 50, &payload);
        assert!(!r.program_result.is_err(), "parity hash failed");
        assert_eq!(r.return_data, blst_hash_g2_serialized(msg).to_vec(), "differs from blst");
        let cu = r.compute_units_consumed;
        lo = lo.min(cu);
        hi = hi.max(cu);
        sum += cu;
    }
    println!(
        "parity hash CU over {} messages: min {lo} / avg {} / max {hi}",
        messages.len(),
        sum / messages.len() as u64
    );
}

// The parity layout's soundness sweep. No flags byte exists, so the flag
// probes drop out, and the other square root is no longer an equally valid
// witness: negating a root half flips its parity, which reads as a branch
// lie and must abort (the steered generator ships exactly those roots).
#[test]
fn witness_g2_parity_soundness() {
    let mollusk = mollusk();

    let witness = bls381_hash::witness::g2::generate_compact_parity(MESSAGE);
    let mut payload = witness.clone();
    payload.extend_from_slice(MESSAGE);

    let good = run(&mollusk, 50, &payload);
    assert!(!good.program_result.is_err(), "honest witness rejected");
    let truth = good.return_data.clone();

    for i in 0..witness.len() {
        let mut bad = payload.clone();
        bad[i] ^= 1;
        let r = run(&mollusk, 50, &bad);
        if !r.program_result.is_err() {
            assert_eq!(r.return_data, truth, "witness byte {i} steered the output");
        }
    }

    // a witness limb at or above the modulus must be rejected by the parser
    for start in [0usize, 48] {
        let mut oob = payload.clone();
        for byte in oob[start..start + 48].iter_mut() {
            *byte = 0xff;
        }
        assert!(
            run(&mollusk, 50, &oob).program_result.is_err(),
            "out-of-range witness accepted at byte {start}"
        );
    }

    // the witness is bound to the message it was generated for
    let mut replay = witness.clone();
    replay.extend_from_slice(b"a different snapshot vote payload");
    assert!(run(&mollusk, 50, &replay).program_result.is_err(), "cross-message replay accepted");

    // the other root of either map is a parity lie, not a valid witness
    for steer in [1u8, 2, 3] {
        let mut lied = bls381_hash::witness::g2::generate_compact_parity_steered(MESSAGE, steer);
        lied.extend_from_slice(MESSAGE);
        assert!(
            run(&mollusk, 50, &lied).program_result.is_err(),
            "flipped root (parity lie {steer:#04b}) accepted"
        );
    }
}

// The compact-blob counterpart of witness_g2_soundness, shared by both
// layouts: no bit of the blob can steer the output, the flag range is
// canonical, witnesses are message-bound, the other square root does not
// steer, and a branch lie with a self-consistent batch still aborts.
fn compact_soundness_sweep(tag: u8, witness: Vec<u8>, steered: fn(&[u8], u8) -> Vec<u8>) {
    let mollusk = mollusk();

    let mut payload = witness.clone();
    payload.extend_from_slice(MESSAGE);

    let good = run(&mollusk, tag, &payload);
    assert!(!good.program_result.is_err(), "honest witness rejected");
    let truth = good.return_data.clone();

    for i in 0..witness.len() {
        let mut bad = payload.clone();
        bad[i] ^= 1;
        let r = run(&mollusk, tag, &bad);
        if !r.program_result.is_err() {
            assert_eq!(r.return_data, truth, "witness byte {i} steered the output");
        }
    }
    for bit in 1..8 {
        let mut bad = payload.clone();
        bad[0] ^= 1 << bit;
        let r = run(&mollusk, tag, &bad);
        if !r.program_result.is_err() {
            assert_eq!(r.return_data, truth, "flag bit {bit} steered the output");
        }
    }

    // non-canonical flags byte
    let mut bad = payload.clone();
    bad[0] = 4;
    assert!(run(&mollusk, tag, &bad).program_result.is_err(), "flags=4 accepted");

    // a witness limb at or above the modulus must be rejected by the parser
    for start in (1..witness.len()).step_by(48) {
        let mut oob = payload.clone();
        for byte in oob[start..start + 48].iter_mut() {
            *byte = 0xff;
        }
        assert!(
            run(&mollusk, tag, &oob).program_result.is_err(),
            "out-of-range witness accepted at byte {start}"
        );
    }

    // the witness is bound to the message it was generated for
    let mut replay = witness.clone();
    replay.extend_from_slice(b"a different snapshot vote payload");
    assert!(run(&mollusk, tag, &replay).program_result.is_err(), "cross-message replay accepted");

    // the other square root is an equally valid witness and must not steer
    let mut alt = bls381_hash::witness::g2::flip_first_root(&witness);
    alt.extend_from_slice(MESSAGE);
    let same = run(&mollusk, tag, &alt);
    assert!(!same.program_result.is_err(), "flipped root rejected");
    assert_eq!(same.return_data, truth, "flipped root changed the point");

    // a branch lie whose batch stays consistent with the lie still has no
    // satisfiable sqrt check: wrong-branch gx is a non-square
    for steer in [1u8, 2, 3] {
        let mut lied = steered(MESSAGE, steer);
        lied.extend_from_slice(MESSAGE);
        assert!(
            run(&mollusk, tag, &lied).program_result.is_err(),
            "steered branch flags {steer:#04b} accepted"
        );
    }
}

#[test]
fn witness_g2_compact_soundness() {
    compact_soundness_sweep(
        48,
        bls381_hash::witness::g2::generate_compact(MESSAGE),
        bls381_hash::witness::g2::generate_compact_steered,
    );
}

#[test]
fn witness_g2_xgcd_soundness() {
    compact_soundness_sweep(
        49,
        bls381_hash::witness::g2::generate_compact_xgcd(MESSAGE),
        bls381_hash::witness::g2::generate_compact_xgcd_steered,
    );
}

// No witness byte can steer the output: every single-bit corruption of the G2
// witness must abort or reproduce the exact same point. Also covers the branch
// flag range, canonical-form parsing, and message binding.
#[test]
fn witness_g2_soundness() {
    let mollusk = mollusk();

    let witness = bls381_hash::witness::g2::generate(MESSAGE);
    let mut payload = witness.clone();
    payload.extend_from_slice(MESSAGE);

    let good = run(&mollusk, 41, &payload);
    assert!(!good.program_result.is_err(), "honest witness rejected");
    let truth = good.return_data.clone();

    for i in 0..witness.len() {
        let mut bad = payload.clone();
        bad[i] ^= 1;
        let r = run(&mollusk, 41, &bad);
        if !r.program_result.is_err() {
            assert_eq!(r.return_data, truth, "witness byte {i} steered the output");
        }
    }

    // a branch flag above 1 is non-canonical
    for flag in [0usize, 97] {
        let mut bad = payload.clone();
        bad[flag] = 2;
        assert!(run(&mollusk, 41, &bad).program_result.is_err(), "flag=2 at {flag} accepted");
    }

    // a witness limb at or above the modulus must be rejected by the parser
    let mut oob = payload.clone();
    for byte in oob[1..49].iter_mut() {
        *byte = 0xff;
    }
    assert!(run(&mollusk, 41, &oob).program_result.is_err(), "out-of-range witness accepted");

    // the witness is bound to the message it was generated for
    let mut replay = witness.clone();
    replay.extend_from_slice(b"a different snapshot vote payload");
    assert!(run(&mollusk, 41, &replay).program_result.is_err(), "cross-message replay accepted");

    // the other square root is an equally valid witness and must not steer
    let mut alt = bls381_hash::witness::g2::flip_first_sqrt(&witness);
    alt.extend_from_slice(MESSAGE);
    let same = run(&mollusk, 41, &alt);
    assert!(!same.program_result.is_err(), "flipped root rejected");
    assert_eq!(same.return_data, truth, "flipped root changed the point");
}

// The G1 (min-sig) counterpart of the G2 soundness sweep.
#[test]
fn witness_g1_soundness() {
    let mollusk = mollusk();

    let witness = bls381_hash::witness::g1::generate(MESSAGE);
    let mut payload = witness.clone();
    payload.extend_from_slice(MESSAGE);

    let good = run(&mollusk, 40, &payload);
    assert!(!good.program_result.is_err(), "honest witness rejected");
    let truth = good.return_data.clone();

    for i in 0..witness.len() {
        let mut bad = payload.clone();
        bad[i] ^= 1;
        let r = run(&mollusk, 40, &bad);
        if !r.program_result.is_err() {
            assert_eq!(r.return_data, truth, "witness byte {i} steered the output");
        }
    }

    for flag in [0usize, 97] {
        let mut bad = payload.clone();
        bad[flag] = 2;
        assert!(run(&mollusk, 40, &bad).program_result.is_err(), "flag=2 at {flag} accepted");
    }

    let mut oob = payload.clone();
    for byte in oob[1..49].iter_mut() {
        *byte = 0xff;
    }
    assert!(run(&mollusk, 40, &oob).program_result.is_err(), "out-of-range witness accepted");

    let mut replay = witness.clone();
    replay.extend_from_slice(b"a different snapshot vote payload");
    assert!(run(&mollusk, 40, &replay).program_result.is_err(), "cross-message replay accepted");

    // the other square root is an equally valid witness and must not steer
    let mut alt = bls381_hash::witness::g1::flip_first_sqrt(&witness);
    alt.extend_from_slice(MESSAGE);
    let same = run(&mollusk, 40, &alt);
    assert!(!same.program_result.is_err(), "flipped root rejected");
    assert_eq!(same.return_data, truth, "flipped root changed the point");
}

#[test]
fn bench_witness_nu_encode() {
    let mollusk = mollusk();
    const DST_G1_NU: &[u8] = b"BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_NU_POP_";
    const DST_G2_NU: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_NU_POP_";

    let payload = bls381_hash::witness::g1::generate_nu(MESSAGE);
    let r1 = run(&mollusk, 44, &payload);
    assert!(!r1.program_result.is_err(), "nu g1: {:?}", r1.program_result);
    let mut pt = blst::blst_p1::default();
    let mut ser = [0u8; 96];
    unsafe {
        blst::blst_encode_to_g1(&mut pt, MESSAGE.as_ptr(), MESSAGE.len(), DST_G1_NU.as_ptr(), DST_G1_NU.len(), std::ptr::null(), 0);
        blst::blst_p1_serialize(ser.as_mut_ptr(), &pt);
    }
    assert_eq!(r1.return_data, ser.to_vec(), "nu g1 differs from blst encode");
    println!("witness-assisted NU encode_to_G1: {} CU ({} witness bytes)", r1.compute_units_consumed, payload.len() - MESSAGE.len());
    let mut bad = payload.clone();
    bad[60] ^= 1;
    assert!(run(&mollusk, 44, &bad).program_result.is_err(), "corrupt nu g1 accepted");

    let payload = bls381_hash::witness::g2::generate_nu(MESSAGE);
    let r2 = run(&mollusk, 45, &payload);
    assert!(!r2.program_result.is_err(), "nu g2: {:?}", r2.program_result);
    let mut pt = blst::blst_p2::default();
    let mut ser = [0u8; 192];
    unsafe {
        blst::blst_encode_to_g2(&mut pt, MESSAGE.as_ptr(), MESSAGE.len(), DST_G2_NU.as_ptr(), DST_G2_NU.len(), std::ptr::null(), 0);
        blst::blst_p2_serialize(ser.as_mut_ptr(), &pt);
    }
    assert_eq!(r2.return_data, ser.to_vec(), "nu g2 differs from blst encode");
    println!("witness-assisted NU encode_to_G2: {} CU ({} witness bytes)", r2.compute_units_consumed, payload.len() - MESSAGE.len());
    let mut bad = payload.clone();
    bad[120] ^= 1;
    assert!(run(&mollusk, 45, &bad).program_result.is_err(), "corrupt nu g2 accepted");
}

// Correctness guard: the bare Montgomery reduction (from_mont) against the
// general multiply, plus the adapted iso-11 chain against Horner. Host-side.
#[test]
fn field_arithmetic_selftest() {
    bls381_hash::witness::g1::iso11_adapted_selftest();
    bls381_hash::witness::g1::redc_selftest();
    bls381_hash::witness::g1::inv_divsteps_selftest();
    bls381_hash::witness::g1::fixed_selftest();
}

#[test]
fn bench_min_pk_verify_end_to_end() {
    use blst::min_pk::{AggregatePublicKey, AggregateSignature, PublicKey, SecretKey, Signature};

    let mollusk = mollusk();

    let keys: Vec<SecretKey> = (0..20u8)
        .map(|i| {
            let ikm = [i + 1; 32];
            SecretKey::key_gen(&ikm, &[]).unwrap()
        })
        .collect();
    let pks: Vec<PublicKey> = keys.iter().map(|s| s.sk_to_pk()).collect();
    let all_refs: Vec<&PublicKey> = pks.iter().collect();
    let agg_all = AggregatePublicKey::aggregate(&all_refs, false)
        .unwrap()
        .to_public_key();

    let witness = bls381_hash::witness::g2::generate(MESSAGE);

    for k in [14usize, 20] {
        let sigs: Vec<Signature> = keys[..k]
            .iter()
            .map(|s| s.sign(MESSAGE, DST_G2, &[]))
            .collect();
        let sig_refs: Vec<&Signature> = sigs.iter().collect();
        let agg_sig = AggregateSignature::aggregate(&sig_refs, false)
            .unwrap()
            .to_signature();

        let mut payload = vec![(20 - k) as u8];
        payload.extend_from_slice(&agg_all.serialize());
        payload.extend_from_slice(&agg_sig.compress());
        for pk in &pks[k..] {
            payload.extend_from_slice(&pk.compress());
        }
        payload.extend_from_slice(&witness);
        payload.extend_from_slice(MESSAGE);

        let result = run(&mollusk, 51, &payload);
        assert!(
            !result.program_result.is_err(),
            "min-pk verify failed at k={k}: {:?}",
            result.program_result
        );
        println!(
            "min-pk end-to-end verify k={k}: {} CU",
            result.compute_units_consumed
        );

        // tampered signature must fail
        let mut bad = payload.clone();
        bad[1 + 96 + 10] ^= 1;
        let rejected = run(&mollusk, 51, &bad);
        assert!(rejected.program_result.is_err(), "forged min-pk verify accepted at k={k}");
    }
}

// The same end-to-end min-pk verify against the 145-byte compact witness
// (tag 54), the 97-byte witness-free-inverse blob (tag 55) and the
// 96-byte parity blob (tag 56).
#[test]
fn bench_min_pk_verify_compact_end_to_end() {
    use blst::min_pk::{AggregatePublicKey, AggregateSignature, PublicKey, SecretKey, Signature};

    let mollusk = mollusk();

    let keys: Vec<SecretKey> = (0..20u8)
        .map(|i| {
            let ikm = [i + 1; 32];
            SecretKey::key_gen(&ikm, &[]).unwrap()
        })
        .collect();
    let pks: Vec<PublicKey> = keys.iter().map(|s| s.sk_to_pk()).collect();
    let all_refs: Vec<&PublicKey> = pks.iter().collect();
    let agg_all = AggregatePublicKey::aggregate(&all_refs, false)
        .unwrap()
        .to_public_key();

    let witness = bls381_hash::witness::g2::generate_compact(MESSAGE);
    let witness_xgcd = bls381_hash::witness::g2::generate_compact_xgcd(MESSAGE);
    let witness_parity = bls381_hash::witness::g2::generate_compact_parity(MESSAGE);

    for k in [14usize, 20] {
        let sigs: Vec<Signature> = keys[..k]
            .iter()
            .map(|s| s.sign(MESSAGE, DST_G2, &[]))
            .collect();
        let sig_refs: Vec<&Signature> = sigs.iter().collect();
        let agg_sig = AggregateSignature::aggregate(&sig_refs, false)
            .unwrap()
            .to_signature();

        let mut head = vec![(20 - k) as u8];
        head.extend_from_slice(&agg_all.serialize());
        head.extend_from_slice(&agg_sig.compress());
        for pk in &pks[k..] {
            head.extend_from_slice(&pk.compress());
        }

        for (tag, label, wit) in [
            (54u8, "compact", &witness),
            (55, "xgcd", &witness_xgcd),
            (56, "parity", &witness_parity),
        ] {
            let mut payload = head.clone();
            payload.extend_from_slice(wit);
            payload.extend_from_slice(MESSAGE);

            let result = run(&mollusk, tag, &payload);
            assert!(
                !result.program_result.is_err(),
                "{label} min-pk verify failed at k={k}: {:?}",
                result.program_result
            );
            println!(
                "{label} min-pk end-to-end verify k={k}: {} CU ({} witness bytes)",
                result.compute_units_consumed,
                wit.len()
            );

            // tampered signature must fail
            let mut bad = payload.clone();
            bad[1 + 96 + 10] ^= 1;
            let rejected = run(&mollusk, tag, &bad);
            assert!(rejected.program_result.is_err(), "forged {label} min-pk verify accepted at k={k}");
        }
    }
}

fn fat_payload(msg: &[u8]) -> Vec<u8> {
    let mut payload = bls381_hash::witness::g2::generate_fat(msg);
    payload.extend_from_slice(msg);
    payload
}

#[test]
fn bench_witness_hash_to_g2_fat() {
    let mollusk = mollusk();

    let payload = fat_payload(MESSAGE);
    let result = run(&mollusk, 60, &payload);
    assert!(
        !result.program_result.is_err(),
        "fat hash_to_g2 failed: {:?}",
        result.program_result
    );
    assert_eq!(result.return_data, blst_hash_g2_serialized(MESSAGE).to_vec(), "differs from blst");
    println!(
        "fat witness-assisted hash_to_G2: {} CU ({} witness bytes)",
        result.compute_units_consumed,
        payload.len() - MESSAGE.len()
    );

    let mut bad = payload.clone();
    bad[200] ^= 1;
    assert!(run(&mollusk, 60, &bad).program_result.is_err(), "corrupt fat witness was accepted");
}

#[test]
fn bench_fat_stage_breakdown() {
    let mollusk = mollusk();
    let names = ["hash_to_field", "x pins", "roots + iso + add", "clear + validate"];
    let witness = bls381_hash::witness::g2::generate_fat(MESSAGE);
    let mut cumulative = [0u64; 4];
    for stage in 0..4u8 {
        let mut payload = vec![stage];
        payload.extend_from_slice(&witness);
        payload.extend_from_slice(MESSAGE);
        let r = run(&mollusk, 62, &payload);
        assert!(!r.program_result.is_err(), "fat stage {stage} failed");
        cumulative[stage as usize] = r.compute_units_consumed;
    }
    println!("fat min-pk G2 stages (flags {:#04b}):", witness[0]);
    for (i, name) in names.iter().enumerate() {
        let delta = if i == 0 { cumulative[0] } else { cumulative[i] - cumulative[i - 1] };
        println!("  {name:<18} {delta:>7} CU (cumulative {})", cumulative[i]);
    }
}

// Default and fat layouts side by side over the same messages: branch
// flags vary, so both the mean and the spread matter.
#[test]
fn bench_fat_message_spread() {
    let mollusk = mollusk();
    let messages: Vec<Vec<u8>> = (0..16u8)
        .map(|i| format!("tapedrive vote payload: epoch {i}, slot {}", 1337 + i as u32).into_bytes())
        .collect();
    for (tag, label) in [(41u8, "default"), (60, "fat")] {
        let (mut lo, mut hi, mut sum) = (u64::MAX, 0u64, 0u64);
        for msg in &messages {
            let mut payload = if tag == 60 {
                bls381_hash::witness::g2::generate_fat(msg)
            } else {
                bls381_hash::witness::g2::generate(msg)
            };
            payload.extend_from_slice(msg);
            let r = run(&mollusk, tag, &payload);
            assert!(!r.program_result.is_err(), "{label} hash failed");
            assert_eq!(r.return_data, blst_hash_g2_serialized(msg).to_vec(), "{label} differs from blst");
            let cu = r.compute_units_consumed;
            lo = lo.min(cu);
            hi = hi.max(cu);
            sum += cu;
        }
        println!(
            "{label} hash CU over {} messages: min {lo} / avg {} / max {hi}",
            messages.len(),
            sum / messages.len() as u64
        );
    }
}

// No witness bit can steer the fat output. The flags byte is canonical,
// witnesses are message-bound, either root of either map gives the same
// point, and a branch lie (the other SSWU candidate, pinned consistently
// with the lied flag) cannot land on E.
#[test]
fn witness_g2_fat_soundness() {
    let mollusk = mollusk();

    let witness = bls381_hash::witness::g2::generate_fat(MESSAGE);
    let mut payload = witness.clone();
    payload.extend_from_slice(MESSAGE);

    let good = run(&mollusk, 60, &payload);
    assert!(!good.program_result.is_err(), "honest witness rejected");
    let truth = good.return_data.clone();

    for i in 0..witness.len() {
        let mut bad = payload.clone();
        bad[i] ^= 1;
        let r = run(&mollusk, 60, &bad);
        if !r.program_result.is_err() {
            assert_eq!(r.return_data, truth, "witness byte {i} steered the output");
        }
    }
    for bit in 2..8 {
        let mut bad = payload.clone();
        bad[0] ^= 1 << bit;
        assert!(run(&mollusk, 60, &bad).program_result.is_err(), "flag bit {bit} accepted");
    }

    for start in (1..witness.len()).step_by(48) {
        let mut oob = payload.clone();
        for byte in oob[start..start + 48].iter_mut() {
            *byte = 0xff;
        }
        assert!(
            run(&mollusk, 60, &oob).program_result.is_err(),
            "out-of-range witness accepted at byte {start}"
        );
    }

    let mut replay = witness.clone();
    replay.extend_from_slice(b"a different snapshot vote payload");
    assert!(run(&mollusk, 60, &replay).program_result.is_err(), "cross-message replay accepted");

    for map in 0..2 {
        let mut alt = bls381_hash::witness::g2::flip_fat_root(&witness, map);
        alt.extend_from_slice(MESSAGE);
        let same = run(&mollusk, 60, &alt);
        assert!(!same.program_result.is_err(), "flipped root {map} rejected");
        assert_eq!(same.return_data, truth, "flipped root {map} changed the point");
    }

    for steer in [1u8, 2, 3] {
        let mut lied = bls381_hash::witness::g2::generate_fat_steered(MESSAGE, steer);
        lied.extend_from_slice(MESSAGE);
        assert!(
            run(&mollusk, 60, &lied).program_result.is_err(),
            "steered branch flags {steer:#04b} accepted"
        );
    }
}

#[test]
fn bench_min_pk_verify_fat_end_to_end() {
    use blst::min_pk::{AggregatePublicKey, AggregateSignature, PublicKey, SecretKey, Signature};

    let mollusk = mollusk();
    let keys: Vec<SecretKey> = (0..20u8)
        .map(|i| SecretKey::key_gen(&[i + 1; 32], &[]).unwrap())
        .collect();
    let pks: Vec<PublicKey> = keys.iter().map(|s| s.sk_to_pk()).collect();
    let all_refs: Vec<&PublicKey> = pks.iter().collect();
    let agg_all = AggregatePublicKey::aggregate(&all_refs, false).unwrap().to_public_key();
    let witness = bls381_hash::witness::g2::generate_fat(MESSAGE);

    for k in [14usize, 20] {
        let sigs: Vec<Signature> = keys[..k].iter().map(|s| s.sign(MESSAGE, DST_G2, &[])).collect();
        let sig_refs: Vec<&Signature> = sigs.iter().collect();
        let agg_sig = AggregateSignature::aggregate(&sig_refs, false).unwrap().to_signature();

        let mut payload = vec![(20 - k) as u8];
        payload.extend_from_slice(&agg_all.serialize());
        payload.extend_from_slice(&agg_sig.compress());
        for pk in &pks[k..] {
            payload.extend_from_slice(&pk.compress());
        }
        payload.extend_from_slice(&witness);
        payload.extend_from_slice(MESSAGE);

        let result = run(&mollusk, 61, &payload);
        assert!(!result.program_result.is_err(), "fat min-pk verify failed at k={k}: {:?}", result.program_result);
        println!("fat min-pk end-to-end verify k={k}: {} CU", result.compute_units_consumed);

        let mut bad = payload.clone();
        bad[1 + 96 + 10] ^= 1;
        assert!(run(&mollusk, 61, &bad).program_result.is_err(), "forged fat min-pk verify accepted at k={k}");
    }
}

// The same verify with the signature and absentee keys uncompressed (tag
// 63): the v1 byte budget buys back every decompress syscall.
#[test]
fn bench_min_pk_verify_fat_uncompressed() {
    use blst::min_pk::{AggregatePublicKey, AggregateSignature, PublicKey, SecretKey, Signature};

    let mollusk = mollusk();
    let keys: Vec<SecretKey> = (0..20u8)
        .map(|i| SecretKey::key_gen(&[i + 1; 32], &[]).unwrap())
        .collect();
    let pks: Vec<PublicKey> = keys.iter().map(|s| s.sk_to_pk()).collect();
    let all_refs: Vec<&PublicKey> = pks.iter().collect();
    let agg_all = AggregatePublicKey::aggregate(&all_refs, false).unwrap().to_public_key();
    let witness = bls381_hash::witness::g2::generate_fat(MESSAGE);

    for k in [14usize, 20] {
        let sigs: Vec<Signature> = keys[..k].iter().map(|s| s.sign(MESSAGE, DST_G2, &[])).collect();
        let sig_refs: Vec<&Signature> = sigs.iter().collect();
        let agg_sig = AggregateSignature::aggregate(&sig_refs, false).unwrap().to_signature();

        let mut payload = vec![(20 - k) as u8];
        payload.extend_from_slice(&agg_all.serialize());
        payload.extend_from_slice(&agg_sig.serialize());
        for pk in &pks[k..] {
            payload.extend_from_slice(&pk.serialize());
        }
        payload.extend_from_slice(&witness);
        payload.extend_from_slice(MESSAGE);

        let result = run(&mollusk, 63, &payload);
        assert!(!result.program_result.is_err(), "uncompressed fat verify failed at k={k}: {:?}", result.program_result);
        println!("fat min-pk verify, uncompressed inputs, k={k}: {} CU", result.compute_units_consumed);

        let mut bad = payload.clone();
        bad[1 + 96 + 100] ^= 1;
        assert!(run(&mollusk, 63, &bad).program_result.is_err(), "tampered uncompressed verify accepted at k={k}");
    }
}

// Byte equality with blst across many messages (both branch bits, all
// sign combinations), the fat layout's broad correctness guard.
#[test]
fn fat_matches_blst_across_messages() {
    let mollusk = mollusk();
    let mut flags_seen = [0usize; 4];
    for i in 0..256u32 {
        let msg = format!("fat sweep message {i}: epoch {}, slot {}", i / 7, 1000 + i).into_bytes();
        let payload = fat_payload(&msg);
        flags_seen[payload[0] as usize] += 1;
        let r = run(&mollusk, 60, &payload);
        assert!(!r.program_result.is_err(), "fat hash failed on message {i}");
        assert_eq!(r.return_data, blst_hash_g2_serialized(&msg).to_vec(), "fat differs from blst on message {i}");
    }
    println!("fat sweep: 256 messages blst-equal, branch flags seen {flags_seen:?}");
    assert!(flags_seen.iter().all(|&n| n > 0), "sweep missed a branch combination");
}

// Empty and identity keys or signatures must never verify. All-zero bytes
// fail the curve syscalls outright (no compression flag, or (0, 0) off the
// curve). The identity encodings (0x40 / 0xc0 then zeros) decode cleanly, and
// an identity key with an identity signature satisfies the pairing for any
// message, so the verifiers screen them on the encoding. Covers every verify
// tag: the stored aggregate as identity, a committee whose absentees cancel
// it, an identity absentee, and zero bytes in each slot.
#[test]
fn verify_rejects_empty_and_identity_points() {
    use blst::min_pk::{AggregatePublicKey, AggregateSignature, PublicKey, SecretKey, Signature};

    let mollusk = mollusk();
    let keys: Vec<SecretKey> = (0..20u8).map(|i| SecretKey::key_gen(&[i + 1; 32], &[]).unwrap()).collect();
    let pks: Vec<PublicKey> = keys.iter().map(|s| s.sk_to_pk()).collect();
    let refs: Vec<&PublicKey> = pks.iter().collect();
    let agg = AggregatePublicKey::aggregate(&refs, false).unwrap().to_public_key();
    let sigs: Vec<Signature> = keys.iter().map(|s| s.sign(MESSAGE, DST_G2, &[])).collect();
    let sig_refs: Vec<&Signature> = sigs.iter().collect();
    let agg_sig = AggregateSignature::aggregate(&sig_refs, false).unwrap().to_signature();

    let mut id_g1 = [0u8; 96];
    id_g1[0] = 0x40;
    let mut id_g1c = [0u8; 48];
    id_g1c[0] = 0xc0;
    let mut id_g2 = [0u8; 192];
    id_g2[0] = 0x40;
    let mut id_g2c = [0u8; 96];
    id_g2c[0] = 0xc0;

    // compressed-signature tags and their witness generators
    type Gen = fn(&[u8]) -> Vec<u8>;
    let tags: [(u8, Gen); 6] = [
        (51, bls381_hash::witness::g2::generate),
        (54, bls381_hash::witness::g2::generate_compact),
        (55, bls381_hash::witness::g2::generate_compact_xgcd),
        (56, bls381_hash::witness::g2::generate_compact_parity),
        (58, |_| Vec::new()),
        (61, bls381_hash::witness::g2::generate_fat),
    ];
    for (tag, gen) in tags {
        let witness = gen(MESSAGE);
        let build = |agg: &[u8], sig: &[u8], absent: &[Vec<u8>]| {
            let mut p = vec![absent.len() as u8];
            p.extend_from_slice(agg);
            p.extend_from_slice(sig);
            for a in absent {
                p.extend_from_slice(a);
            }
            p.extend_from_slice(&witness);
            p.extend_from_slice(MESSAGE);
            p
        };
        let honest = build(&agg.serialize(), &agg_sig.compress(), &[]);
        assert!(!run(&mollusk, tag, &honest).program_result.is_err(), "tag {tag}: honest verify rejected");

        let every: Vec<Vec<u8>> = pks.iter().map(|pk| pk.compress().to_vec()).collect();
        let cases: [(&str, Vec<u8>); 6] = [
            ("identity aggregate + identity signature", build(&id_g1, &id_g2c, &[])),
            ("absentees cancel the aggregate + identity signature", build(&agg.serialize(), &id_g2c, &every)),
            ("identity absentee", build(&agg.serialize(), &agg_sig.compress(), &[id_g1c.to_vec()])),
            ("zero aggregate + identity signature", build(&[0u8; 96], &id_g2c, &[])),
            ("zero signature", build(&agg.serialize(), &[0u8; 96], &[])),
            ("zero absentee", build(&agg.serialize(), &agg_sig.compress(), &[vec![0u8; 48]])),
        ];
        for (label, payload) in cases {
            assert!(run(&mollusk, tag, &payload).program_result.is_err(), "tag {tag}: {label} accepted");
        }
    }

    // tag 63: uncompressed signature and absentee keys
    let witness = bls381_hash::witness::g2::generate_fat(MESSAGE);
    let build = |agg: &[u8], sig: &[u8], absent: &[Vec<u8>]| {
        let mut p = vec![absent.len() as u8];
        p.extend_from_slice(agg);
        p.extend_from_slice(sig);
        for a in absent {
            p.extend_from_slice(a);
        }
        p.extend_from_slice(&witness);
        p.extend_from_slice(MESSAGE);
        p
    };
    assert!(!run(&mollusk, 63, &build(&agg.serialize(), &agg_sig.serialize(), &[])).program_result.is_err());
    let every: Vec<Vec<u8>> = pks.iter().map(|pk| pk.serialize().to_vec()).collect();
    let cases: [(&str, Vec<u8>); 6] = [
        ("identity aggregate + identity signature", build(&id_g1, &id_g2, &[])),
        ("absentees cancel the aggregate + identity signature", build(&agg.serialize(), &id_g2, &every)),
        ("identity absentee", build(&agg.serialize(), &agg_sig.serialize(), &[id_g1.to_vec()])),
        ("zero aggregate + identity signature", build(&[0u8; 96], &id_g2, &[])),
        ("zero signature", build(&agg.serialize(), &[0u8; 192], &[])),
        ("zero absentee", build(&agg.serialize(), &agg_sig.serialize(), &[vec![0u8; 96]])),
    ];
    for (label, payload) in cases {
        assert!(run(&mollusk, 63, &payload).program_result.is_err(), "tag 63: {label} accepted");
    }
}

// The blst-backed fat generator must match the portable one byte for byte:
// the verifier accepts one blob per message, so any drift is a broken witness.
#[test]
fn fat_blst_generator_matches_portable() {
    let mut flags_seen = [0usize; 4];
    for i in 0..512u32 {
        let msg = format!("blst witness check {i}: epoch {}, slot {}", i / 3, 9000 + i).into_bytes();
        let portable = bls381_hash::witness::g2::generate_fat(&msg);
        let fast = bls381_hash::witness::g2::generate_fat_blst(&msg);
        assert_eq!(fast, portable, "blst generator differs on message {i}");
        flags_seen[portable[0] as usize] += 1;
    }
    for msg in [&b""[..], MESSAGE, &[0xffu8; 300][..]] {
        assert_eq!(
            bls381_hash::witness::g2::generate_fat_blst(msg),
            bls381_hash::witness::g2::generate_fat(msg)
        );
    }
    assert!(flags_seen.iter().all(|&n| n > 0), "sweep missed a branch combination");

    // and the program accepts it
    let mollusk = mollusk();
    let mut payload = bls381_hash::witness::g2::generate_fat_blst(MESSAGE);
    payload.extend_from_slice(MESSAGE);
    let r = run(&mollusk, 60, &payload);
    assert!(!r.program_result.is_err());
    assert_eq!(r.return_data, blst_hash_g2_serialized(MESSAGE).to_vec());
}

fn refused_cleanly(r: &InstructionResult) -> bool {
    format!("{:?}", r.program_result) == "Failure(InvalidInstructionData)"
}

// Out-of-range stages are refused, the last stage is the full validated
// hash, and short or empty payloads get a clean error on every entry point
// (a panic would surface as an unknown error instead).
#[test]
fn fat_stage_bounds_and_short_payloads() {
    let mollusk = mollusk();
    let witness = bls381_hash::witness::g2::generate_fat(MESSAGE);
    let mut full = witness.clone();
    full.extend_from_slice(MESSAGE);
    let hash = run(&mollusk, 60, &full);
    assert!(!hash.program_result.is_err());

    let staged = |stage: u8| {
        let mut p = vec![stage];
        p.extend_from_slice(&full);
        run(&mollusk, 62, &p)
    };
    assert_eq!(staged(3).return_data, hash.return_data, "stage 3 is the full hash");
    for stage in [4u8, 5, 200, 255] {
        assert!(refused_cleanly(&staged(stage)), "stage {stage} not refused");
    }

    for tag in [51u8, 54, 55, 56, 58, 60, 61, 62, 63] {
        assert!(refused_cleanly(&run(&mollusk, tag, &[])), "tag {tag}: empty payload");
    }
    for tag in [51u8, 54, 55, 56, 58, 61, 63] {
        // one absentee announced, none present
        assert!(refused_cleanly(&run(&mollusk, tag, &[1u8; 150])), "tag {tag}: short payload");
    }
}
