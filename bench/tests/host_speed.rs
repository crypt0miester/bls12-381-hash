//! Host (off-chain) speed: hash-to-G2 in each implementation, this crate's
//! witness generators, the field multiply, and min-pk sign/verify. Run with
//! cargo test --release --test host_speed -- --ignored --nocapture

use std::hint::black_box;
use std::time::Instant;

const DST: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_";
const MESSAGE: &[u8] = b"tapedrive vote payload: epoch 42, slot 1337, snapshot root cafebabe";

/// Best of five batch means, in microseconds
fn time(label: &str, iters: u32, mut f: impl FnMut()) -> f64 {
    for _ in 0..(iters / 4).max(1) {
        f();
    }
    let mut best = f64::MAX;
    for _ in 0..5 {
        let t = Instant::now();
        for _ in 0..iters {
            f();
        }
        best = best.min(t.elapsed().as_nanos() as f64 / iters as f64);
    }
    println!("  {label:<52} {:>9.2} us", best / 1000.0);
    best / 1000.0
}

#[test]
#[ignore = "host timing, run explicitly"]
fn host_speed() {
    use blst::min_pk::{AggregatePublicKey, AggregateSignature, PublicKey, SecretKey, Signature};

    println!("hash_to_G2 (RO, min-pk DST)");
    time("blst blst_hash_to_g2", 2_000, || {
        let mut p = blst::blst_p2::default();
        unsafe {
            blst::blst_hash_to_g2(
                &mut p,
                MESSAGE.as_ptr(),
                MESSAGE.len(),
                DST.as_ptr(),
                DST.len(),
                std::ptr::null(),
                0,
            )
        };
        black_box(p);
    });
    time("blstrs G2Projective::hash_to_curve", 2_000, || {
        black_box(blstrs::G2Projective::hash_to_curve(black_box(MESSAGE), DST, &[]));
    });
    time("solana-bls-signatures HashedMessage::new", 2_000, || {
        black_box(solana_bls_signatures::HashedMessage::new(black_box(MESSAGE)));
    });
    time("zkcrypto bls12_381 hash_to_curve", 200, || {
        use bls12_381::hash_to_curve::{ExpandMsgXmd, HashToCurve};
        black_box(<bls12_381::G2Projective as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::hash_to_curve(
            black_box(MESSAGE),
            DST,
        ));
    });

    println!("this crate: witness generation (host)");
    time("generate_fat (837 B)", 200, || {
        black_box(bls381_hash::witness::g2::generate_fat(black_box(MESSAGE)));
    });
    time("generate (530 B default)", 200, || {
        black_box(bls381_hash::witness::g2::generate(black_box(MESSAGE)));
    });
    time("generate_compact_parity (96 B)", 200, || {
        black_box(bls381_hash::witness::g2::generate_compact_parity(black_box(MESSAGE)));
    });

    println!("field multiply");
    let n = 100_000u64;
    let per = time("this crate mont_mul (ps30, portable)", 5, || {
        black_box(bls381_hash::probe::mont_mul_loop(black_box(n)));
    });
    println!("  {:<52} {:>9.4} us", "  per mont_mul", per / n as f64);
    let per = time("blst blst_fp_mul (asm)", 5, || {
        let mut a = blst::blst_fp::default();
        a.l[0] = 7;
        let b = a;
        for _ in 0..n {
            let x = a;
            unsafe { blst::blst_fp_mul(&mut a, &x, &b) };
        }
        black_box(a);
    });
    println!("  {:<52} {:>9.4} us", "  per blst_fp_mul", per / n as f64);

    println!("blst field ops a witness generator would use (public C API)");
    {
        let mut x = blst::blst_fp2::default();
        x.fp[0].l[0] = 7;
        x.fp[1].l[0] = 3;
        let mut sq = blst::blst_fp2::default();
        unsafe { blst::blst_fp2_sqr(&mut sq, &x) };
        time("blst_fp2_inverse", 2_000, || {
            let mut out = blst::blst_fp2::default();
            unsafe { blst::blst_fp2_inverse(&mut out, black_box(&x)) };
            black_box(out);
        });
        time("blst_fp2_sqrt (a square)", 1_000, || {
            let mut out = blst::blst_fp2::default();
            assert!(unsafe { blst::blst_fp2_sqrt(&mut out, black_box(&sq)) });
            black_box(out);
        });
    }

    println!("min-pk sign / verify");
    let sks: Vec<SecretKey> = (0..20u8).map(|i| SecretKey::key_gen(&[i + 1; 32], &[]).unwrap()).collect();
    let pks: Vec<PublicKey> = sks.iter().map(|s| s.sk_to_pk()).collect();
    let sigs: Vec<Signature> = sks.iter().map(|s| s.sign(MESSAGE, DST, &[])).collect();
    time("blst sign", 1_000, || {
        black_box(sks[0].sign(black_box(MESSAGE), DST, &[]));
    });
    time("blst verify (one signer)", 500, || {
        assert_eq!(
            sigs[0].verify(true, black_box(MESSAGE), DST, &[], &pks[0], true),
            blst::BLST_ERROR::BLST_SUCCESS
        );
    });
    let refs: Vec<&PublicKey> = pks.iter().collect();
    let sig_refs: Vec<&Signature> = sigs.iter().collect();
    let agg_sig = AggregateSignature::aggregate(&sig_refs, false).unwrap().to_signature();
    time("blst fast_aggregate_verify (20 signers)", 500, || {
        assert_eq!(
            agg_sig.fast_aggregate_verify(true, black_box(MESSAGE), DST, &refs),
            blst::BLST_ERROR::BLST_SUCCESS
        );
    });
    let agg_pk = AggregatePublicKey::aggregate(&refs, false).unwrap().to_public_key();
    time("blst verify against a stored 20-key aggregate", 500, || {
        assert_eq!(
            agg_sig.verify(true, black_box(MESSAGE), DST, &[], &agg_pk, false),
            blst::BLST_ERROR::BLST_SUCCESS
        );
    });

    let kp = solana_bls_signatures::Keypair::derive(&[9u8; 32]).unwrap();
    let ssig = kp.sign(MESSAGE);
    time("solana-bls-signatures sign", 1_000, || {
        black_box(kp.sign(black_box(MESSAGE)));
    });
    time("solana-bls-signatures verify (one signer)", 500, || {
        assert!(kp.verify(&ssig, black_box(MESSAGE)).is_ok());
    });
}
