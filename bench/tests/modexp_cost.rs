use mollusk_svm::{program::loader_keys::LOADER_V3, result::InstructionResult, Mollusk};
use solana_instruction::Instruction;
use solana_pubkey::Pubkey;

const ID: Pubkey = Pubkey::new_from_array([7u8; 32]);

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

#[test]
fn modexp_cost_by_exponent() {
    let mollusk = mollusk();
    let p_be = hex::decode(
        "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",
    )
    .unwrap();

    // measure the entrypoint overhead with the cheapest possible call
    for (label, exp_val) in [
        ("exp=1", vec![1u8]),
        ("exp=2 (square)", vec![2u8]),
        ("exp=3", vec![3u8]),
        ("exp=65537", vec![1u8, 0, 1]),
        (
            "exp=(p+1)/4 [sqrt]",
            hex::decode(
                "0680447a8e5ff9a692c6e9ed90d2eb35d91dd2e13ce144afd9cc34a83dac3d8907aaffffac54ffffee7fbfffffffeaab",
            )
            .unwrap(),
        ),
    ] {
        let mut payload = vec![0u8; 48];
        payload[47] = 5; // base = 5
        let mut exp = vec![0u8; 48 - exp_val.len()];
        exp.extend_from_slice(&exp_val);
        payload.extend_from_slice(&exp);
        payload.extend_from_slice(&p_be);
        let r = run(&mollusk, 20, &payload);
        println!(
            "big_mod_exp 48B {label}: rc={:?} {} CU",
            r.return_data.first(),
            r.compute_units_consumed
        );
    }

    // baseline: entrypoint with trivial loop tag for overhead subtraction
    let base = run(&mollusk, 21, &0u64.to_le_bytes());
    println!("entrypoint overhead (tag 21, count=0): {} CU", base.compute_units_consumed);
}
