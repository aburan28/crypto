//! LOCAL-ONLY fast island screener (NOT part of any submission).
//! Faithfully replicates eval_circuit's Fiat-Shamir seed + run_tests, but
//! early-exits on the first failing batch. Prints `CLEAN tof=<avg> q=<n>`
//! when all 9024 shots pass, else `FAIL batch=<b> cf=<..> pg=<..> ag=<..>`.
use alloy_primitives::U256;
use quantum_ecc::circuit::{analyze_ops, BitId, Op, OperationType, QubitId, QubitOrBit, RegisterId};
use quantum_ecc::sim::Simulator;
use quantum_ecc::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
use std::fs;

const OPS_PATH: &str = "ops.bin";
const MAGIC: &[u8; 8] = b"QECCOPS1";
const OP_BYTES: usize = 56;
const NUM_TESTS: usize = 9024;

fn op_kind_from_u32(v: u32) -> Option<OperationType> {
    use OperationType::*;
    Some(match v {
        0 => Neg, 1 => Register, 2 => AppendToRegister, 3 => BitInvert, 4 => BitStore0,
        5 => BitStore1, 6 => X, 7 => Z, 8 => CX, 9 => CZ, 10 => Swap, 11 => R, 12 => Hmr,
        13 => CCX, 14 => CCZ, 15 => PushCondition, 16 => PopCondition, 17 => DebugPrint,
        _ => return None,
    })
}
fn read_u64(b: &[u8], off: usize) -> u64 { u64::from_le_bytes(b[off..off + 8].try_into().unwrap()) }

fn load_ops(path: &str) -> Vec<Op> {
    let bytes = fs::read(path).expect("read ops.bin");
    let n = u64::from_le_bytes(bytes[MAGIC.len()..MAGIC.len() + 8].try_into().unwrap()) as usize;
    let mut ops = Vec::with_capacity(n);
    let mut off = MAGIC.len() + 8;
    for _ in 0..n {
        let kind = op_kind_from_u32(u32::from_le_bytes(bytes[off..off + 4].try_into().unwrap())).unwrap();
        ops.push(Op {
            kind,
            q_control2: QubitId(read_u64(&bytes, off + 8)),
            q_control1: QubitId(read_u64(&bytes, off + 16)),
            q_target: QubitId(read_u64(&bytes, off + 24)),
            c_target: BitId(read_u64(&bytes, off + 32)),
            c_condition: BitId(read_u64(&bytes, off + 40)),
            r_target: RegisterId(read_u64(&bytes, off + 48)),
        });
        off += OP_BYTES;
    }
    ops
}

fn secp256k1() -> WeierstrassEllipticCurve {
    WeierstrassEllipticCurve {
        modulus: U256::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16).unwrap(),
        a: U256::from(0), b: U256::from(7),
        gx: U256::from_str_radix("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16).unwrap(),
        gy: U256::from_str_radix("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16).unwrap(),
        order: U256::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16).unwrap(),
    }
}

fn fiat_shamir_seed(ops: &[Op]) -> sha3::Shake256Reader {
    let mut h = Shake256::default();
    h.update(b"quantum_ecc-fiat-shamir-v2");
    h.update(&(ops.len() as u64).to_le_bytes());
    for op in ops {
        h.update(&[op.kind as u8]);
        h.update(&op.q_control2.0.to_le_bytes());
        h.update(&op.q_control1.0.to_le_bytes());
        h.update(&op.q_target.0.to_le_bytes());
        h.update(&op.c_target.0.to_le_bytes());
        h.update(&op.c_condition.0.to_le_bytes());
        h.update(&op.r_target.0.to_le_bytes());
    }
    h.finalize_xof()
}

fn main() {
    let ops = load_ops(OPS_PATH);
    let (total_qubits, num_bits, _nr, regs) = analyze_ops(ops.iter());
    let curve = secp256k1();
    let mut xof = fiat_shamir_seed(&ops);

    // Cheaply read all xof input bytes first (positions the RNG for Simulator::new
    // exactly as eval does), but DEFER the expensive scalar-mults to per-batch.
    // The eval `continue` filters trigger with prob ~2^-128 (t.x==o.x / identity),
    // so screening assumes none fire (n==NUM_TESTS); every CLEAN hit is reconfirmed
    // by the real eval_circuit before use.
    let mut ks = Vec::with_capacity(NUM_TESTS);
    for _ in 0..NUM_TESTS {
        let mut rb = [[0u8; 32]; 2];
        xof.read(&mut rb[0]);
        xof.read(&mut rb[1]);
        ks.push((U256::from_le_bytes(rb[0]), U256::from_le_bytes(rb[1])));
    }
    let n = NUM_TESTS;

    let mut sim = Simulator::new(total_qubits as usize, num_bits as usize, &mut xof);
    const BATCH: usize = 64;
    let num_batches = (n + BATCH - 1) / BATCH;
    for batch in 0..num_batches {
        let bs = BATCH.min(n - batch * BATCH);
        let cond_mask: u64 = if bs == 64 { u64::MAX } else { (1u64 << bs) - 1 };
        // Lazily materialize this batch's inputs.
        let mut bt = [(U256::ZERO, U256::ZERO); BATCH];
        let mut bo = [(U256::ZERO, U256::ZERO); BATCH];
        let mut be = [(U256::ZERO, U256::ZERO); BATCH];
        for shot in 0..bs {
            let (k1, k2) = ks[batch * BATCH + shot];
            let t = curve.mul(curve.gx, curve.gy, k1);
            let o = curve.mul(curve.gx, curve.gy, k2);
            bt[shot] = t; bo[shot] = o;
            be[shot] = curve.add(t.0, t.1, o.0, o.1);
        }
        sim.clear_for_shot();
        for shot in 0..bs {
            sim.set_register(&regs[0], bt[shot].0, shot);
            sim.set_register(&regs[1], bt[shot].1, shot);
            sim.set_register(&regs[2], bo[shot].0, shot);
            sim.set_register(&regs[3], bo[shot].1, shot);
        }
        sim.apply_iter(ops.iter());
        let mut cf = 0usize;
        for shot in 0..bs {
            let gx = sim.get_register(&regs[0], shot);
            let gy = sim.get_register(&regs[1], shot);
            if gx != be[shot].0 || gy != be[shot].1 { cf += 1; }
        }
        let pg = if (sim.phase & cond_mask) != 0 { 1 } else { 0 };
        for register in &regs {
            for qb in register {
                if let QubitOrBit::Qubit(q) = *qb { *sim.qubit_mut(q) = 0; }
            }
        }
        let mut ag = 0usize;
        for q in 0..total_qubits {
            if (sim.qubit(QubitId(q)) & cond_mask) != 0 { ag = 1; break; }
        }
        if cf != 0 || pg != 0 || ag != 0 {
            println!("FAIL batch={batch} cf={cf} pg={pg} ag={ag}");
            return;
        }
    }
    let avg_tof = sim.stats.toffoli_gates as f64 / n.max(1) as f64;
    println!("CLEAN tof={:.0} q={}", avg_tof, total_qubits);
}
