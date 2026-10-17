//! **Binary elliptic curves**: `y² + xy = x³ + a·x² + b` over
//! `F_{2^m}`.
//!
//! This is the standard Weierstrass form for characteristic-2
//! fields used by NIST B-163, B-233, B-283, B-409, B-571 (SEC 2,
//! FIPS 186-4 Appendix D.1.3).
//!
//! ## Group law (Putranto et al. §3)
//!
//! For affine points `P₁ = (x₁, y₁)`, `P₂ = (x₂, y₂)`:
//!
//! **Negation**: `−P₁ = (x₁, y₁ + x₁)`.
//!
//! **Point Addition (P₁ ≠ ±P₂)**:
//! ```text
//! λ  = (y₁ + y₂) / (x₁ + x₂)
//! x₃ = λ² + λ + x₁ + x₂ + a
//! y₃ = λ·(x₁ + x₃) + x₃ + y₁
//! ```
//!
//! **Point Doubling (P₁ ≠ −P₁)**:
//! ```text
//! λ  = x₁ + y₁/x₁
//! x₃ = λ² + λ + a
//! y₃ = x₁² + (λ + 1)·x₃
//! ```
//!
//! ## ECPM
//!
//! Scalar multiplication `[k]·P` via left-to-right double-and-add.
//! In the paper, this is the operation whose **quantum-circuit
//! resource estimation** is the main contribution; we provide the
//! **classical reference** here.

use super::f2m::{F2mElement, IrreduciblePoly};
use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};
use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

/// A binary elliptic curve over `F_{2^m}`.
#[derive(Clone, Debug)]
pub struct BinaryCurve {
    pub m: u32,
    pub irreducible: IrreduciblePoly,
    /// Coefficient `a` of `y² + xy = x³ + a·x² + b`.
    pub a: F2mElement,
    /// Coefficient `b` of the curve equation.
    pub b: F2mElement,
    /// A base point of large prime order (the generator).
    pub generator: BinaryPoint,
    /// Order of the generator (a prime number, as `BigUint`).
    pub order: BigUint,
    /// Cofactor of the curve (`h = #E / n`).
    pub cofactor: BigUint,
}

/// An affine point on a binary elliptic curve, or the point at
/// infinity (`O`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BinaryPoint {
    Infinity,
    Affine { x: F2mElement, y: F2mElement },
}

impl BinaryCurve {
    /// Certicom binary challenge parameters and the published target Q.
    /// Source: https://www.certicom.com/en/curves-list (retrieved 2026-10-16;
    /// downloaded HTML SHA-256 34ca84c9ef1021ad2f7a46726cdfd68b209a362fbd0e82a2e55903151d2ef857).
    /// The source repeats `a` for ECC2K-163; its Koblitz coefficient b is 1.
    pub fn certicom_challenge(name: &str) -> Option<(Self, BinaryPoint)> {
        let (m, terms, a, b, px, py, qx, qy, n, h): (u32, &[u32], &str, &str, &str, &str, &str, &str, &str, &str) =
            match name.to_ascii_lowercase().as_str() {
                "ecc2k-130" => (131, &[0, 1, 2, 13],
                    "0000000000000000000000000000000000", "0000000000000000000000000000000001", "051C99BFA6F18DE467C80C23B98C7994AA", "042EA2D112ECEC71FCF7E000D7EFC978BD",
                    "06C997F3E7F2C66A4A5D2FDA13756A37B1", "04A38D11829D32D347BD0C0F584D546E9A", "200000000000000004D4FDD5703A3F269", "04"),
                "ecc2-131" => (131, &[0, 1, 2, 13],
                    "07EBCB7EECC296A1C4A1A14F2C9E44352E", "00610B0A57C73649AD0093BDD622A61D81", "00439CBC8DC73AA981030D5BC57B331663", "014904C07D4F25A16C2DE036D60B762BD4",
                    "0602339C5DB0E9C694AC8908528C51C440", "04F7B99169FA1A0F2737813742B1588CB8", "0400000000000000026ABB991FE311FE83", "02"),
                "ecc2k-163" => (163, &[0, 1, 2, 8],
                    "000000000000000000000000000000000000000001", "1", "02091945E42080CD9CBCF14A7107D8BC55CDD65EA9", "063315693833774294A39CF6F8C175D02B8E6A5587",
                    "007530EE864EDCF4A31C85AA17C197FFF5CAFECAE1", "075DB1E80D7C4A92C7BBB79EAE3EC545F8A31CFA6B", "04000000000000000000020108A2E0CC0D99F8A5EF", "02"),
                "ecc2-163" => (163, &[0, 1, 2, 8],
                    "025C4BEAC8074B8C2D9DF63AF91263EB8229B3C967", "00C9517D06D5240D3CFF38C74B20B6CD4D6F9DD4D9", "023A2E99904996E8679B50FF1E49ADD8BD2388F387", "05FCBFE4098477C9D187EA1CF615C7E91529E73BA2",
                    "0438D8B3821C8E9264637F2FC74F8007B21210F0F2", "073FCEA8D5E247CE367368F006EBD5B32FDF4286D2", "0400000000000000000001E60FC8821CC74DAEAFC1", "02"),
                "ecc2-191" => (191, &[0, 9],
                    "2011D237B30025CFD61F8D6F66741CA6CEDD313F9BD4431F", "2E45EF571F00786F67B0081B9495A3D95462F5DE0AA185EC", "52F3E1ECFCF3B0A6E0945D6B92A8FFC42C280545C0F8321B", "00F19F7A172A1F7191614DB1DAE7CABC01E85806D42E2627",
                    "7911BC300972E2D3385D6D528440710DD73F281472B961E3", "492A2DED9E740B3603F39DABAA0070F99674C57F55AD4F1D", "40000000000000000000000004A20E90C39067C893BBB9A5", "02"),
                "ecc2k-238" => (239, &[0, 36],
                    "000000000000000000000000000000000000000000000000000000000000", "000000000000000000000000000000000000000000000000000000000001", "78CFAE9860F96A9D48261B7A9CF3C5B94190888B184BC3F83B416B56992E", "2DEBA096F37DEDBE7D13C09F6C03680B2E5371A08FEC01E37775E31CC7C9",
                    "721DD31C6766532EC44C73BCBF54B2969400FA6C4CC1556E89AC1FF8622A", "347D6F86DDDEDD7F138FC33BEEAB4FFA029A2EB2F011A5BC8C64A773B519", "2000000000000000000000000000005A79FEC67CB6E91F1C1DA800E478A5", "04"),
                "ecc2-238" => (239, &[0, 36],
                    "23AACC83DE3549BEAE1B601EDA88BF45517D5ADFED437F71784B5DD911DF", "790408F2EEDAF392B012EDEFB3392F30F4327C0CA3F31FC383C422AA8C16", "4BC630D209D98B810FA42A95FDB77F87DBE68CC076110AADE0B2AFC47B2B", "194C8DA12F2763E3E9AA84F4F17D56FF10FCF20F300FD1278FFB976A91C2",
                    "26966E5DDE492AF6CEDAF8BB8BCADA02BFEF2944815F142F57450CA145BF", "4F2AC9ED29887CB53ACBB29AD031B4A24C2CFA8A32D5B164E89FA283966B", "2000000000000000000000000000000F4D42FFE1492A4993F1CAD666E447", "04"),
                "ecc2k-358" => (359, &[0, 68],
                    "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", "7035ADDF5EAAB0CFBA9588332603F870D11AA41A6B57ADA3A1CA18F9C7169185A1F191E96293EE6E041C8E918C", "7122A49D2EF45E68D11C0C0ED86B8AE422F54E49BDA0AC15177F3824D35809AD4FAEF7B7AB4810E1DBEC5ABFF0",
                    "0FE3828D45A36EC0A9D04BB8C03586A41ED098F4BC55718A7604A005DB050E58F6DC2A451CB9D83448D7B1DDA6", "5322AFB6DA476F2B9C3ABF31C22978C318AE358F6B3373CEAC599189A35DD6C121C23170D710563F436CAAFDC8", "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF6F61406970B03AE00FB512716A4151CA0EAC9D1B0E107", "02"),
                "ecc2-353" => (359, &[0, 68],
                    "6962E67E43FB6D66C46C7FDDF5A22E2C8857537C428FD691A26D2F9BAEB871C276DFA9928F01634B1459CF431E", "2472E2D0197C49363F1FE7F5B6DB075D52B6947D135D8CA445805D39BC345626089687742B6329E70680231988", "0EF6A87D27D25EA944A629664A9D34664C6A79DD86360D8AB57FAE11DF30CB6E78C7E01BA3ADC081FA76048092", "1DA4ED719174E356FC0F968E80D2974965BB90805DBB1D8FE57333DFA50D163B88465BB5FB55046F949AE27F51",
                    "4BD3C6E97BFF49C1B3FD667684834FBE8D17E22F94F795480BD5166B04A267F94033D61FBA43FDF62758982462", "605F4BC8221F266F724DEF750069B240A67AC4ADAC15386354D7E1700C5F5F3FB6769F26130D73A8D721A019E4", "01AF286BCA1AF286BCA1AF286BCA1AF286BCA1AF286BC9FB8F6B85C556892C20A7EB964FE7719E74F490758D3B", "4C"),
                _ => return None,
            };
        let irr = IrreduciblePoly {
            degree: m,
            low_terms: terms.to_vec(),
        };
        let point = |x: &str, y: &str| BinaryPoint::Affine {
            x: F2mElement::from_hex(x, m),
            y: F2mElement::from_hex(y, m),
        };
        let curve = Self {
            m,
            irreducible: irr,
            a: F2mElement::from_hex(a, m),
            b: F2mElement::from_hex(b, m),
            generator: point(px, py),
            order: BigUint::parse_bytes(n.as_bytes(), 16).expect("Certicom subgroup order"),
            cofactor: BigUint::parse_bytes(h.as_bytes(), 16).expect("Certicom cofactor"),
        };
        Some((curve, point(qx, qy)))
    }

    /// **Low-security historical binary curves** from SECG and IKE/Oakley.
    ///
    /// Includes the 56-bit, 64-bit, and pre-NIST IKE EC2N groups used in
    /// WTLS/WAP and IPsec-era interoperability.  These are retained for
    /// catalog/audit coverage and attack demos only.
    pub fn legacy_low_security_family() -> Vec<(&'static str, Self)> {
        vec![
            ("sect113r1", Self::sect113r1()),
            ("sect113r2", Self::sect113r2()),
            ("sect131r1", Self::sect131r1()),
            ("sect131r2", Self::sect131r2()),
            ("ike-oakley-group3", Self::ike_oakley_group3()),
            ("sect163k1", Self::sect163k1()),
            ("sect163r1", Self::sect163r1()),
            ("sect163r2", Self::sect163r2()),
        ]
    }

    /// **sect113 family** from SEC 2 / SECG, WTLS/WAP-era 56-bit curves.
    pub fn sect113_family() -> Vec<(&'static str, Self)> {
        vec![
            ("sect113r1", Self::sect113r1()),
            ("sect113r2", Self::sect113r2()),
        ]
    }

    /// **sect131 family** from SEC 2 / SECG, WTLS/WAP-era 64-bit curves.
    pub fn sect131_family() -> Vec<(&'static str, Self)> {
        vec![
            ("sect131r1", Self::sect131r1()),
            ("sect131r2", Self::sect131r2()),
        ]
    }

    /// **sect163 family** from SEC 2 / NIST, historically deployed at
    /// roughly 80-bit generic ECDLP security.
    ///
    /// These binary-field curves are retained for legacy audits,
    /// cryptanalytic demonstrations, and regression tests. They are not
    /// suitable for new production systems.
    pub fn sect163_family() -> Vec<(&'static str, Self)> {
        vec![
            ("sect163k1", Self::sect163k1()),
            ("sect163r1", Self::sect163r1()),
            ("sect163r2", Self::sect163r2()),
        ]
    }

    /// **sect113r1** (SEC 2 v1 §3.2.1).
    ///
    /// Binary random curve over `F_{2^113}` with reduction polynomial
    /// `z¹¹³ + z⁹ + 1`; roughly 56-bit generic ECDLP security.
    pub fn sect113r1() -> Self {
        let m = 113;
        let irr = IrreduciblePoly::deg_113();
        let a = F2mElement::from_hex("003088250CA6E7C7FE649CE85820F7", m);
        let b = F2mElement::from_hex("00E8BEE4D3E2260744188BE0E9C723", m);
        let gx = F2mElement::from_hex("009D73616F35F4AB1407D73562C10F", m);
        let gy = F2mElement::from_hex("00A52830277958EE84D1315ED31886", m);
        let order = BigUint::parse_bytes(b"0100000000000000D9CCEC8A39E56F", 16).unwrap();
        Self {
            m,
            irreducible: irr,
            a,
            b,
            generator: BinaryPoint::Affine { x: gx, y: gy },
            order,
            cofactor: BigUint::from(2u32),
        }
    }

    /// **sect113r2** (SEC 2 v1 §3.2.2).
    ///
    /// Binary random curve over `F_{2^113}` with reduction polynomial
    /// `z¹¹³ + z⁹ + 1`; roughly 56-bit generic ECDLP security.
    pub fn sect113r2() -> Self {
        let m = 113;
        let irr = IrreduciblePoly::deg_113();
        let a = F2mElement::from_hex("00689918DBEC7E5A0DD6DFC0AA55C7", m);
        let b = F2mElement::from_hex("0095E9A9EC9B297BD4BF36E059184F", m);
        let gx = F2mElement::from_hex("01A57A6A7B26CA5EF52FCDB8164797", m);
        let gy = F2mElement::from_hex("00B3ADC94ED1FE674C06E695BABA1D", m);
        let order = BigUint::parse_bytes(b"010000000000000108789B2496AF93", 16).unwrap();
        Self {
            m,
            irreducible: irr,
            a,
            b,
            generator: BinaryPoint::Affine { x: gx, y: gy },
            order,
            cofactor: BigUint::from(2u32),
        }
    }

    /// **sect131r1** (SEC 2 v1 §3.3.1).
    ///
    /// Binary random curve over `F_{2^131}` with reduction polynomial
    /// `z¹³¹ + z⁸ + z³ + z² + 1`; roughly 64-bit generic ECDLP security.
    pub fn sect131r1() -> Self {
        let m = 131;
        let irr = IrreduciblePoly::deg_131();
        let a = F2mElement::from_hex("07A11B09A76B562144418FF3FF8C2570B8", m);
        let b = F2mElement::from_hex("0217C05610884B63B9C6C7291678F9D341", m);
        let gx = F2mElement::from_hex("0081BAF91FDF9833C40F9C181343638399", m);
        let gy = F2mElement::from_hex("078C6E7EA38C001F73C8134B1B4EF9E150", m);
        let order = BigUint::parse_bytes(b"0400000000000000023123953A9464B54D", 16).unwrap();
        Self {
            m,
            irreducible: irr,
            a,
            b,
            generator: BinaryPoint::Affine { x: gx, y: gy },
            order,
            cofactor: BigUint::from(2u32),
        }
    }

    /// **sect131r2** (SEC 2 v1 §3.3.2).
    ///
    /// Binary random curve over `F_{2^131}` with reduction polynomial
    /// `z¹³¹ + z⁸ + z³ + z² + 1`; roughly 64-bit generic ECDLP security.
    pub fn sect131r2() -> Self {
        let m = 131;
        let irr = IrreduciblePoly::deg_131();
        let a = F2mElement::from_hex("03E5A88919D7CAFCBF415F07C2176573B2", m);
        let b = F2mElement::from_hex("04B8266A46C55657AC734CE38F018F2192", m);
        let gx = F2mElement::from_hex("0356DCD8F2F95031AD652D23951BB366A8", m);
        let gy = F2mElement::from_hex("0648F06D867940A5366D9E265DE9EB240F", m);
        let order = BigUint::parse_bytes(b"0400000000000000016954A233049BA98F", 16).unwrap();
        Self {
            m,
            irreducible: irr,
            a,
            b,
            generator: BinaryPoint::Affine { x: gx, y: gy },
            order,
            cofactor: BigUint::from(2u32),
        }
    }

    /// **IKE/Oakley Group 3** (RFC 2409 §6.3), EC2N over `GF(2^155)`.
    ///
    /// RFC 2409 publishes the generator by x-coordinate (`0x7b`) and
    /// permits either of the two y-coordinates satisfying the curve
    /// equation.  We use `y = 0x1c8`, one valid solution under the
    /// polynomial basis `z¹⁵⁵ + z⁶² + 1`, so the group is usable by the
    /// existing affine arithmetic while preserving the RFC generator x.
    pub fn ike_oakley_group3() -> Self {
        let m = 155;
        let irr = IrreduciblePoly::deg_155_oakley_group3();
        let a = F2mElement::zero(m);
        let b = F2mElement::from_hex("07338F", m);
        let gx = F2mElement::from_hex("7B", m);
        let gy = F2mElement::from_hex("01C8", m);
        let order = BigUint::parse_bytes(b"0800000000000000000057DB5698537193AEF944", 16).unwrap();
        Self {
            m,
            irreducible: irr,
            a,
            b,
            generator: BinaryPoint::Affine { x: gx, y: gy },
            order,
            // RFC 2409 gives a group order, not a separate prime-subgroup
            // order and cofactor pair.
            cofactor: BigUint::one(),
        }
    }

    /// **NIST sect163k1** (Koblitz curve over `F_{2^163}`).
    ///
    /// Curve equation: `y² + xy = x³ + x² + 1` (a = 1, b = 1).
    /// Reduction polynomial: `z¹⁶³ + z⁷ + z⁶ + z³ + 1`.
    /// Per FIPS 186-4 Appendix D.1.3 (and SEC 2 v2.0 §3.1).
    pub fn sect163k1() -> Self {
        let m = 163;
        let irr = IrreduciblePoly::deg_163();
        let a = F2mElement::one(m); // a = 1
        let b = F2mElement::one(m); // b = 1
                                    // Generator coordinates from SEC 2.
        let gx = F2mElement::from_hex("02fe13c0537bbc11acaa07d793de4e6d5e5c94eee8", m);
        let gy = F2mElement::from_hex("0289070fb05d38ff58321f2e800536d538ccdaa3d9", m);
        // Order n (SEC 2): 4000000000000000000020108A2E0CC0D99F8A5EF
        let order = BigUint::parse_bytes(b"4000000000000000000020108A2E0CC0D99F8A5EF", 16).unwrap();
        Self {
            m,
            irreducible: irr,
            a,
            b,
            generator: BinaryPoint::Affine { x: gx, y: gy },
            order,
            cofactor: BigUint::from(2u32),
        }
    }

    /// **sect163r1** (SEC 2 v2 §3.1.2).
    ///
    /// Older random binary curve over `F_{2^163}`. Legacy/insecure by
    /// modern standards; included so binary-curve attack tooling can
    /// cover historical Internet-standard parameters.
    pub fn sect163r1() -> Self {
        let m = 163;
        let irr = IrreduciblePoly::deg_163();
        let a = F2mElement::from_hex("07B6882CAAEFA84F9554FF8428BD88E246D2782AE2", m);
        let b = F2mElement::from_hex("0713612DCDDCB40AAB946BDA29CA91F73AF958AFD9", m);
        let gx = F2mElement::from_hex("0369979697AB43897789566789567F787A7876A654", m);
        let gy = F2mElement::from_hex("00435EDB42EFAFB2989D51FEFCE3C80988F41FF883", m);
        let order =
            BigUint::parse_bytes(b"03FFFFFFFFFFFFFFFFFFFF48AAB689C29CA710279B", 16).unwrap();
        Self {
            m,
            irreducible: irr,
            a,
            b,
            generator: BinaryPoint::Affine { x: gx, y: gy },
            order,
            cofactor: BigUint::from(2u32),
        }
    }

    /// **NIST B-163 / sect163r2** (FIPS 186-4 §D.1.3, SEC 2 v2 §3.1.3).
    ///
    /// Random binary curve over `F_{2^163}`, roughly 80-bit generic ECDLP
    /// security. Present only for historical/audit and attack-demo use.
    pub fn sect163r2() -> Self {
        let m = 163;
        let irr = IrreduciblePoly::deg_163();
        let a = F2mElement::one(m);
        let b = F2mElement::from_hex("020A601907B8C953CA1481EB10512F78744A3205FD", m);
        let gx = F2mElement::from_hex("03F0EBA16286A2D57EA0991168D4994637E8343E36", m);
        let gy = F2mElement::from_hex("00D51FBC6C71A0094FA2CDD545B11C5C0C797324F1", m);
        let order =
            BigUint::parse_bytes(b"040000000000000000000292FE77E70C12A4234C33", 16).unwrap();
        Self {
            m,
            irreducible: irr,
            a,
            b,
            generator: BinaryPoint::Affine { x: gx, y: gy },
            order,
            cofactor: BigUint::from(2u32),
        }
    }

    /// Build a curve from hex parameter strings.  All of `a, b, gx, gy` are
    /// parsed as `F_{2^m}` elements and `n` as the subgroup order; this keeps
    /// the many standardized-curve constructors below terse and uniform.
    fn from_hex_parts(
        m: u32,
        irr: IrreduciblePoly,
        a: &str,
        b: &str,
        gx: &str,
        gy: &str,
        n: &str,
        h: u32,
    ) -> Self {
        Self {
            m,
            irreducible: irr,
            a: F2mElement::from_hex(a, m),
            b: F2mElement::from_hex(b, m),
            generator: BinaryPoint::Affine {
                x: F2mElement::from_hex(gx, m),
                y: F2mElement::from_hex(gy, m),
            },
            order: BigUint::parse_bytes(n.as_bytes(), 16).expect("valid order hex"),
            cofactor: BigUint::from(h),
        }
    }

    /// **sect193r1** (SEC 2 v2 §3.2.1). Random binary curve over `F_{2^193}`.
    pub fn sect193r1() -> Self {
        Self::from_hex_parts(
            193,
            IrreduciblePoly::deg_193(),
            "0017858feb7a98975169e171f77b4087de098ac8a911df7b01",
            "00fdfb49bfe6c3a89facadaa7a1e5bbc7cc1c2e5d831478814",
            "01f481bc5f0ff84a74ad6cdf6fdef4bf6179625372d8c0c5e1",
            "0025e399f2903712ccf3ea9e3a1ad17fb0b3201b6af7ce1b05",
            "01000000000000000000000000c7f34a778f443acc920eba49",
            2,
        )
    }

    /// **sect193r2** (SEC 2 v2 §3.2.2). Random binary curve over `F_{2^193}`.
    pub fn sect193r2() -> Self {
        Self::from_hex_parts(
            193,
            IrreduciblePoly::deg_193(),
            "0163f35a5137c2ce3ea6ed8667190b0bc43ecd69977702709b",
            "00c9bb9e8927d4d64c377e2ab2856a5b16e3efb7f61d4316ae",
            "00d9b67d192e0367c803f39e1a7e82ca14a651350aae617e8f",
            "01ce94335607c304ac29e7defbd9ca01f596f927224cdecf6c",
            "010000000000000000000000015aab561b005413ccd4ee99d5",
            2,
        )
    }

    /// **NIST K-233 / sect233k1** (FIPS 186-4, SEC 2). Koblitz, `a=0, b=1`.
    pub fn sect233k1() -> Self {
        Self::from_hex_parts(
            233,
            IrreduciblePoly::deg_233(),
            "0",
            "1",
            "017232ba853a7e731af129f22ff4149563a419c26bf50a4c9d6eefad6126",
            "01db537dece819b7f70f555a67c427a8cd9bf18aeb9b56e0c11056fae6a3",
            "8000000000000000000000000000069d5bb915bcd46efb1ad5f173abdf",
            4,
        )
    }

    /// **NIST B-233 / sect233r1** (FIPS 186-4, SEC 2). Random binary curve.
    pub fn sect233r1() -> Self {
        Self::from_hex_parts(
            233,
            IrreduciblePoly::deg_233(),
            "1",
            "0066647ede6c332c7f8c0923bb58213b333b20e9ce4281fe115f7d8f90ad",
            "00fac9dfcbac8313bb2139f1bb755fef65bc391f8b36f8f8eb7371fd558b",
            "01006a08a41903350678e58528bebf8a0beff867a7ca36716f7e01f81052",
            "1000000000000000000000000000013e974e72f8a6922031d2603cfe0d7",
            2,
        )
    }

    /// **sect239k1** (SEC 2 v2 §3.4.1). Koblitz over `F_{2^239}`, `a=0, b=1`.
    pub fn sect239k1() -> Self {
        Self::from_hex_parts(
            239,
            IrreduciblePoly::deg_239(),
            "0",
            "1",
            "29a0b6a887a983e9730988a68727a8b2d126c44cc2cc7b2a6555193035dc",
            "76310804f12e549bdb011c103089e73510acb275fc312a5dc6b76553f0ca",
            "2000000000000000000000000000005a79fec67cb6e91f1c1da800e478a5",
            4,
        )
    }

    /// **NIST K-283 / sect283k1** (FIPS 186-4, SEC 2). Koblitz, `a=0, b=1`.
    pub fn sect283k1() -> Self {
        Self::from_hex_parts(
            283,
            IrreduciblePoly::deg_283(),
            "0",
            "1",
            "503213f78ca44883f1a3b8162f188e553cd265f23c1567a16876913b0c2ac2458492836",
            "1ccda380f1c9e318d90f95d07e5426fe87e45c0e8184698e45962364e34116177dd2259",
            "1ffffffffffffffffffffffffffffffffffe9ae2ed07577265dff7f94451e061e163c61",
            4,
        )
    }

    /// **NIST B-283 / sect283r1** (FIPS 186-4, SEC 2). Random binary curve.
    pub fn sect283r1() -> Self {
        Self::from_hex_parts(
            283,
            IrreduciblePoly::deg_283(),
            "1",
            "27b680ac8b8596da5a4af8a19a0303fca97fd7645309fa2a581485af6263e313b79a2f5",
            "5f939258db7dd90e1934f8c70b0dfec2eed25b8557eac9c80e2e198f8cdbecd86b12053",
            "3676854fe24141cb98fe6d4b20d02b4516ff702350eddb0826779c813f0df45be8112f4",
            "3ffffffffffffffffffffffffffffffffffef90399660fc938a90165b042a7cefadb307",
            2,
        )
    }

    /// **NIST K-409 / sect409k1** (FIPS 186-4, SEC 2). Koblitz, `a=0, b=1`.
    pub fn sect409k1() -> Self {
        Self::from_hex_parts(
            409,
            IrreduciblePoly::deg_409(),
            "0",
            "1",
            "060f05f658f49c1ad3ab1890f7184210efd0987e307c84c27accfb8f9f67cc2c460189eb5aaaa62ee222eb1b35540cfe9023746",
            "1e369050b7c4e42acba1dacbf04299c3460782f918ea427e6325165e9ea10e3da5f6c42e9c55215aa9ca27a5863ec48d8e0286b",
            "7ffffffffffffffffffffffffffffffffffffffffffffffffffe5f83b2d4ea20400ec4557d5ed3e3e7ca5b4b5c83b8e01e5fcf",
            4,
        )
    }

    /// **NIST B-409 / sect409r1** (FIPS 186-4, SEC 2). Random binary curve.
    pub fn sect409r1() -> Self {
        Self::from_hex_parts(
            409,
            IrreduciblePoly::deg_409(),
            "1",
            "021a5c2c8ee9feb5c4b9a753b7b476b7fd6422ef1f3dd674761fa99d6ac27c8a9a197b272822f6cd57a55aa4f50ae317b13545f",
            "15d4860d088ddb3496b0c6064756260441cde4af1771d4db01ffe5b34e59703dc255a868a1180515603aeab60794e54bb7996a7",
            "061b1cfab6be5f32bbfa78324ed106a7636b9c5a7bd198d0158aa4f5488d08f38514f1fdf4b4f40d2181b3681c364ba0273c706",
            "10000000000000000000000000000000000000000000000000001e2aad6a612f33307be5fa47c3c9e052f838164cd37d9a21173",
            2,
        )
    }

    /// **NIST K-571 / sect571k1** (FIPS 186-4, SEC 2). Koblitz, `a=0, b=1`.
    pub fn sect571k1() -> Self {
        Self::from_hex_parts(
            571,
            IrreduciblePoly::deg_571(),
            "0",
            "1",
            "26eb7a859923fbc82189631f8103fe4ac9ca2970012d5d46024804801841ca44370958493b205e647da304db4ceb08cbbd1ba39494776fb988b47174dca88c7e2945283a01c8972",
            "349dc807f4fbf374f4aeade3bca95314dd58cec9f307a54ffc61efc006d8a2c9d4979c0ac44aea74fbebbb9f772aedcb620b01a7ba7af1b320430c8591984f601cd4c143ef1c7a3",
            "20000000000000000000000000000000000000000000000000000000000000000000000131850e1f19a63e4b391a8db917f4138b630d84be5d639381e91deb45cfe778f637c1001",
            4,
        )
    }

    /// **NIST B-571 / sect571r1** (FIPS 186-4, SEC 2). Random binary curve.
    pub fn sect571r1() -> Self {
        Self::from_hex_parts(
            571,
            IrreduciblePoly::deg_571(),
            "1",
            "2f40e7e2221f295de297117b7f3d62f5c6a97ffcb8ceff1cd6ba8ce4a9a18ad84ffabbd8efa59332be7ad6756a66e294afd185a78ff12aa520e4de739baca0c7ffeff7f2955727a",
            "303001d34b856296c16c0d40d3cd7750a93d1d2955fa80aa5f40fc8db7b2abdbde53950f4c0d293cdd711a35b67fb1499ae60038614f1394abfa3b4c850d927e1e7769c8eec2d19",
            "37bf27342da639b6dccfffeb73d69d78c6c27a6009cbbca1980f8533921e8a684423e43bab08a576291af8f461bb2a8b3531d2f0485c19b16e2f1516e23dd3c1a4827af1b8ac15b",
            "3ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffe661ce18ff55987308059b186823851ec7dd9ca1161de93d5174d66e8382e9bb2fe84e47",
            2,
        )
    }

    /// Tiny test curve over `F_{2^8}` with `a = 1, b = 1`.  Used
    /// for unit-testing without depending on the full sect163k1
    /// parameters.
    pub fn test_curve_f256() -> Self {
        let m = 8;
        let irr = IrreduciblePoly::deg_8();
        let a = F2mElement::one(m);
        let b = F2mElement::one(m);
        // Pick a generator by trial (we'll find one valid point).
        let _gx = F2mElement::from_bit_positions(&[1], m); // x = z
                                                           // Find y satisfying y² + zy = z³ + z² + 1 ...
                                                           // Just brute-force search for a valid point.
        let (gen, ord) = find_any_generator(m, &irr, &a, &b);
        Self {
            m,
            irreducible: irr,
            a,
            b,
            generator: gen,
            order: ord,
            cofactor: BigUint::one(),
            // Note: cofactor unknown precisely without full point
            // counting; for the test curve it's fine.
        }
    }

    /// Verify a point satisfies the curve equation.
    pub fn is_on_curve(&self, p: &BinaryPoint) -> bool {
        match p {
            BinaryPoint::Infinity => true,
            BinaryPoint::Affine { x, y } => {
                // y² + xy = x³ + a·x² + b
                let y_sq = y.square(&self.irreducible);
                let xy = x.mul(y, &self.irreducible);
                let lhs = y_sq.add(&xy);
                let x_sq = x.square(&self.irreducible);
                let x_cu = x_sq.mul(x, &self.irreducible);
                let ax_sq = self.a.mul(&x_sq, &self.irreducible);
                let rhs = x_cu.add(&ax_sq).add(&self.b);
                lhs == rhs
            }
        }
    }

    /// Validate a non-identity point as a member of the curve's declared
    /// prime-order subgroup.
    ///
    /// [`scalar_mul`] is deliberately a low-level arithmetic primitive: it
    /// does not call this method and therefore accepts formula-compatible
    /// points from other equations. Protocols that process attacker-selected
    /// points must perform this check (or an equivalent, proved validation)
    /// before multiplying by a secret scalar.
    ///
    /// On curves where [`Self::has_fast_4e_subgroup_check`] holds
    /// (`y² + xy = x³ + 1`, odd `m`, cofactor 4: sect233k1, sect239k1,
    /// sect283k1, sect409k1, sect571k1) the `[n]P = O` check is replaced by
    /// the closed-form [`Self::is_in_4e`] test, which costs one field
    /// multiplication, a half-trace and three traces instead of a full
    /// scalar multiplication.  Debug builds still run `[n]P = O` and assert
    /// that both agree.  Every other curve takes the `[n]P = O` path.
    ///
    /// Not constant-time on either path: the inputs are public points.
    pub fn is_valid_public_point(&self, p: &BinaryPoint) -> bool {
        match p {
            BinaryPoint::Infinity => false,
            BinaryPoint::Affine { x, y } if x.m_value() != self.m || y.m_value() != self.m => false,
            BinaryPoint::Affine { .. } => {
                if self.has_fast_4e_subgroup_check() {
                    let fast = self.is_in_4e(p) == Some(true);
                    debug_assert_eq!(
                        fast,
                        self.is_on_curve(p)
                            && matches!(scalar_mul(self, p, &self.order), BinaryPoint::Infinity),
                        "closed-form 4E test disagrees with [n]P = O"
                    );
                    return fast;
                }
                self.is_on_curve(p)
                    && matches!(scalar_mul(self, p, &self.order), BinaryPoint::Infinity)
            }
        }
    }

    /// Whether this is the Koblitz curve `y² + xy = x³ + 1` over
    /// `F_{2^m}` with `m` odd, the shape on which [`Self::is_in_4e`]
    /// applies.
    fn is_odd_degree_koblitz_a0(&self) -> bool {
        self.m % 2 == 1
            && self.irreducible.degree == self.m
            && self.a.m_value() == self.m
            && self.a.is_zero()
            && self.b == F2mElement::one(self.m)
    }

    /// Whether [`Self::is_valid_public_point`] decides subgroup membership
    /// with the closed-form [`Self::is_in_4e`] test instead of `[n]P = O`.
    ///
    /// True exactly when the curve is `y² + xy = x³ + 1` over `F_{2^m}`
    /// with `m` odd, the declared cofactor is 4, and the declared subgroup
    /// order equals the actual curve order divided by 4.  For odd `m` the
    /// 2-Sylow subgroup of `E(F_{2^m})` is exactly `Z/4` (`(1, 0)` has
    /// order 4 and `Tr(1) = 1`, so it has no half), so with
    /// `#E = 4·n`, `n` odd, `4E` is the order-`n` subgroup.  The actual
    /// order is computed from the Koblitz Frobenius trace recurrence and
    /// cached per field degree.  Inconsistent declarations keep the
    /// `[n]P = O` path.
    pub fn has_fast_4e_subgroup_check(&self) -> bool {
        self.is_odd_degree_koblitz_a0()
            && self.cofactor == BigUint::from(4u32)
            && self.order == *koblitz_a0_4e_order(self.m)
    }

    /// **Closed-form `4E` membership** on `y² + xy = x³ + 1` over
    /// `F_{2^m}`, `m` odd:
    ///
    /// ```text
    /// P ∈ 4E  ⇔  P = O,  or  [Tr(x) = 0  and  Tr(y) + Tr(x·H(x)) = 0]
    /// ```
    ///
    /// where `Tr` is the absolute trace and `H` the half-trace
    /// ([`F2mElement::trace`], [`F2mElement::half_trace`], both table
    /// driven).  Returns `None` when the curve is not of that shape
    /// (`a ≠ 0`, `b ≠ 1` or `m` even); otherwise `Some(membership)`.  A
    /// point whose coordinates are not in `F_{2^m}` or that is not on the
    /// curve is reported as `Some(false)`.
    ///
    /// *Derivation.* Doubling `Q = (u, v)`, `u ≠ 0`, gives
    /// `λ = u + v/u`, `x = λ² + λ`, `y = u² + (λ + 1)·x`.  So `P = (x, y)`
    /// is halvable iff `λ² + λ = x` is solvable, iff `Tr(x) = 0`, with
    /// roots `H(x)` and `H(x) + 1` (the two halves `Q`, `Q + (0, 1)`).
    /// A half has `u² = y + (λ + 1)·x`, so
    /// `Tr(u) = Tr(y) + Tr(x·H(x))` when `Tr(x) = 0`, and the half is
    /// itself halvable iff that is `0`.  The 2-torsion point `(0, 1)`
    /// needs no special case: `Tr(0) = 0` and `Tr(1) = m mod 2 = 1`, so it
    /// is rejected.  The order-4 points `(1, 0)`, `(1, 1)` fail
    /// `Tr(x) = 0`.
    ///
    /// *Evidence.* Re-derived and checked in this repository against
    /// `[#E/4]P = O` on every point for `m = 5, 7, …, 23` and on 10⁶
    /// points each over ECC2K-130's field (`m = 131`) and on sect233k1
    /// (`RESEARCH_ISOGENY_FIBERS_ECDLP.md` Prop. 14.2 and §17.2;
    /// `research/isogeny_fibers_20261010/round4/16_2_fast_4e/`).
    ///
    /// *Prior art; no novelty is claimed.* The halvability criterion
    /// `Tr(x) = Tr(a)` is classical (Knudsen 1999; the point-halving
    /// literature), and halving-based validation of points on binary
    /// curves with cofactor 2 or 4, avoiding a scalar multiplication, is
    /// known: US Patent 7,257,709 / US 2004/0114760, "Method and apparatus
    /// for performing validation of elliptic curve public keys";
    /// B. King, "A point compression method for elliptic curves defined
    /// over GF(2^n)", PKC 2004, LNCS 2947; P. Eagle, S. Galbraith,
    /// J. Ong, "Point compression for Koblitz elliptic curves", IACR
    /// ePrint 2009/086 (two-level halving for `#E = 4ℓ`).  The closed
    /// form here only inlines the half's trace,
    /// `Tr(u) = Tr(u²) = Tr(y) + Tr(x·H(x))`, so no square root or second
    /// half-trace is needed.
    ///
    /// Not constant-time (table lookups indexed by the coordinates,
    /// early exit on `Tr(x)`): intended for public points only.
    pub fn is_in_4e(&self, p: &BinaryPoint) -> Option<bool> {
        if !self.is_odd_degree_koblitz_a0() {
            return None;
        }
        Some(match p {
            BinaryPoint::Infinity => true,
            BinaryPoint::Affine { x, y } => {
                if x.m_value() != self.m || y.m_value() != self.m || !self.is_on_curve(p) {
                    return Some(false);
                }
                let irr = &self.irreducible;
                if x.trace(irr) != 0 {
                    return Some(false);
                }
                let h = x
                    .half_trace(irr)
                    .expect("m is odd, so the half-trace exists");
                // Tr is linear: Tr(y) + Tr(x·H(x)) = Tr(y + x·H(x)).
                y.add(&x.mul(&h, irr)).trace(irr) == 0
            }
        })
    }
}

/// The order of `4E` on `y² + xy = x³ + 1` over an odd-degree binary field.
///
/// Over `F_2`, `#E = 4`, so the Frobenius traces satisfy `t_0 = 2`,
/// `t_1 = -1`, and `t_m = -t_{m-1} - 2t_{m-2}`.  Thus
/// `#E(F_{2^m}) = 2^m + 1 - t_m`.  For odd `m` its 2-part is exactly 4.
/// The cache depends only on the field degree, not its polynomial basis.
fn koblitz_a0_4e_order(m: u32) -> Arc<BigUint> {
    type Cache = Mutex<HashMap<u32, Arc<BigUint>>>;
    static CACHE: OnceLock<Cache> = OnceLock::new();
    let mut orders = CACHE
        .get_or_init(|| Mutex::new(HashMap::new()))
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    orders
        .entry(m)
        .or_insert_with(|| {
            let (mut previous, mut trace) = (BigInt::from(2u32), -BigInt::one());
            for _ in 1..m {
                let next = -&trace - (&previous << 1u32);
                previous = trace;
                trace = next;
            }
            let count = (BigInt::one() << m) + BigInt::one() - trace;
            Arc::new(count.to_biguint().expect("positive Koblitz curve order") >> 2u32)
        })
        .clone()
}

/// Negation: `−(x, y) = (x, y + x)`.
pub fn point_neg(p: &BinaryPoint) -> BinaryPoint {
    match p {
        BinaryPoint::Infinity => BinaryPoint::Infinity,
        BinaryPoint::Affine { x, y } => BinaryPoint::Affine {
            x: x.clone(),
            y: y.add(x),
        },
    }
}

/// **Point Addition** on a binary curve (Putranto et al. §3, formula
/// reproduced in the module header).
///
/// Cases:
/// - One operand is `O` → return the other.
/// - `P₁ = P₂` → defer to [`point_double`].
/// - `P₁ = −P₂` → return `O`.
/// - Otherwise: λ = (y₁+y₂)/(x₁+x₂), x₃ = λ² + λ + x₁ + x₂ + a,
///   y₃ = λ(x₁+x₃) + x₃ + y₁.
pub fn point_add(curve: &BinaryCurve, p1: &BinaryPoint, p2: &BinaryPoint) -> BinaryPoint {
    match (p1, p2) {
        (BinaryPoint::Infinity, q) | (q, BinaryPoint::Infinity) => q.clone(),
        (BinaryPoint::Affine { x: x1, y: y1 }, BinaryPoint::Affine { x: x2, y: y2 }) => {
            if x1 == x2 {
                // x₁ = x₂.  Either P + P (double) or P + (-P) = O.
                if y1.add(y2) == *x1 {
                    // y₂ = y₁ + x₁ ⇒ P₂ = -P₁.
                    return BinaryPoint::Infinity;
                }
                return point_double(curve, p1);
            }
            let irr = &curve.irreducible;
            let num = y1.add(y2);
            let den = x1.add(x2);
            let lambda = num.mul(&den.flt_inverse(irr).expect("non-zero"), irr);
            let lambda_sq = lambda.square(irr);
            // x₃ = λ² + λ + x₁ + x₂ + a
            let x3 = lambda_sq.add(&lambda).add(x1).add(x2).add(&curve.a);
            // y₃ = λ·(x₁ + x₃) + x₃ + y₁
            let y3 = lambda.mul(&x1.add(&x3), irr).add(&x3).add(y1);
            BinaryPoint::Affine { x: x3, y: y3 }
        }
    }
}

/// **Point Doubling** on a binary curve.
///
/// `[2]·(x₁, y₁)`:
/// - If `x₁ = 0` → result is `O`.
/// - Otherwise: λ = x₁ + y₁/x₁, x₃ = λ² + λ + a, y₃ = x₁² + (λ + 1)·x₃.
pub fn point_double(curve: &BinaryCurve, p: &BinaryPoint) -> BinaryPoint {
    match p {
        BinaryPoint::Infinity => BinaryPoint::Infinity,
        BinaryPoint::Affine { x: x1, y: y1 } => {
            let irr = &curve.irreducible;
            if x1.is_zero() {
                return BinaryPoint::Infinity;
            }
            let x1_inv = x1.flt_inverse(irr).expect("non-zero x₁");
            // λ = x₁ + y₁/x₁
            let lambda = x1.add(&y1.mul(&x1_inv, irr));
            let lambda_sq = lambda.square(irr);
            let x3 = lambda_sq.add(&lambda).add(&curve.a);
            // y₃ = x₁² + (λ + 1)·x₃
            let lambda_plus_one = lambda.add(&F2mElement::one(curve.m));
            let y3 = x1.square(irr).add(&lambda_plus_one.mul(&x3, irr));
            BinaryPoint::Affine { x: x3, y: y3 }
        }
    }
}

/// **ECPM**: scalar multiplication `[k]·P` via left-to-right
/// double-and-add.  This is the operation analysed in Putranto
/// et al. for quantum-resource estimation; classical version here
/// for verification.
pub fn scalar_mul(curve: &BinaryCurve, p: &BinaryPoint, k: &BigUint) -> BinaryPoint {
    if k.is_zero() {
        return BinaryPoint::Infinity;
    }
    let bits = k.bits();
    let mut result = BinaryPoint::Infinity;
    for i in (0..bits).rev() {
        result = point_double(curve, &result);
        if k.bit(i) {
            result = point_add(curve, &result, p);
        }
    }
    result
}

/// Brute-force search for any valid point and its order, on a tiny
/// test curve.  Only feasible for `m ≤ 16`.
fn find_any_generator(
    m: u32,
    irr: &IrreduciblePoly,
    a: &F2mElement,
    b: &F2mElement,
) -> (BinaryPoint, BigUint) {
    let curve_for_check = BinaryCurve {
        m,
        irreducible: irr.clone(),
        a: a.clone(),
        b: b.clone(),
        generator: BinaryPoint::Infinity,
        order: BigUint::one(),
        cofactor: BigUint::one(),
    };
    for xi in 1..(1u64 << m) {
        let x = F2mElement::from_biguint(&BigUint::from(xi), m);
        // Find y satisfying y² + xy = x³ + ax² + b.
        // Brute force: try every possible y.
        for yi in 1..(1u64 << m) {
            let y = F2mElement::from_biguint(&BigUint::from(yi), m);
            let pt = BinaryPoint::Affine { x: x.clone(), y };
            if curve_for_check.is_on_curve(&pt) {
                // Determine order by repeated doubling.
                let ord = point_order(&curve_for_check, &pt);
                if ord >= BigUint::from(2u32) {
                    return (pt, ord);
                }
            }
        }
    }
    (BinaryPoint::Infinity, BigUint::one())
}

fn point_order(curve: &BinaryCurve, p: &BinaryPoint) -> BigUint {
    let mut k = BigUint::one();
    let mut acc = p.clone();
    while k < BigUint::from(1u32 << 10) {
        acc = point_add(curve, &acc, p);
        k += 1u32;
        if matches!(acc, BinaryPoint::Infinity) {
            return k;
        }
    }
    // Couldn't find order in budget; return placeholder.
    k
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn certicom_challenge_points_and_subgroups() {
        for name in [
            "ECC2K-130",
            "ECC2-131",
            "ECC2K-163",
            "ECC2-163",
            "ECC2-191",
            "ECC2K-238",
            "ECC2-238",
            "ECC2K-358",
            "ECC2-353",
        ] {
            let (curve, target) = BinaryCurve::certicom_challenge(name).unwrap();
            assert!(curve.is_on_curve(&curve.generator), "{name} generator");
            assert!(curve.is_on_curve(&target), "{name} target");
            assert!(
                matches!(
                    scalar_mul(&curve, &curve.generator, &curve.order),
                    BinaryPoint::Infinity
                ),
                "{name} generator order"
            );
            assert!(
                matches!(
                    scalar_mul(&curve, &target, &curve.order),
                    BinaryPoint::Infinity
                ),
                "{name} target order"
            );
        }
    }

    /// Historical binary-curve generators are on their curves.
    #[test]
    fn legacy_binary_generators_are_on_curve() {
        for (name, curve) in BinaryCurve::legacy_low_security_family() {
            assert!(
                curve.is_on_curve(&curve.generator),
                "generator off curve for {}",
                name
            );
        }
    }

    #[test]
    fn low_security_binary_base_points_match_published_constants() {
        let cases = [
            (
                "sect113r1",
                BinaryCurve::sect113r1(),
                "009D73616F35F4AB1407D73562C10F",
                "00A52830277958EE84D1315ED31886",
            ),
            (
                "sect113r2",
                BinaryCurve::sect113r2(),
                "01A57A6A7B26CA5EF52FCDB8164797",
                "00B3ADC94ED1FE674C06E695BABA1D",
            ),
            (
                "sect131r1",
                BinaryCurve::sect131r1(),
                "0081BAF91FDF9833C40F9C181343638399",
                "078C6E7EA38C001F73C8134B1B4EF9E150",
            ),
            (
                "sect131r2",
                BinaryCurve::sect131r2(),
                "0356DCD8F2F95031AD652D23951BB366A8",
                "0648F06D867940A5366D9E265DE9EB240F",
            ),
            (
                "ike-oakley-group3",
                BinaryCurve::ike_oakley_group3(),
                "7B",
                "01C8",
            ),
        ];
        for (name, curve, gx, gy) in cases {
            let (x, y) = match curve.generator {
                BinaryPoint::Affine { x, y } => (x, y),
                BinaryPoint::Infinity => panic!("{} generator is infinity", name),
            };
            assert_eq!(
                x,
                F2mElement::from_hex(gx, curve.m),
                "Gx mismatch for {}",
                name
            );
            assert_eq!(
                y,
                F2mElement::from_hex(gy, curve.m),
                "Gy mismatch for {}",
                name
            );
        }
    }

    #[test]
    fn low_security_binary_catalog_security_categories() {
        for (name, curve) in BinaryCurve::sect113_family() {
            assert_eq!(curve.m, 113, "{} should be over F_2^113", name);
            assert_eq!(curve.order.bits() / 2, 56, "{} should be ~56-bit", name);
            assert_eq!(curve.cofactor, BigUint::from(2u32), "{} cofactor", name);
        }
        for (name, curve) in BinaryCurve::sect131_family() {
            assert_eq!(curve.m, 131, "{} should be over F_2^131", name);
            assert!(
                (64..=65).contains(&(curve.order.bits() / 2)),
                "{} should be ~64-bit",
                name
            );
            assert_eq!(curve.cofactor, BigUint::from(2u32), "{} cofactor", name);
        }

        let oakley = BinaryCurve::ike_oakley_group3();
        assert_eq!(oakley.m, 155);
        assert_eq!(oakley.irreducible.low_terms, vec![0, 62]);
        assert!(
            (76..=78).contains(&(oakley.order.bits() / 2)),
            "Oakley Group 3 should be roughly 76-bit"
        );
    }

    /// Adding `P + (-P) = O`.
    #[test]
    fn point_plus_neg_is_infinity() {
        let curve = BinaryCurve::sect163k1();
        let neg_g = point_neg(&curve.generator);
        let sum = point_add(&curve, &curve.generator, &neg_g);
        assert!(matches!(sum, BinaryPoint::Infinity));
    }

    /// `2P` via doubling equals `P + P` via addition.
    #[test]
    fn double_equals_self_add() {
        let curve = BinaryCurve::sect163k1();
        let g = curve.generator.clone();
        let g2_dbl = point_double(&curve, &g);
        let g2_add = point_add(&curve, &g, &g);
        assert_eq!(g2_dbl, g2_add);
    }

    /// `2P` is on the curve.
    #[test]
    fn double_stays_on_curve() {
        let curve = BinaryCurve::sect163k1();
        let g2 = point_double(&curve, &curve.generator);
        assert!(curve.is_on_curve(&g2));
    }

    /// Scalar multiplication: `[1]·G = G`, `[2]·G = 2G`, etc.
    #[test]
    fn scalar_mul_small_consistency() {
        let curve = BinaryCurve::sect163k1();
        let g = curve.generator.clone();
        let g1 = scalar_mul(&curve, &g, &BigUint::one());
        assert_eq!(g1, g);
        let g2_sm = scalar_mul(&curve, &g, &BigUint::from(2u32));
        let g2_dbl = point_double(&curve, &g);
        assert_eq!(g2_sm, g2_dbl);
        let g3_sm = scalar_mul(&curve, &g, &BigUint::from(3u32));
        let g3_add = point_add(&curve, &g2_dbl, &g);
        assert_eq!(g3_sm, g3_add);
    }

    /// `[5]G = 2·(2G) + G` via mixed PA/PD.  Demonstrates ECPM
    /// matches a hand-computed sequence.
    #[test]
    fn scalar_mul_matches_handcomputed() {
        let curve = BinaryCurve::sect163k1();
        let g = curve.generator.clone();
        let g2 = point_double(&curve, &g);
        let g4 = point_double(&curve, &g2);
        let g5 = point_add(&curve, &g4, &g);
        let g5_sm = scalar_mul(&curve, &g, &BigUint::from(5u32));
        assert_eq!(g5, g5_sm);
    }

    #[test]
    fn sect163_base_points_match_sec2_constants() {
        let cases = [
            (
                "sect163k1",
                BinaryCurve::sect163k1(),
                "02fe13c0537bbc11acaa07d793de4e6d5e5c94eee8",
                "0289070fb05d38ff58321f2e800536d538ccdaa3d9",
            ),
            (
                "sect163r1",
                BinaryCurve::sect163r1(),
                "0369979697AB43897789566789567F787A7876A654",
                "00435EDB42EFAFB2989D51FEFCE3C80988F41FF883",
            ),
            (
                "sect163r2",
                BinaryCurve::sect163r2(),
                "03F0EBA16286A2D57EA0991168D4994637E8343E36",
                "00D51FBC6C71A0094FA2CDD545B11C5C0C797324F1",
            ),
        ];
        for (name, curve, gx, gy) in cases {
            let (x, y) = match curve.generator {
                BinaryPoint::Affine { x, y } => (x, y),
                BinaryPoint::Infinity => panic!("{} generator is infinity", name),
            };
            assert_eq!(
                x,
                F2mElement::from_hex(gx, curve.m),
                "Gx mismatch for {}",
                name
            );
            assert_eq!(
                y,
                F2mElement::from_hex(gy, curve.m),
                "Gy mismatch for {}",
                name
            );
        }
    }

    #[test]
    fn sect163_nominal_security_category_is_80_bit() {
        for (name, curve) in BinaryCurve::sect163_family() {
            assert_eq!(curve.m, 163, "{} should be over F_2^163", name);
            assert_eq!(
                curve.order.bits() / 2,
                81,
                "{} has expected ~80-bit order",
                name
            );
            assert_eq!(curve.cofactor, BigUint::from(2u32), "{} cofactor", name);
        }
    }

    #[test]
    fn sect163_order_times_generator_is_identity() {
        for (name, curve) in BinaryCurve::sect163_family() {
            let n_g = scalar_mul(&curve, &curve.generator, &curve.order);
            assert!(
                matches!(n_g, BinaryPoint::Infinity),
                "n*G should be infinity for {}",
                name
            );
        }
    }

    #[test]
    fn low_security_order_times_generator_is_identity() {
        for (name, curve) in [
            ("sect113r1", BinaryCurve::sect113r1()),
            ("sect113r2", BinaryCurve::sect113r2()),
            ("sect131r1", BinaryCurve::sect131r1()),
            ("sect131r2", BinaryCurve::sect131r2()),
            ("ike-oakley-group3", BinaryCurve::ike_oakley_group3()),
        ] {
            let n_g = scalar_mul(&curve, &curve.generator, &curve.order);
            assert!(
                matches!(n_g, BinaryPoint::Infinity),
                "n*G should be infinity for {}",
                name
            );
        }
    }

    #[test]
    fn sect113r1_public_point_validation_accepts_generator() {
        let curve = BinaryCurve::sect113r1();
        assert!(curve.is_valid_public_point(&curve.generator));
    }

    #[test]
    fn sect113r1_public_point_validation_rejects_same_a_singular_point() {
        let curve = BinaryCurve::sect113r1();
        let singular_point = BinaryPoint::Affine {
            x: curve.a.clone(),
            y: F2mElement::zero(curve.m),
        };

        let singular_curve = BinaryCurve {
            b: F2mElement::zero(curve.m),
            ..curve.clone()
        };
        assert!(singular_curve.is_on_curve(&singular_point));
        assert!(!curve.is_on_curve(&singular_point));
        assert!(!curve.is_valid_public_point(&singular_point));
    }

    #[test]
    fn sect113r1_public_point_validation_rejects_nominal_order_two_point() {
        let curve = BinaryCurve::sect113r1();
        let order_two = BinaryPoint::Affine {
            x: F2mElement::zero(curve.m),
            // In characteristic two, sqrt(b) = b^(2^(m-1)).
            y: curve.b.square_k_times(curve.m - 1, &curve.irreducible),
        };

        assert!(curve.is_on_curve(&order_two));
        assert!(matches!(
            scalar_mul(&curve, &order_two, &BigUint::from(2u32)),
            BinaryPoint::Infinity
        ));
        assert_eq!(scalar_mul(&curve, &order_two, &curve.order), order_two);
        assert!(!curve.is_valid_public_point(&order_two));
    }

    #[test]
    fn public_point_validation_rejects_malformed_coordinate_width_without_panicking() {
        let curve = BinaryCurve::sect113r1();
        let (x, y) = match &curve.generator {
            BinaryPoint::Affine { x, y } => (x, y),
            BinaryPoint::Infinity => unreachable!("sect113r1 generator is affine"),
        };
        let malformed = BinaryPoint::Affine {
            x: F2mElement::from_biguint(&x.to_biguint(), curve.m + 1),
            y: F2mElement::from_biguint(&y.to_biguint(), curve.m + 1),
        };

        assert!(!curve.is_valid_public_point(&malformed));
    }

    // ── closed-form 4E test (`is_in_4e`) ──────────────────────────

    /// `[order]P = O` on the curve: the reference the closed form replaces.
    fn subgroup_by_scalar(curve: &BinaryCurve, p: &BinaryPoint) -> bool {
        curve.is_on_curve(p) && matches!(scalar_mul(curve, p, &curve.order), BinaryPoint::Infinity)
    }

    /// `y² + xy = x³ + 1` over `F_2[z]/irr`, cofactor 4, the given order
    /// and generator.
    fn koblitz_a0(irr: IrreduciblePoly, order: BigUint, generator: BinaryPoint) -> BinaryCurve {
        let m = irr.degree;
        BinaryCurve {
            m,
            irreducible: irr,
            a: F2mElement::zero(m),
            b: F2mElement::one(m),
            generator,
            order,
            cofactor: BigUint::from(4u32),
        }
    }

    fn affine(x: F2mElement, y: F2mElement) -> BinaryPoint {
        BinaryPoint::Affine { x, y }
    }

    /// A point with abscissa `x != 0` on `y² + xy = x³ + 1`, `m` odd, if
    /// one exists: `y = w·x` with `w² + w = x + 1/x²`.
    fn lift_x(curve: &BinaryCurve, x: &F2mElement, flip: bool) -> Option<BinaryPoint> {
        let irr = &curve.irreducible;
        let xi = x.flt_inverse(irr)?;
        let c = x.add(&xi.square(irr));
        if c.trace(irr) != 0 {
            return None;
        }
        let w = c.half_trace(irr).expect("odd m");
        let mut y = w.mul(x, irr);
        if flip {
            y = y.add(x);
        }
        let p = affine(x.clone(), y);
        assert!(curve.is_on_curve(&p));
        Some(p)
    }

    /// `T = (0, 1)` (order 2) and `P4 = (1, 0)` (order 4, `2·P4 = T`).
    fn two_part_points(m: u32) -> (BinaryPoint, BinaryPoint) {
        (
            affine(F2mElement::zero(m), F2mElement::one(m)),
            affine(F2mElement::one(m), F2mElement::zero(m)),
        )
    }

    /// `#E(F_{2^n})` for `y² + xy = x³ + 1`: `2^n + 1 − s_n` with
    /// `s_0 = 2, s_1 = −1, s_{k+1} = −s_k − 2·s_{k−1}`.
    fn koblitz_a0_order(n: u32) -> u64 {
        let (mut s0, mut s1) = (2i64, -1i64);
        for _ in 1..n {
            let s2 = -s1 - 2 * s0;
            s0 = s1;
            s1 = s2;
        }
        ((1i64 << n) + 1 - s1) as u64
    }

    /// Every point of `E(F_{2^n})`, n = 5, 7, 9, 11, 13: the closed form,
    /// `is_valid_public_point` (fast path) and `[#E/4]P = O` agree, and
    /// all four cosets of `4E` (read off from `[#E/4]P ∈ {O, P4, T, −P4}`)
    /// are visited, each with `#E/4` points.
    #[test]
    fn is_in_4e_matches_scalar_check_on_every_point_small_odd_m() {
        let fields = [
            (5u32, vec![0u32, 2]),
            (7, vec![0, 1]),
            (9, vec![0, 4]),
            (11, vec![0, 2]),
            (13, vec![0, 1, 3, 4]),
        ];
        for (m, low) in fields {
            let irr = IrreduciblePoly {
                degree: m,
                low_terms: low,
            };
            let n_e = koblitz_a0_order(m);
            assert_eq!(n_e % 8, 4, "2-part of #E must be Z/4 for odd m");
            let r = BigUint::from(n_e / 4);
            let mut curve = koblitz_a0(irr, r.clone(), BinaryPoint::Infinity);
            assert!(curve.has_fast_4e_subgroup_check());

            let mut points = vec![BinaryPoint::Infinity, two_part_points(m).0];
            for v in 1..(1u64 << m) {
                let x = F2mElement::from_words(&[v], m);
                if let Some(p) = lift_x(&curve, &x, false) {
                    points.push(point_neg(&p));
                    points.push(p);
                }
            }
            assert_eq!(points.len() as u64, n_e, "point count, m = {m}");
            curve.generator = points
                .iter()
                .find(|p| !matches!(p, BinaryPoint::Infinity) && subgroup_by_scalar(&curve, p))
                .expect("4E has affine points")
                .clone();

            let (t, p4) = two_part_points(m);
            let cosets = [BinaryPoint::Infinity, p4.clone(), t, point_neg(&p4)];
            let mut per_coset = [0u64; 4];
            for p in &points {
                let rp = scalar_mul(&curve, p, &r);
                let k = cosets
                    .iter()
                    .position(|c| *c == rp)
                    .expect("[#E/4]P is 4-torsion");
                per_coset[k] += 1;
                let truth = k == 0;
                assert_eq!(curve.is_in_4e(p), Some(truth), "m = {m}, P = {p:?}");
                let affine_truth = truth && !matches!(p, BinaryPoint::Infinity);
                assert_eq!(curve.is_valid_public_point(p), affine_truth, "m = {m}");
            }
            assert_eq!(per_coset, [n_e / 4; 4], "coset sizes, m = {m}");
        }
    }

    /// Random points from each coset `Q, Q + T, Q + P4, Q − P4` (`Q = 4P`)
    /// plus the uniform `P` itself, on the ECC2K-130 curve
    /// (`F_2[z]/(z^131 + z^8 + z^3 + z^2 + 1)`) and sect233k1: the fast
    /// path agrees with `[n]P = O` on every one.
    #[test]
    fn is_in_4e_matches_scalar_check_on_random_coset_points_ecc2k130_sect233k1() {
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};

        let ecc2k130_order =
            BigUint::parse_bytes(b"680564733841876926932320129493409985129", 10).unwrap();
        let mut ecc2k130 = koblitz_a0(
            IrreduciblePoly::deg_131(),
            ecc2k130_order,
            BinaryPoint::Infinity,
        );
        let mut rng = StdRng::seed_from_u64(0x4e_4e);
        for (name, curve) in [("ECC2K-130", &mut ecc2k130)] {
            loop {
                let x = F2mElement::from_words(&[rng.gen(), rng.gen(), rng.gen()], curve.m);
                if let Some(p) = lift_x(curve, &x, rng.gen()) {
                    let g = point_double(curve, &point_double(curve, &p));
                    if !matches!(g, BinaryPoint::Infinity) {
                        curve.generator = g;
                        break;
                    }
                }
            }
            assert!(subgroup_by_scalar(curve, &curve.generator), "{name}");
        }
        let sect233k1 = BinaryCurve::sect233k1();

        for (name, curve) in [("ECC2K-130", &ecc2k130), ("sect233k1", &sect233k1)] {
            assert!(curve.has_fast_4e_subgroup_check(), "{name}");
            let m = curve.m;
            let nw = m.div_ceil(64) as usize;
            let (t, p4) = two_part_points(m);
            let p4_neg = point_neg(&p4);
            let mut per_coset = [0u32; 4];
            let mut uniform_accepted = 0u32;
            let mut sampled = 0u32;
            while sampled < 40 {
                let words: Vec<u64> = (0..nw).map(|_| rng.gen()).collect();
                let x = F2mElement::from_words(&words, m);
                let Some(p) = lift_x(curve, &x, rng.gen()) else {
                    continue;
                };
                sampled += 1;
                let truth = subgroup_by_scalar(curve, &p);
                assert_eq!(curve.is_valid_public_point(&p), truth, "{name} uniform");
                assert_eq!(curve.is_in_4e(&p), Some(truth), "{name} uniform");
                uniform_accepted += truth as u32;

                let q = point_double(curve, &point_double(curve, &p));
                let coset_points = [
                    q.clone(),
                    point_add(curve, &q, &p4),
                    point_add(curve, &q, &t),
                    point_add(curve, &q, &p4_neg),
                ];
                for (k, c) in coset_points.iter().enumerate() {
                    let truth = k == 0 && !matches!(c, BinaryPoint::Infinity);
                    assert_eq!(subgroup_by_scalar(curve, c), k == 0, "{name} coset {k}");
                    assert_eq!(curve.is_valid_public_point(c), truth, "{name} coset {k}");
                    assert_eq!(curve.is_in_4e(c), Some(k == 0), "{name} coset {k}");
                    per_coset[k] += 1;
                }
            }
            assert_eq!(per_coset, [40; 4], "{name}");
            // About a quarter of uniform points lie in 4E; both outcomes occur.
            assert!(
                uniform_accepted > 0 && uniform_accepted < 40,
                "{name}: {uniform_accepted}/40 uniform points accepted"
            );
        }
    }

    /// The library's other odd-degree a = 0, cofactor-4 Koblitz curves take
    /// the fast path, accept their generator and reject it shifted into the
    /// three other cosets, in agreement with `[n]P = O`.
    #[test]
    fn is_in_4e_on_library_koblitz_cofactor4_generators() {
        for (name, curve) in [
            ("sect239k1", BinaryCurve::sect239k1()),
            ("sect283k1", BinaryCurve::sect283k1()),
            ("sect409k1", BinaryCurve::sect409k1()),
            ("sect571k1", BinaryCurve::sect571k1()),
        ] {
            assert!(curve.has_fast_4e_subgroup_check(), "{name}");
            let (t, p4) = two_part_points(curve.m);
            let g = curve.generator.clone();
            assert!(curve.is_valid_public_point(&g), "{name}");
            assert!(subgroup_by_scalar(&curve, &g), "{name}");
            for shift in [&p4, &t] {
                let c = point_add(&curve, &g, shift);
                assert!(!curve.is_valid_public_point(&c), "{name}");
                assert_eq!(curve.is_in_4e(&c), Some(false), "{name}");
            }
        }
    }

    /// O, T = (0, 1), the order-4 points, off-curve and wrong-width points.
    #[test]
    fn is_in_4e_edge_cases() {
        let curve = BinaryCurve::sect233k1();
        let m = curve.m;
        let (t, p4) = two_part_points(m);
        assert!(curve.is_on_curve(&t) && curve.is_on_curve(&p4));

        // O is in 4E, but is not a valid public point.
        assert_eq!(curve.is_in_4e(&BinaryPoint::Infinity), Some(true));
        assert!(!curve.is_valid_public_point(&BinaryPoint::Infinity));

        // T = (0, 1): Tr(0) = 0 and Tr(1) = m mod 2 = 1, rejected.
        assert_eq!(curve.is_in_4e(&t), Some(false));
        assert!(!curve.is_valid_public_point(&t));
        assert!(!subgroup_by_scalar(&curve, &t));

        // (1, 0), (1, 1): order 4, Tr(1) = 1.
        for p in [p4.clone(), point_neg(&p4)] {
            assert_eq!(curve.is_in_4e(&p), Some(false));
            assert!(!curve.is_valid_public_point(&p));
        }

        // Off the curve: the generator with y + 1.
        let (gx, gy) = match &curve.generator {
            BinaryPoint::Affine { x, y } => (x.clone(), y.clone()),
            BinaryPoint::Infinity => unreachable!(),
        };
        let off = affine(gx.clone(), gy.add(&F2mElement::one(m)));
        assert!(!curve.is_on_curve(&off));
        assert_eq!(curve.is_in_4e(&off), Some(false));
        assert!(!curve.is_valid_public_point(&off));

        // Coordinates of the wrong width.
        let wide = affine(
            F2mElement::from_biguint(&gx.to_biguint(), m + 2),
            F2mElement::from_biguint(&gy.to_biguint(), m + 2),
        );
        assert_eq!(curve.is_in_4e(&wide), Some(false));
        assert!(!curve.is_valid_public_point(&wide));
    }

    /// A malformed declared order must not enable the fast path even when
    /// the shape and cofactor match.  The point remains in 4E, but it is not
    /// killed by the declared scalar, so public validation must reject it.
    #[test]
    fn inconsistent_declared_order_keeps_the_scalar_subgroup_check() {
        let curve = BinaryCurve::sect233k1();
        let generator = curve.generator.clone();
        assert!(curve.has_fast_4e_subgroup_check());
        assert_eq!(curve.is_in_4e(&generator), Some(true));
        for order in [BigUint::one(), &curve.order + 2u32] {
            let misdeclared = BinaryCurve {
                order,
                ..curve.clone()
            };
            assert_eq!(misdeclared.cofactor, BigUint::from(4u32));
            assert!(!misdeclared.has_fast_4e_subgroup_check());
            assert_eq!(misdeclared.is_in_4e(&generator), Some(true));
            assert!(!subgroup_by_scalar(&misdeclared, &generator));
            assert!(!misdeclared.is_valid_public_point(&generator));
        }
    }

    /// Curves outside the preconditions keep the `[n]P = O` path.
    #[test]
    fn curves_outside_preconditions_take_the_scalar_path() {
        // a = 1 (sect163k1, cofactor 2), random curves, even m.
        for (name, curve) in [
            ("sect163k1", BinaryCurve::sect163k1()),
            ("sect233r1", BinaryCurve::sect233r1()),
            ("sect113r1", BinaryCurve::sect113r1()),
        ] {
            assert!(!curve.has_fast_4e_subgroup_check(), "{name}");
            assert_eq!(curve.is_in_4e(&curve.generator), None, "{name}");
            assert!(curve.is_valid_public_point(&curve.generator), "{name}");
        }
        let even = BinaryCurve {
            a: F2mElement::zero(8),
            ..BinaryCurve::test_curve_f256()
        };
        assert!(!even.has_fast_4e_subgroup_check());
        assert_eq!(even.is_in_4e(&BinaryPoint::Infinity), None);

        // a = 0, b = 1, odd m, but declared cofactor 2: the closed form
        // would still say G ∈ 4E; with a deliberately wrong order (n + 2)
        // the [n]P = O path rejects G, which shows that path is taken.
        let k233 = BinaryCurve::sect233k1();
        let misdeclared = BinaryCurve {
            cofactor: BigUint::from(2u32),
            order: &k233.order + 2u32,
            ..k233.clone()
        };
        assert!(!misdeclared.has_fast_4e_subgroup_check());
        assert_eq!(misdeclared.is_in_4e(&misdeclared.generator), Some(true));
        assert!(!misdeclared.is_valid_public_point(&misdeclared.generator));
        assert!(k233.is_valid_public_point(&k233.generator));
    }

    /// `[n]·G = O` where `n` is the generator's order.  This is
    /// the critical Lagrange-theorem test that distinguishes a
    /// real elliptic-curve implementation from one with subtle
    /// bugs.  Expensive (~thousands of point ops at sect163k1
    /// scale); ignored by default but can be run via --ignored.
    #[test]
    #[ignore]
    fn scalar_mul_by_order_is_identity() {
        let curve = BinaryCurve::sect163k1();
        let n_g = scalar_mul(&curve, &curve.generator, &curve.order);
        assert!(matches!(n_g, BinaryPoint::Infinity));
    }
}
