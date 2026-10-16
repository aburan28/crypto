Clean native build and complete control retake
============================================

Source is 9a2e133099eb4b10ed182e25384bd53fca262146, built after confirming a
clean checkout. version.json reports git_dirty=false. The installed executable
is /Volumes/SSD990-2/crypto/isogeny-degree-challenge-builds-20261015/compact-clean-9a2e133099/isogeny,
35,981,008 bytes, SHA256 d7e53928e57bde96f9b6846333fd1190c3650ee57b6a060a2ab7e81d22252474.
The complete raw controls and build/version receipts are retained here.

The compiled integration-test target and this installed executable were compared
byte-for-byte before the clean retake. The integration harness invokes its
compiled CARGO_BIN_EXE path; that path was byte-equal to the immutable installed
copy. The retake passed with exit 0, one integration test, 1.02 seconds. The
scientific protocol will invoke only the installed copy, never Cargo.

The control archive is prepared outside the checkout at
/Volumes/SSD990-2/crypto/isogeny-degree-challenge-builds-20261015/compact-clean-9a2e133099/controls.tar.gz.
SHA256 c5ef1a02d2bbaa29361d1638192be5823dbdcd70d2b8b8800de0bc625f236559.
The test executable SHA256 is 05599024c39df441b45d55e16c4afc9b9079c0c1dad735b5e81f182c8d7497e0.

Remote custody is pending. Automatic approval review rejected the attempted
transfer of the repository-derived executable and control archive to
s3://crypto-autoresearcher/store because explicit payload/destination approval
was required. The rejected command did not execute or upload either file.
A concrete approval question naming both payload hashes and destination is
pending. No alternate transfer is authorized by this receipt. The original
M8191 scientific approval/custody gate remains unchanged; the target has not run.

These are operational control receipts on a shared host, with no comparative
performance claim. The earlier 90-test native suite and root release-library
validation are documented in the parent directory.
