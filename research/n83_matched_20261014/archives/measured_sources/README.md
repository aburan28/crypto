# Byte-identical relevant measured sources

These12 files copy the relevant worker/core sources and Cargo inputs from
measurement-source freeze c6e8daafb5a59d27e1be88a2f61f0bd6ee83ef03. The
machine protocol and build receipt pin their original-path SHA256 values.
The old speed worker is an excluded singleton equivalence control; the new
rho worker and stock IC worker are the measured methods. Cargo.toml is named
Cargo.toml.snapshot so this archive does not create another Cargo package.

The full source revision remains reproducible from Git. These selected copies
and executable/source bindings are not a proof of the complete compiled
dependency closure. Later delivery integration must not retime or relabel
these frozen binaries as measurements of a different source tree.
