# Exact affine-cone substitution profile

Keep the complete binary norm predicate, primary variables, constants and
AND dependencies unchanged. Treat every primary bit and every reachable
AND output as a formal independent atom, plus a constant-one atom. By
induction over the original DAG, each XOR node is an affine F2-form in
these atoms: constants/input/AND nodes have their obvious forms, and XOR
forms add by symmetric difference. No nonlinear term is linearized away;
the AND relation on its two affine input values is retained exactly.

Replace an affine AND input by a signed atom/constant if its support has
at most one nonconstant atom. Otherwise introduce one shared variable for
its canonical nonconstant form; constant-one offset is represented by a
literal sign. Its exact definition is one long XOR equation. Require each
original final-zero form directly with its original constant parity.
Constant false constraints must remain explicit; the profile records any
constant contradictions and never silently treats them as satisfaction.

For any assignment of primary variables, original reachable AND values
are uniquely determined in topological order. Each substituted input form
uses only prior AND atoms, hence its introduced variable is uniquely fixed
and the transformed AND relations reproduce precisely those values.
Conversely a transformed model reconstructs every removed XOR node by its
form and therefore satisfies the original predicate. Projection onto the
primary variables is a bijection; no norm-primary root can be lost or added.
This is classical affine variable elimination, not a novelty claim.

The worker counts formal long-XOR rows and variables and profiles their
incidence components. The old graph removes only explicit constants/output
units and is BEFORE further native unit propagation; its maxima do not
prove native rejection. Compressed maxima are formal bounds assuming long
rows are preserved; a concrete native interface may insert internal cuts,
so this is NOT native matrix admission or a SAT solve. Test the actual
long-row ABI, configuration, model mapping and all tiny roots before a
hard search. The existing pd_add_xor interface rejects rows longer than3.

Prediction: the compressed formal components fit the existing16384-row and
32768-column limits on the actual GF83 norm circuit. On refusal, redesign
representation instead of allocating a bigger matrix without evidence.
The profile retains the full nonlinear AND count and reports longest row,
total nonzero coefficients, all variables and constant contradictions.

Independent controls retain73216 literal-norm assignment checks and test
every16384 tiny primary assignment against full DAG evaluation for every
reachable affine form, using independent dot-product parity evaluation.
90seconds/4GiB/oneCPU, complete serialization checked and late flags demoted;
full supervisor lifetime charged. No root search, original tuple, speedup,
matched F4/F5/F6 or qualified8x claim is made by this representation profile.
