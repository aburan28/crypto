# Kernel-first continuation protocol

Frozen before target construction, 2026-10-09. The four full-Hecke probes in
PROTOCOL.md continue unchanged. This separate route attempts one specified
Frobenius eigenline, using the existing native Kohel/kernel-map algorithms.
It does not replace or satisfy the two-map acceptance rule of the full probes.

Exact screening through 65537 selected P-224 degree 1471, eigenvalue 554 of
order 5 (other eigenvalue 1178, order 490), and P-192 degree 10453, eigenvalue
270 of order 3 (other eigenvalue 8274, order 402). Use seed 1, 1800 seconds per
construction and the same native 8 GiB RSS cap. Run one construction at a time,
after the preregistered full-Hecke jobs have finished. Each selection and every
outcome stays preserved.

Build a native degree-5 or degree-3 extension over the source prime field.
Check its defining polynomial by Rabin's irreducibility criterion. Derive the
extension group order from the published source trace using the standard Lucas
recurrence. Obtain nonzero prime-degree torsion, verify its order and selected
Frobenius eigenvalue, and form its cyclic kernel. Every kernel coefficient must
lie in the base field. Construct the rational map with the existing API.

Field arithmetic must pass independent small-field and large-prime checks,
including inverse/Frobenius relations and rejection of reducible moduli.
Accept each map only after the existing walker's independent kernel, codomain,
exact rational-map and fresh public subgroup transport checks. The one-map
record and replay must explicitly retain degree coverage PARTIAL: one of the
two eigenlines is constructed, with the other still unresolved. A verified
individual map is a map-level result, without a claim of complete degree
enumeration or an ECDLP-cost change. Register any new models before citing them.
