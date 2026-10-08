# Stage 201 implementation clarification

The protocol phrase “updates its row end to the same exact maximum used by the
current algebra” means that the candidate preserves the control's block-wide
row-end metadata after applying a shorter combination. The skipped suffix is
known zero, but retaining the control end ensures that later pivot selection,
table shapes, logical-XOR accounting, and every downstream schedule remain
identical. Stage 201 measures only the table-preparation and row-application
words directly removed by exact combination ends.

For one combination, the candidate first writes through the maximum of its
source combination end and the added pivot end, then scans backward across
that newly written range to retain the actual last non-zero word. The scan is
runtime work and is included in process wall and CPU measurements; it is not
misreported as a word XOR. Only table words that the control writes but the
candidate never writes enter `full_m4ri_trimmed_table_words`.
