# Repeat-signature timing clarification

The frozen prose in `probe_01/protocol.json` and `probe_02/protocol.json` says that
repeat-signature checking occurs outside arm clocks. The archived driver source
shows that these checks occur **inside validation time and total time**. Their
cost was charged. The raw observations, formulas and results require no numerical
correction and remain unchanged.

The current discovery/full protocol wording describes the actual timer placement.
Serialization and destruction of returned diagnostic records remain outside the
arm clocks and inside worker process receipts. No observer time is subtracted.
