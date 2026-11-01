# Pass rule: transcription of the source values

Written 2026-10-31, before the script first ran.

**NBS rows:** every value in `data/constants/nbs.yaml` must equal, character for character, the same cell of NIST's
digitized tables (`NBS_Tables_with_ErratumValues.csv`, doi:10.18434/M32124), whose SHA-256 must match the checksum
NIST publishes beside it. The CSV stays in the git-ignored `cache/nbs/`; the summary records the comparison. Each row
was also read against the scanned page of the 1982 tables (pp. 2-37 and 2-38).

**Equations:** every correlation that prints check values must reproduce them within the tolerance stored with it:
two parts per million for the IAPWS equations (they print six or more figures) and one unit of the last printed
digit for Battino et al.'s Table 1 column c.

**Verdict:** passed when every value matches.
