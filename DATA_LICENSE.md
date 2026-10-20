# Data source and license notice

## Original GazeBaseVR data

**Creators:** Dillon Lohr, Samantha Aziz, Lee Friedman and Oleg Komogortsev.

**Dataset:** *GazeBaseVR*, Figshare, Version 2.
[Versioned data record](https://doi.org/10.6084/m9.figshare.21308391.v2).

**License:** [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/),
as stated on the original Figshare record. The license permits reuse and redistribution
with appropriate attribution, a link to the license and an indication of changes.

**Data descriptor:** Lohr et al. (2023), *GazeBaseVR, a large-scale, longitudinal,
binocular eye-tracking dataset collected in virtual reality*, Scientific Data.
https://doi.org/10.1038/s41597-023-02075-5

The release asset `texran.zip` contains the 2,008 original task CSV files used here:
1,004 reading (TEX) and 1,004 random-saccade (RAN) recordings. The selection of the
two tasks is the subset operation; the contents of the original CSV files were not
modified. Individual CSV hashes were checked against the synchronized source files.
The archive does not include every task in the complete GazeBaseVR dataset.

These recordings are public, anonymized human eye-movement data. The present study
is a secondary analysis; no new participants were recruited and no new recordings collected.
The data descriptor documents the original collection and its ethics procedures.

## Derived analysis material

The derived outputs were produced by the extraction and analysis specifications
identified in the manuscript and repository. `derived_results.tar.gz` preserves the
source provenance paths; `DERIVED_FILE_MANIFEST.csv` records each retained file's
archive-relative path, size, SHA-256 digest and origin. Identical-content entries may
use standard tar hard links, preserving the retained paths after extraction.

The derived outputs' public availability does not remove the need to cite the original
data creators and data descriptor. This notice documents the upstream data license;
it does not grant an additional license for the authors' original code. All rights
reserved for the original code; contact the authors for licensing.
