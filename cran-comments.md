## R CMD check results

0 errors | 0 warnings | 0 notes

Tested on:
- macOS Sequoia 15.7.4, R 4.5.2 (aarch64-apple-darwin20)
- Windows (win-builder, R-devel)

## Win-builder notes

The win-builder check produced 2 notes:

1. **Possibly misspelled words in DESCRIPTION**: The flagged words (Biomechanics,
   energetics, Higham, Santo, Tytell) are all either standard scientific terminology
   or author surnames from the book chapter citation in the Description field.
   They are spelled correctly.

2. **Non-standard file at top level**: `cran-comments.md` has been added to
   `.Rbuildignore` and will not be included in the package tarball.

## Notes

This is the initial CRAN submission.
