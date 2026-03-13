# Zebrafish good frame ranges

Start and end frame indices for the usable portions of each zebrafish
swimming trial in `zfishdata`. Some trials contain multiple blocks of
good frames.

## Usage

``` r
zfish_goodframes
```

## Format

A data frame with 8 rows and 4 columns:

- File:

  File name of the corresponding trial in `zfishdata`

- Start:

  First usable frame index

- End:

  Last usable frame index

- Block:

  Block number within the trial
