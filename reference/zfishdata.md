# Zebrafish keypoint tracking data

Keypoint tracking data for zebrafish (*Danio rerio*) swimming at several
speeds. Keypoints are tracked using DeepLabCut and include body
landmarks and fin tips.

## Usage

``` r
zfishdata
```

## Format

A data frame with 1056 rows and 46 columns. Most columns follow the
pattern `<keypoint>.x`, `<keypoint>.y`, and `<keypoint>.score` for each
tracked keypoint. Additional columns:

- `fn` Source file name

- `frame_idx` Frame index within the trial

- `id` Fish identifier

- `speed` Swimming speed in cm/s

- `datetime` Date and time of the trial

- `instance.score` Overall instance detection score

- `track` Track identifier
