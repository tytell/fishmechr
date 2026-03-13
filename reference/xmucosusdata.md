# Prickleback tracking data

Tracking data for a a swimming rock prickleback, *Xiphister mucosus*.
The data was tracked using Sleap (https://sleap.ai/) and comes out in
the following format. `frame_idx` is the frame number, and each point
along the body is identified with the point name and `.x`, `.y`, the
coordinate, and `.score`, which is a measure of the estimated accuracy
of the point. All of the points together are also given a score
(`instance.score`).

## Usage

``` r
xmucosusdata
```

## Format

A data frame with 711 rows and 27 columns. Columns follow the pattern
`<keypoint>.x`, `<keypoint>.y`, and `<keypoint>.score` for each tracked
keypoint, plus `frame_idx`, `instance.score`, and `track`. Keypoints
are: `Snout`, `BP1`–`BP6`, and `Tail`.
