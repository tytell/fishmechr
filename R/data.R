#' Lamprey midline data
#'
#' Midline tracking data for a sea lamprey (_Petromyzon marinus_) swimming
#' steadily. Contains 20 evenly-spaced points along the body from head to tail,
#' across 80 frames.
#'
#' @format A data frame with 1600 rows and 5 columns:
#'   \describe{
#'     \item{t}{Time in seconds}
#'     \item{frame}{Frame number}
#'     \item{point}{Point index along the body (1 = head, 20 = tail)}
#'     \item{mxmm}{x coordinate in mm}
#'     \item{mymm}{y coordinate in mm}
#'   }
"lampreydata"

#' Fish body width profiles
#'
#' Body width (diameter) as a function of fractional arc length for an eel
#' (_Anguilla rostrata_) and a sea lamprey (_Petromyzon marinus_). Used for
#' estimating the center of mass via [get_midline_center_df()].
#'
#' @format A data frame with 20 rows and 3 columns:
#'   \describe{
#'     \item{s}{Arc length as a fraction of body length (0 = head, 1 = tail)}
#'     \item{eelwidth}{Body width of the eel as a fraction of body length}
#'     \item{ammowidth}{Body width of the lamprey as a fraction of body length}
#'   }
"fishwidth"

#' Prickleback tracking data
#'
#' Tracking data for a a swimming rock prickleback, _Xiphister mucosus_.
#' The data was tracked using Sleap (https://sleap.ai/) and comes out in the following format.
#' `frame_idx` is the frame number, and each point along the body is identified with the point
#' name and `.x`, `.y`, the coordinate, and `.score`, which is a measure of the estimated
#' accuracy of the point. All of the points together are also given a score (`instance.score`).
#'
#' @format A data frame with 711 rows and 27 columns. Columns follow the
#'   pattern `<keypoint>.x`, `<keypoint>.y`, and `<keypoint>.score` for each
#'   tracked keypoint, plus `frame_idx`, `instance.score`, and `track`.
#'   Keypoints are: `Snout`, `BP1`--`BP6`, and `Tail`.
"xmucosusdata"

#' Zebrafish body shape
#'
#' Width and height profile of a zebrafish body as a function of fractional
#' arc length. Used for estimating volumes and centers of mass.
#'
#' @format A data frame with 20 rows and 3 columns:
#'   \describe{
#'     \item{s}{Arc length as a fraction of body length (0 = head, 1 = tail)}
#'     \item{width}{Lateral body width as a fraction of body length}
#'     \item{height}{Dorso-ventral body height as a fraction of body length}
#'   }
"zebrafish_shape"

#' Zebrafish good frame ranges
#'
#' Start and end frame indices for the usable portions of each zebrafish
#' swimming trial in `zfishdata`. Some trials contain multiple blocks of
#' good frames.
#'
#' @format A data frame with 8 rows and 4 columns:
#'   \describe{
#'     \item{File}{File name of the corresponding trial in `zfishdata`}
#'     \item{Start}{First usable frame index}
#'     \item{End}{Last usable frame index}
#'     \item{Block}{Block number within the trial}
#'   }
"zfish_goodframes"

#' Zebrafish keypoint tracking data
#'
#' Keypoint tracking data for zebrafish (_Danio rerio_) swimming at several
#' speeds. Keypoints are tracked using DeepLabCut and include body landmarks
#' and fin tips.
#'
#' @format A data frame with 1056 rows and 46 columns. Most columns follow
#'   the pattern `<keypoint>.x`, `<keypoint>.y`, and `<keypoint>.score` for
#'   each tracked keypoint. Additional columns:
#'   * `fn` Source file name
#'   * `frame_idx` Frame index within the trial
#'   * `id` Fish identifier
#'   * `speed` Swimming speed in cm/s
#'   * `datetime` Date and time of the trial
#'   * `instance.score` Overall instance detection score
#'   * `track` Track identifier
"zfishdata"
