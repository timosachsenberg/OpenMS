# Preliminary selection and calibration

One combined review-derived candidate per backend was tested against preserved originals; this is not an exhaustive search of optimization subsets. All candidates pass their affected class tests, focused regressions, benchmark checks and huge-array digest check. The preliminary matrix has 288 successful measurements: two inputs, six modes, 1/8 readers, three repetitions. It is separate from the final five-repetition results.

| Backend | Original stream 1 (s), median of 5 | Candidate stream 1 (s), median of 3 | Original stream peak at 8 (GiB) | Candidate peak at 8 (GiB) |
|---|---:|---:|---:|---:|
| xerces | 2.926 | 2.628 | 0.14 | 0.14 |
| libxml2 | 5.536 | 3.234 | 0.13 | 0.13 |
| expat | 6.776 | 4.665 | 0.15 | 0.15 |
| pugixml | 5.902 | 2.594 | 25.81 | 9.86 |

Retain the combined native candidates for full reruns: large streaming improvements and the pugixml RSS reduction justify them; unchanged fragment-query paths show small fluctuations that are not attributed to new parser algorithms. Correctness fixes remain regardless of noise in timing.

**Control caution:** old Xerces single-reader stream range 2.602–2.949 s overlaps candidate range 2.491–2.647 s. A fresh unchanged-core calibration produced 2.785 s (2.730–2.807). The unchanged control uses its original ABI-compatible benchmark executable; source of the timed harness and compiler flags match the reviewed build. There is no confirmed algorithmic reason for a full-stream Xerces gain. Do not attribute this shift to skipping metadata-only data or the unused UTF8 hook. Temporal scheduling/frequency and binary layout effects are not separated by these experiments.

Full results compare all four freshly run reviewed implementations, with the optimized Xerces control carrying shared fixes. Before/after ratios remain observations across sessions, with min/max retained; they are not causal confidence intervals. Calibration results are preserved separately and not mixed into any final medians.
