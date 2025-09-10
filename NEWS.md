# skipTrack 0.2.0

## New Features

*`skipTrack.fit()` can now incorporate longitudinal data in X matrix. Users are no
longer limited to one row of data per individual in X, but can now have one row 
of data per cycle if desired.

## Minor improvements and fixes

*`skipTrack.results()` now produces a helpful error message if burnIn is less than fit length.

*Fit now sets all starting cij values to 1 which is much more logical and mixing occurs 
much quicker.

# skipTrack 0.1.2

*Patch to align with parameter name changes from `genMCMCDiag`. 

*No functionality changes.

# skipTrack 0.1.1

* Minor edits to paper and readme based on JOSS review. No functional package changes.

# skipTrack 0.1.0

# skipTrack 0.0.1

* Initial CRAN submission.
