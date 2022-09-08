# Astrometric_Offset
The code estimates the astrometric offset between two fields.

These codes calculate the mean astrometric offset
between the 2 matched fields. They also plot the
offsets in a Sky-plot and create the histograms
of the offsets. 

The code uses Topcat in order to perform the match.

The repository contains:
1) Astrometric_Corr.py: This is the main code.
                          It uses the other code to
                          invoke Topcat and perform the match.
2) Topcat_Match.py      : This code calls Topcat.
