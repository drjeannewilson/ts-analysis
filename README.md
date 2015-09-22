# ts-analysis
Analysis Macros



recoAnalysis.C, recoAnalysis2.C
  Takes the output of the full reconstruction and produces some plots showing the reconstruction resolutions. These are resolutions for CCQE events >1m from the wall, but don't have the other cuts which improve things, so are worse than generally reported values but good for quickly comparing results after running the reconstruction after some change. The recoAnalysis2 compares two different versions distinguished by a suffix applied to the reconstruction output files (which is set in the batch_reco script described above).

MakeSelections.C, MakeSelections.h
  Takes the output of the full reconstruction and creates files with the selections in the format that Raj uses for his sensitivity studies. Most of the variables are copied from the NEUT output files that I don't understand, the rest (all the variables actually set by this code) should be straightforward.

selections.C, selections.h
  Takes the output of the full reconstruction and putting it into a more compact and useful format for analysis. The input is all the reconstruction files produced from the production and the output is a single file with a summary of the reconstruction for all events. Note that be default it includes the events from both horn currents with no entry in the tree to distinguish them (the first 400,000 will be FHC, the rest RHC) or you can just change the loop over files to only include one or the other. See the .C file for details.

selection_test.C
  Produces a load of plots that are useful when optimising the seleciton cuts. Things like how the resolutions change depending on towall and dwall cuts, how the purity and efficiencies change depending the wall and PID cuts.

selectioncounts.C
  Calculates and outputs tables of numbers for events passing the various selection cuts

misCountEvents.C
  Takes all the reconstructed events and selects out just those which I was investigating as having bad ring counting - basically it takes CC1pi events which are reconstructed as single ring and pass the other cuts for the 1ring e sample. It puts these events into new files to work with.

hough_test.C
  Produces a few plots that help look at what's going on with the ring counting. The plots are of the Hough space after the transform (top) and what the rings look like when mapped onto a sphere from the reconstructed vertex at the point where the Hough transform is done (bottom) i.e. this is the image of the rings that is the input to the Hough transform - it's really just the event display but with the TITUS shape mapped onto a sphere. Left is where ROOT does a projection of the sphere and right is just theta vs phi, rotated so the middle is looking in the beam direction. The markers on the top right plot are for the true particle directions (triangle, circle, X and + for muons, pions, electrons, gammas) and a star for the largest peak of the Hough space. Ideally, the centres of the rings on the bottom plots should correspond to peaks in the top plots, and these would be where the markers show the true particle directions. If the Hough transform is failing then the peaks would be in the wrong place, or wouldn't be sharp peaks, so multiple peaks would merge together into one. If it looks good, but the event is still reconstructed as the wrong number of rings, then it's not the Hough transform itself which is bad but either the way we assign hits to a ring or the way we decide if a second peak corresponds to a genuine ring.i
