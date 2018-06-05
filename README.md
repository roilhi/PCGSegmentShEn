# PCGSegmentShEn
Routines to find cardiac sound segmentation from Shannons Envelope and Matching Pursuit

- The file SegmentPCGFromEnvelope_copy.m file is a script which contains an example of implementing the function to take the Shannons Envelope. It requires the installation of the MPTK matching pursuit package to perform the decomposition with the mpdecomp function. The dictionary is an .xml file which contains the features of the atoms to be used: 32 and 64 samples size and a FFT of 1024 points. This file shows as an output a function NewSeg which points the possible onsets or start of cardiac sound events (marked as 1) and the offsets or end of cardiac sound events (marked as -1). 

- The file ShannonSmoothEnvelope.m is a function which calculates the Shannon's envelope based on Varghees[1] work. The threshold can be modified according to its amplitude and the sampling frequency is the other input parameter.

- The file SegmentPCGFromEnvelope.m is a similar script than the copy, it has more plots to show which possible onsets and offsets have to be eliminated according to HS criteria. 

- The xml file gabor_32_64_FFT1024.xml contains 2 blocks of Gabor atoms to perform the MP decomposition using the MPTK toolbox. 
