# Clinical_fingerprinting
Sample code for reproducing the analysis reported in Sorrentino et al. "Clinical connectome fingerprints of cognitive decline", NeuroImage 2021.

The code compares two functional connectome (FC) acquisitions/sessions for each subject (defined as test and retest)
contained in the sample dataset (data_test_FC) and give as result an "identifiabilty matrix" (Amico & Goñi, Scientific Reports 2018). Its main diagonal (top left to bottom right) highlights the self-identifiability (Iself) of the sample, i.e. the similarity between the test and retest FC of the same subject. The remaining elements consists in the comparison between test and retest FC of different subjects (Iothers).

Clinical connetome fingerprinting extends this idea by comparing individual connectomes between different populations (i.e. Controls and Mild Cognitive Impairment, see Sorrentino et al. for details). Functional Connectomes matrices are obtained through Phase Linearity Measurement (PLM) (Baselice et al., 2019), here only reported in the Alpha band.

*PLEASE NOTE: For privacy reasons, we could not use the same data as in Sorrentino et al. 
Hence, the sample results reported here do not correspond to the original manuscript, instead they are merely illustrative of the methodology.
Please contact giuseppe.sorrentino@uniparthenope.it for (reasonably) requesting the original sample connectomes of the manuscript.
The data used in this code come from healthy MEG connectomes obtained from the Human Connectome Project dataset (see Sareen et al., NeuroImage 2021, for details).*

Authors: Emahnuel TROISI LOPEZ, Pierpaolo SORRENTINO & Enrico AMICO
 version 1.0. July 01, 2021

**PLEASE CITE US! If you are using this code for your research, please kindly cite us:
Authors: Sorrentino P, Rucco R, Lardone A, Liparoti M, Lopez ET, Cavaliere C, Soricelli A, Jirsa V, Sorrentino G, Amico E.
Title: Clinical connectome fingerprints of cognitive decline.
Published on: NeuroImage - 2021 Jun 9, p. 118253.
Doi: doi.org/10.1016/j.neuroimage.2021.118253**
