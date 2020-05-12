# CUDAGECKO

CUDAGECKO is a GPU implementation for the large-scale, seed-and-extend GECKO algorithm. 

## Installation

Clone the repository and make install each folder, or run ./INSTALL in the repository folder

## Use

cudagecko/bin/gpu_cuda_workflow -query inputQueryFasta -ref inputReferenceFasta

For a list of parameters, use --help
Some useful parameters:

1. -len <min len for an HSP to be reported, default: 32> NOTE: length runs in multiples of 32, so using 50 for instance will yield HSPs of 64 or more
2. -dev <ID of the gpu device to use, default: 0>
3. -factor <float between 0 and 1 to select the percentage of GPU memory reserved for words, default 0.125> NOTE: lower only if running out of memory for hits

## Extracting alignments

To get the alignments reported in the CSV, simply run the "gpu_get_alignments.sh" script in the binary folder such as this:

cudagecko/bin/gpu_get_alignments.sh csvFile inputQueryFasta inputReferenceFasta > output-alignments

## Visualization

You can also explore the sequence comparison using the GECKO-MGV server available here: https://pistacho.ac.uma.es/
Simply upload your CSV to the portal.
