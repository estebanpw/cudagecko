# CUDAGECKO

CUDAGECKO is a GPU implementation for the large-scale, seed-and-extend GECKO algorithm. 

## Requirements

CUDAGECKO requires a CUDA capable device that supports warp shuffle instructions.

Known to work with CUDA 10.0 and maxwell devices with compute capabilitiy 5.2

## Installation

Clone the repository and run:

`make install -C cudagecko/cuda`

## Use

`cudagecko/bin/gpu_cuda_workflow -query inputQueryFasta -ref inputReferenceFasta -dev 0`

For a list of parameters, use --help

Some useful parameters:

1. `-len`        <min len for an HSP to be reported, default: 32> NOTE: length runs in multiples of 32, so using 50 for instance will yield HSPs of 64 or more
2. `-dev`        <ID of the gpu device to use, default: 0>
3. `-factor`     <float between 0 and 1 to select the percentage of GPU memory reserved for words, default 0.125> NOTE: lower only if running out of memory for hits>
4. `-max_freq`   <only works in --sensitive Maximum frequency per hit (default: unlimited) (fast mode can skip highly repeated seeds)>
5. `--fast`      <Runs in fast mode as opposed to sensitive which is the default (this mode skips lots of repetitive seeds)>
6. `--sensitive` <Runs in sensitive mode (default). This mode is exhaustive, ALL seeds are calculated.> 
7. `--hyperfast` <Runs in hyperfast mode. This mode will match every word in one sequence to a different word in the other sequence. Only use to detect main syntenies.>

## Unattended execution

In case you want to run several comparisons in an unattended fashion (such as all vs all, see next section), it is recommended that you use the `gpu_cuda_wrapper.sh` script in the bin folder. This script will automatically lower the factor in case too many hits are found (and therefore will restart the execution of failed comparisons automatically). See FAQ at the end for more details.

Usage:

`cudagecko/bin/gpu_cuda_wrapper.sh <inputQueryFasta> <inputReferenceFasta> <minLen> <deviceID>`

Example usage:

`cudagecko/bin/gpu_cuda_wrapper.sh HOMSA.Chr.X.fasta MUSMU.Chr.X.fasta 256 0`

## Massive all vs all execution

In case you want to run GPUGECKO for all vs all experiments, use the following script:

`cudagecko/bin/gpu_two_species.sh <genomesDirectory1> <genomesDirectory2> <minLen> <device>`

This will compare all sequences in <genomesDirectory1> to all sequences of <genomesDirectory2>

## Extracting alignments

To get the alignments reported in the CSV, run the gpu_get_alignments tool in the binary folder such as this:

`cudagecko/bin/gpu_get_alignments csvFile inputQueryFasta inputReferenceFasta > output-alignments`

## Visualization

You can also explore the sequence comparison using the GECKO-MGV server available here: https://pistacho.ac.uma.es/
Simply upload your CSV to the portal.

## FAQ

- ERROR: `Reached maximum limit of hits (max 146459306)`. Use a smaller factor value (-factor, e.g. half of the current value). Default is 0.125 so you can lower in halves each time (0.125, 0.0625, 0.03, etc.)


- ERROR: `Not enough memory in pool`. This error will often appear in conjunction with the previous one. Use a smaller factor as well (see above).
