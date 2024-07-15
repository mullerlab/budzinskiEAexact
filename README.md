# Accompanying code to "An exact mathematical description of computation with transient spatiotemporal dynamics in a complex-valued neural network"

This repository contains the source code from

Budzinski\*, Busch\*, Mestern, Martin, Liboni, Pasini, Mináč, Coleman, Inoue, and Muller (2024) An exact mathematical description of computation with transient spatiotemporal dynamics in a complex-valued neural network. *Communications Physics* \*equal contribution

This paper introduces a complex-valued neural network (cv-NN) for understanding spatiotemporal computation. The source code in this repository reproduces Figures 2 (``cvnn_xor_gate.m``) and 3 (``cvnn_memory_task.m``) from the manuscript. The repository also provides an example of designing arbitrary target states for computation in the cv-NN (``cvnn_general_framework.m``). Installation only requires adding the files to the MATLAB path (done automatically where needed). Each script requires only a few minutes to run on a standard desktop workstation.

## Testing

Tested on MATLAB (R2021a, R2023a, and R2024a) under macOS.

## Dependencies

The helper and analysis functions required to run this code are included in ./helper_functions/ and ./analysis/, respectively, and are added to the path as needed.
