Power Scaling Laws and Near-Field Behaviors of Massive MIMO and Intelligent Reflecting Surfaces
==================

This is a code package is related to the follow scientific article:

Emil Björnson and Luca Sanguinetti, “[Power Scaling Laws and Near-Field Behaviors of Massive MIMO and Intelligent Reflecting Surfaces](https://arxiv.org/pdf/2002.04960),” IEEE Open Journal of the Communications Society, to appear.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. *We encourage you to also perform reproducible research!*


## Abstract of Article

The use of large arrays might be the solution to the capacity problems in wireless communications. The signal-to-noise ratio (SNR) grows linearly with the number of array elements N when using Massive MIMO receivers and half-duplex relays. Moreover, intelligent reflecting surfaces (IRSs) have recently attracted attention since these can relay signals to achieve an SNR that grows as N^2, which seems like a major benefit. In this paper, we use a deterministic propagation model for a planar array of arbitrary size, to demonstrate that the mentioned SNR behaviors, and associated power scaling laws, only apply in the far-field. They cannot be used to study the regime where N → ∞. We derive an exact channel gain expression that captures three essential near-field behaviors and use it to revisit the power scaling laws. We derive new finite asymptotic SNR limits but also conclude that these are unlikely to be approached in practice. We further prove that an IRS-aided setup cannot achieve a higher SNR than an equal-sized Massive MIMO setup, despite its faster SNR growth. We quantify analytically how much larger the IRS must be to achieve the same SNR. Finally, we show that an optimized IRS does not behave as an "anomalous" mirror but can vastly outperform that benchmark.

## Content of Code Package

The article contains 8 simulation figures, numbered 2, 3, 6, 7, 8, 9, 10, 11. simulationFigureX generates Figure X. The package also contains two Matlab functions that are used by some of the scripts.

See each file for further documentation.


## Acknowledgements

E. Björnson was supported by ELLIIT and the Wallenberg AI, Autonomous Systems and Software Program (WASP). L. Sanguinetti was partially supported by the University of Pisa under the PRA 2018-2019 Research Project CONCEPT, and by the Italian Ministry of Education and Research (MIUR) in the framework of the CrossLab project (Departments of Excellence).


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
