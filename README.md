# dijet-resonance-search

Dijet Resonance Search Project

Overview

This project aims to investigate dijet invariant mass distributions resulting from proton-proton collisions, with the ultimate goal of searching for potential new physics signals such as heavy resonance particles. By analyzing the properties of jets produced in high-energy collisions, the study focuses on refining methods to accurately calculate invariant masses and improve sensitivity for possible resonance signals.

Project Goals

○ Develop and validate simulation tools for proton-proton collisions, including jet clustering and invariant mass calculation.

○ Prepare the framework for signal injection (e.g., hypothetical Z′-like resonances) to evaluate analysis sensitivity.

○ Integrate and optimize ROOT-based macro tools for data processing, visualization, and analysis.

Current Status: Demonstration and Initial Validation

This repository includes a functional simulation along with ROOT macros that process the simulation output. The macros generate event-wise comparisons of invariant masses for the two highest-energy jets (leading and subleading). In addition, the dijet invariant masses are computed using both the four-vector summation method and an alternative collider-experiment formula for cross-validation.

This initial demonstration confirms the reliability of the calculation methods and the analysis workflow, providing foundational plots and data samples to support further development and review.

Repository Structure

src/ — Source code files for the simulation and ROOT macro.

data/ — Sample output data files from both the simulation and the ROOT macro.

plots/ — Graphical outputs illustrating invariant mass comparisons.

Next Steps

Planned development includes increasing event statistics, implementing signal injection with defined resonance parameters, and performing sensitivity studies based on statistical analysis metrics such as signal-to-background ratio or likelihood methods. The analysis framework will be finalized with a consistent invariant mass definition, and results will be compiled in a detailed report.
