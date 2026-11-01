Source code accompanying the manuscript "On the stripe width distribution of
regular arid vegetation patterns", by Karl Kästner.

Run the batch script to reproduce all tables and figures of the manuscript at
once. The script fetches required library functions from other github repositories.
A package which includes the dependencies is available for download at Zenodo.

All scripts, except those running the Rietkerk and Tlidi simulation, terminate
nearly instantly. The scripts sw_rietkerk_simulate and sw_tlidi_simulate run
several hundred simulations of which each takes approximately 30 min on a single
core of a AMD Ryzen 7250. Simulations can be run in parallel by simply running
the scripts in parallel Matlab instances on the same computer, or several
computers with shared output directory. The scripts takes automatically care
to run each simulation only once by reserving the output file before starting
the model run.

Test system:
MATLAB Version: 25.2.0.2998904 (R2025b)
Linux / Ubuntu 26.04 LTS
Lenovo ThinkPad E14 Gen 7

The source code is licenced under GPL 3. The licence does not extent over the
input satellite image "pattern-sw-utm.png" copyright by Google Maps, Airbus,
Maxar Technologies, 2026. Displayed and analyzed under the fair use policy.

