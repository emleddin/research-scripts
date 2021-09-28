# Berendsen
This directory contains scripts that use the Berendsen thermostat.

`mineq.sh` is a script for looping through minimization to equilibration using
CPUs. It calls `mdin.1` through `mdin.10`.
`dynamicsgpu.sh` uses GPUs for production. It calls `mdin.11`.

You can read more about what's happening in `mdin.1` to `mdin.11` at:
[https://emleddin.github.io/comp-chem-website/AMBERguide-example-files.html#NVT](https://emleddin.github.io/comp-chem-website/AMBERguide-example-files.html#NVT)

