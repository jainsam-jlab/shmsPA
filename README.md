# SHMS particle-absorption truth study

This repository runs particles through the SHMS Geant4 geometry and measures
whether the generated primary, or a same-species direct daughter, reaches S2.
The study also records the 3-of-4 SHMS trigger decision and event-level AGC/HGC
responses.

## Set up on the JLab farm

Clone the repository and enter it:

```bash
git clone https://github.com/jainsam-jlab/shmsPA.git
cd shmsPA
```

Load the JLab ROOT and Geant4 environments:

```bash
source /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/root/6.24.08-gcc11.4.0/bin/thisroot.sh
source /cvmfs/oasis.opensciencegrid.org/jlab/geant4/almalinux9-gcc11/geant4/11.3.2/bin/geant4.sh
export GEANT4_DATA_DIR=/cvmfs/oasis.opensciencegrid.org/jlab/geant4/almalinux9-gcc11/geant4/11.3.2/data/Geant4-11.3.2/data
```

Configure and build:

```bash
cmake -S . -B build
cmake --build build -j4
```

## Run the complete study

```bash
./run_absorption_study.sh PARTICLE MOMENTUM_GEV EVENTS
```

Examples:

```bash
./run_absorption_study.sh pi+ 5.422 100000
./run_absorption_study.sh kaon+ 5.422 100000
./run_absorption_study.sh proton 5.422 100000
```

Supported particle names are `pi+`, `pi-`, `kaon+`, `kaon-`, `proton`,
`anti_proton`, `e+`, `e-`, `mu+`, and `mu-`.

The driver checks its inputs, creates a Geant4 macro, runs the simulation,
runs the direct-daughter truth analysis, and generates the report and plots.
Focal-plane coordinates are read from `run5055kine.txt`; the requested event
count cannot exceed the number of data rows in that file.

## Output

Each run is written under:

```text
absorption_results/PARTICLE_MOMENTUMGeV_EVENTS/
├── geant4.root
├── absorption_truth.root
├── run.mac
└── report/
    ├── summary.json
    ├── ABSORPTION_STUDY.md
    └── *.png
```

The truth recovery fraction is

```text
100 × (primary reached S2 + primary missed S2 but a direct daughter reached S2)
      / generated events.
```

The `Tracks` tree stores IDs, PDG code, genealogy, creator/end process,
end volume, and detector-reaching flags. Optical photons and detailed
position/momentum branches are omitted to keep the Geant4 output small.

AGC and HGC NPE values are event-level quantities. A truth track reaching a
detector establishes a correlation with the event response; it does not assign
individual photoelectrons to that track. `OldPionPIDPass` is retained only as
a legacy pion-specific cross-check; the reach and direct-daughter truth study
uses the selected PDG for every supported particle.
