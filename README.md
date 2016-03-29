16SpeB
===

Example run for 16SpeB with *Mycoplasma hominis*

1. Compile palign
---

    ./compile_palign

Note: This will require g++ on your system. It's pre-installed on Mac OS X and
most Linux distribution. For installation of g++ on Debian/Ubuntu:

    sudo apt-get install g++

2. Extracting species
---

We already added *Mycoplasma hominis* data within the `relevant_species` folder.

3. Filtering regions
---

    ./isolate_multiregions.pl --midstr='ACTCCTACGGGAGGCAGCA' --rightstr='GTCGTCAGCTCGTGYYG' --rightcoord=1061 --midcoord=338 --right-trimneg=258 --right-trimpos=0 --mid-trimneg=270 --mid-trimpos=0 --longlen-trimneg=270 --longlen-trimpos=1000 --offsetleft=100 --offsetright=100 --fsa=relevant_species/Mycoplasma_hominis.fsa --longlen-fsa=test_Mycoplasma_hominis/Mycoplasma_hominis_raw_16s.fsa --mid-fsa=test_Mycoplasma_hominis/Mycoplasma_hominis_raw_v2.fsa --right-fsa=test_Mycoplasma_hominis/Mycoplasma_hominis_raw_v6.fsa

4. Removing exact duplicates within intra-species
---

    ./filter_duplicates_within_intra_species.pl --fsa=test_Mycoplasma_hominis/Mycoplasma_hominis_raw_16s.fsa --species=Mycoplasma_hominis --nobadwords --auxin=test_Mycoplasma_hominis/Mycoplasma_hominis_raw_v6.fsa --auxout=test_Mycoplasma_hominis/Mycoplasma_hominis_nondup_v6.fsa

5. Computing intra-species PID
---

    ./compute_intra_species_pid_distrib.pl --iters=100 --dir=test_Mycoplasma_hominis --genus=Mycoplasma
