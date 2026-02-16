<p align="center">
  <img src="docs/images/LCBC_logo.png" width="500">
</p>

---

## Repository associated with the paper

# Cortical thickness changes precede high levels of amyloid by at least seven years

**Manuscript:**  
[Preprint (bioRxiv)](https://www.biorxiv.org/content/10.1101/2025.08.14.670398v1)  
_Status: In Review_

---

## Running the simulated example

The following scripts will run with provided simulated input:

### 01-prepSlopesYearsBeforeAB_simulated.r

**Purpose**  
Simulate longitudinal cortical thickness trajectories with GAMMs using MRIs at a given distance from Aβ+.  
Random intercepts and slopes are generated to emulate subject-specific variability around a nonlinear age trajectory.  
The simulated results are saved for subsequent analyses.

**Instructions**  
Set the following inputs to see examples of time cutoffs, group-matching procedures, and filtering based on scanner field strength:

- `yearsBeforeAB` — choose distance from AB+ in years  
- `matchFollow` — group-matching procedure  
- `singleT` — filter on scanner field strength


### 02-prepSlopesCompare_simulate.r

**Purpose**  
Load simulated GAMM random slopes and intercepts and test cortex-wide models comparing thickness effects between Aβ+ converter and Aβ− groups.

**Instructions**  
Set the following inputs to match and load the simulated results from `01-prepSlopesYearsBeforeAB_simulated.r`.  
The below example is provided:

```r
analysnum = 1
cohortSelect = "adnibacslcbc"
yearsBeforeAB = 1
matchFollow = 0
procType = "LONG"
predAB = 0
singleT = 0
analysname = paste0(
  "analysnum",
  analysnum,
  "matchF",
  matchFollow,
  "proc",
  procType,
  "predAB",
  predAB,
  "singleT",
  singleT,
  "reproduce1"
)
```

<br>
<hr>


## The following scripts reproduce paper results but require access to restricted individual-level data and are _not executable_. 

### 01-prepSlopesYearsBeforeAB.r

**Purpose**  
Compute longitudinal cortical thickness trajectories with GAMMs using MRIs at a given distance from Aβ+ combined with the full longitudinal sample.  
Random intercepts and slopes are estimated at each distance, and the results are saved for subsequent analyses.


---

### 02-prepSlopesCompare.r

**Purpose**  
Load GAMM random slopes and intercepts and test cortex-wide models comparing thickness effects between Aβ+ converter and Aβ− groups at given distances from observed/predicted Aβ+.

---

### 03-cortmaps_resample_nTime.r

**Purpose**  
Run resampling-based robustness analysis

---

### 04-cortmaps_summarize.r

**Purpose**  
Summarize results from resampling-based robustness check

---

### 06-regionalWildBootstrap.r

**Purpose**  
Run regional wild bootstrap resampling (guard against group differences in variance)

---

### 07-empiricalPmaps.r

**Purpose**  
Compute empirical p-value maps using null models from regional wild bootstrap resampling

---

### 08-rankAmyloidOrder.r

**Purpose**  
Calculate AB deposition maps and rank order of regional deposition

---

### 09-spinTestAB.r

**Purpose**  
Run spatial correlation between effect size maps at each cutoff and 3 Aβ deposition maps

---

### 10-rankThicknesChange.r

**Purpose**  
Compute the rank order of thickness changes across all converter MRIs and correlate with the rank order of regional Aβ deposition. 
Uses linear time-to-Aβ trajectories or SILA input, and computes time-to-Aβ thickness trajectories in regions of high-v-low Aβ.

---

## Master submit script

**NB!** `reproduce.sh` is the master script enabling submission to HPC for all analyses that are batched (those with an associated *batch script for HPC submission)
All other analyses were run locally.