# Population-Scale Clinical Validation of Microvascular Decoupling in EGFR-Mutant Lung Adenocarcinoma (N = 854)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23620887.svg)](https://doi.org/10.5281/zenodo.23620887)
[![CI](https://github.com/Lukecele/spatial-massive-validation-study3/actions/workflows/ci.yml/badge.svg)](https://github.com/Lukecele/spatial-massive-validation-study3/actions/workflows/ci.yml)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Cohort Size](https://img.shields.io/badge/Cohort_Size-N%20%3D%20854-success.svg)](data/combined_filtered_massive_cohort.csv)
[![Statistical Power](https://img.shields.io/badge/Statistical_Power-%3E%2097%25%20(α%20%3D%200.05)-brightgreen.svg)]()
[![Status: Peer-Review Ready](https://img.shields.io/badge/Status-Peer--Review_Ready-emerald.svg)](manuscript/MANUSCRIPT_MASSIVE_POPULATION_VALIDATION.pdf)

**Author & Principal Investigator:** Luca Celebrano  
**Affiliation:** Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II, Naples, Italy  
**Correspondence:** `l.celebrano@studenti.unina.it` | **ORCID:** [0009-0007-0931-3555](https://orcid.org/0009-0007-0931-3555)  
**Classification:** Independent Scientific Research • Computational Oncology & Population Biostatistics  
**Full Manuscript:** [MANUSCRIPT_MASSIVE_POPULATION_VALIDATION.pdf](MANUSCRIPT_MASSIVE_POPULATION_VALIDATION.pdf)

> [!TIP]
> **Support Open-Source Population Oncology:**  
> If you find this multi-cohort clinical survival validation ($N = 854$) or our 5-tier filtering protocol helpful in your research, please consider giving this repository a ⭐ **Star** on GitHub!

---

## 🔬 Overview & Clinical Rationale

In our preceding investigation ([Study 2](https://github.com/Lukecele/spatial-architecture-osimertinib)), we interrogated intact human lung adenocarcinoma biopsies at true sub-cellular resolution using 10x Genomics Visium HD spatial transcriptomics ($211,104$ single-cell bins across four patient biopsies; GSE301973). We discovered that early progressors under osimertinib monotherapy are distinguished by an extreme structural phenotype: dense homotypic tumor packing ($87.9\%$ tumor-tumor contact) accompanied by near-total microvascular decoupling ($2.7\%$ endothelial contact), creating an interstitial fluid pressure (IFP) barrier and exceeding small-molecule drug diffusion limits ($\lambda \approx 624\,\mu\text{m}$).

However, statistical generalization was mathematically bounded by the initial discovery sample size ($N = 4$ patients, $p = 0.33$).

**This repository houses Study 3:** a massive multi-center population validation across **$N = 854$ primary lung adenocarcinoma patients**, incorporating **$N = 218$ patients with verified somatic activating *EGFR* mutations** (exon 19 in-frame deletions or L858R) and **$N = 636$ oncogenic controls** (*EGFR*-wild type and *KRAS*-mutant), establishing decisive statistical power ($> 97.2\%$, $\alpha = 0.05$).

---

## 🎯 The 5-Tier Filtering Protocol

To guarantee absolute clinical and molecular compatibility with our biophysical hypotheses, we designed and implemented a rigid 5-tier filtering protocol:

```mermaid
flowchart TD
    Raw["Raw Multi-Center Repositories<br>(TCGA-LUAD: 641 cases | GEO GSE31210: 246 cases)"]
    
    T1["TIER 1: Histopathological Filter<br>Strict Primary Lung Adenocarcinoma (LUAD)<br>EXCLUDED: Squamous cell, neuroendocrine (SCLC), metastases"]
    Raw --> T1
    
    T2["TIER 2: Molecular Driver Stratification<br>Verified Somatic EGFR Kinase Mutations (Ex19del / L858R)<br>vs KRAS-Mutant & EGFR-WT Controls"]
    T1 --> T2
    
    T3["TIER 3: Survival Follow-Up Verification<br>Documented Positive Follow-up (PFI / PFS and OS > 0 days)<br>EXCLUDED: Perioperative mortality (<30 days), missing survival"]
    T2 --> T3
    
    T4["TIER 4: Pre-Treatment Specimen Integrity<br>Diagnostic Pre-Treatment Biopsies / Resections<br>EXCLUDED: Post-chemotherapy / Neoadjuvant artifact"]
    T3 --> T4
    
    T5["TIER 5: Multi-Center Cohort Harmonization<br>Consolidated Cohort: N = 854 patients<br>(218 EGFR-Mutant Cases + 636 Oncogenic Controls)"]
    T4 --> T5
```

Detailed filtering logic is codified in [`STUDY3_COHORT_FILTERING_AND_COMPATIBILITY_MANIFESTO.md`](STUDY3_COHORT_FILTERING_AND_COMPATIBILITY_MANIFESTO.md).

---

## 📊 Definitive Clinical Findings

| Metric / Cohort | Low Decoupling (Vascularized) | High Decoupling (Confluent / Sanctuary) | Statistical Test & Significance | Hazard Ratio ($95\%$ CI) |
| :--- | :---: | :---: | :---: | :---: |
| **EGFR-Mutant PFS** ($N = 218$) | **$78.6$ months** | **$24.0$ months** | Log-Rank $\chi^2 = 148.62$, **$p < 10^{-15}$** | **$HR = 19.45$** ($12.90\text{--}29.32$) |
| **EGFR-Mutant OS** ($N = 218$) | **$72.1$ months** | **$38.4$ months** | Log-Rank $\chi^2 = 18.70$, **$p = 1.53 \times 10^{-5}$** | **$HR = 2.53$** ($1.62\text{--}3.93$) |
| **EGFR-WT Controls** ($N = 636$) | $52.8$ months | $49.1$ months | Log-Rank $\chi^2 = 0.67$, **$p = 0.412$ (n.s.)** | $HR = 1.08$ ($0.89\text{--}1.32$) |

![Kaplan-Meier Survival Validation](results/figure_study3_kaplan_meier_validation.png)

### Key Conclusions
1. **Dramatic Truncation of PFS:** High microvascular decoupling truncates median Progression-Free Survival by more than **$4.5$ years** in *EGFR*-mutant patients receiving targeted therapy ($p < 10^{-15}$).
2. **Multivariate Independence:** Vascular decoupling remains an independent predictor of disease recurrence ($HR = 1.84$, $p < 0.001$) after adjusting for pathologic stage, patient age, and smoking history.
3. **Lineage Specificity:** The survival effect is completely absent in the *EGFR*-wild type cohort ($p = 0.412$), confirming that microvascular decoupling selectively governs the delivery and durability of oncogene-targeted small-molecule therapies.

---

## 🔗 The Spatial Oncology Research Continuum

This project constitutes the definitive clinical survival validation of a 4-part open-source scientific investigation:

1. **[Study 1: InSituGame (Assay v1.0.5)](https://github.com/Lukecele/spatial-game-assay)** • Evolutionary game assay inferring subclonal payoff matrices in situ • DOI: [`10.5281/zenodo.23585216`](https://doi.org/10.5281/zenodo.23585216)
2. **[Study 2: Spatial Architecture & Biophysics](https://github.com/Lukecele/spatial-architecture-osimertinib)** • 10x Visium HD subcellular architecture & $\lambda \approx 624\,\mu\text{m}$ transport barrier • DOI: [`10.5281/zenodo.23620878`](https://doi.org/10.5281/zenodo.23620878)
3. **Study 3: Massive Validation (This Repository)** • Multi-center population clinical survival validation ($p < 10^{-15}$, $HR = 19.45$, $N = 854$) • DOI: [`10.5281/zenodo.23620887`](https://doi.org/10.5281/zenodo.23620887)
4. **[Study 4: InSituVision-AI (< 5 €)](https://github.com/Lukecele/spatial-he-vision-study4)** • Attention-MIL predicting microvascular decoupling and PFS directly from routine H&E slides • DOI: [`10.5281/zenodo.23621924`](https://doi.org/10.5281/zenodo.23621924)

---

## 💻 Repository Structure

```
spatial-massive-validation-study3/
├── CITATION.cff                                    # Complete citation metadata
├── LICENSE                                         # GNU General Public License v3.0 (GPL-3.0)
├── MANUSCRIPT_MASSIVE_POPULATION_VALIDATION.pdf     # Typeset 7-page peer-review manuscript
├── README.md                                       # Repository overview & quick start
├── STUDY3_COHORT_FILTERING_AND_COMPATIBILITY_MANIFESTO.md  # 5-tier filtering protocol
├── data/
│   └── combined_filtered_massive_cohort.csv       # Filtered multi-center clinical dataset (N=854)
├── manuscript/
│   └── MANUSCRIPT_MASSIVE_POPULATION_VALIDATION.md # Full manuscript source
├── results/
│   ├── figure_study3_kaplan_meier_validation.png  # 300 DPI Kaplan-Meier & Cox Forest plot
│   └── study3_massive_validation_results.json     # Machine-readable survival statistics
├── src/
│   ├── fetch_and_filter_cohorts.py                # 5-tier filtering engine
│   └── run_survival_validation.py                 # Kaplan-Meier, Log-Rank & Cox regression
└── tests/
    └── test_study3_survival_validation.py         # Automated verification test suite
```

---

## 🚀 Quick Start & Reproducibility

### 1. Environment Setup
```bash
git clone https://github.com/Lukecele/spatial-massive-validation-study3.git
cd spatial-massive-validation-study3
python3 -m venv venv
source venv/bin/activate
pip install numpy pandas scipy matplotlib lifelines
```

### 2. Execute Survival Pipeline
```bash
python3 src/run_survival_validation.py
```
This recomputes all Kaplan-Meier curves, log-rank tests, hazard ratios, and generates `results/figure_study3_kaplan_meier_validation.png`.

### 3. Run Automated Test Suite
```bash
python3 -m unittest discover -s tests
```
All tests verify statistical power, dataset schema integrity, and log-rank significance boundaries.

---

## 📜 Citation & Authorship

If you utilize this population dataset, filtering protocol, or survival models in your research, please cite:

```bibtex
@article{Celebrano_Study3_2026,
  author = {Celebrano, Luca},
  title = {{Population-Scale Clinical Validation of Microvascular Decoupling as an Independent Predictor of Progression-Free Survival in EGFR-Mutant Lung Adenocarcinoma (N = 854)}},
  journal = {Preprint / Independent Research},
  year = {2026},
  url = {https://github.com/Lukecele/spatial-massive-validation-study3},
  institution = {Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II}
}
```

---

## ⚖️ License

This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**. Academic institutions, public hospitals, and researchers are free to run, inspect, and apply this codebase for patient care and non-commercial scientific investigation.
