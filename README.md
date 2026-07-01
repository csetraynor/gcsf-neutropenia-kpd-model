# Kinetic-Pharmacodynamic Stan Models for Docetaxel-Induced Neutropenia

This repository contains two Bayesian kinetic-pharmacodynamic (KPD) models, implemented in [Stan](https://mc-stan.org/) with [Torsten](https://github.com/metrumresearchgroup/Torsten) extensions, for describing absolute neutrophil count (ANC) dynamics following docetaxel chemotherapy:

- **KPD model without G-CSF** — a Friberg-style transit-compartment model with virtual kinetics for docetaxel-driven myelosuppression.
- **KPD model with G-CSF** — extends the above to incorporate granulocyte colony-stimulating factor (G-CSF) stimulation of neutrophil proliferation and maturation.

These models are associated with a poster presented at **PAGE 34 (2026)**.

---

## Repository Contents

```
models/
├── modKPD.stan        # KPD model without G-CSF
└── modKPDGCSF.stan    # KPD model with G-CSF stimulation
LICENSE                # Apache License 2.0
gcsf-neutropenia-kpd-model.Rproj   # RStudio project file
README.md
```

The repository currently contains the core Stan model definitions. Data preparation, fitting scripts, posterior updating pipelines, and downstream machine learning components described in the associated abstract are not included at present.

---

## Scientific Context

Neutropenia is a common dose-limiting toxicity of cytotoxic chemotherapy. These models describe individual ANC trajectories following docetaxel treatment using a **Friberg-style semi-mechanistic framework** with:

- A proliferative cell compartment with feedback regulation
- Three transit compartments representing neutrophil maturation
- A circulating neutrophil compartment
- A virtual kinetic-pharmacodynamic (KPD) component representing docetaxel-driven proliferation inhibition

The **G-CSF-augmented model** (`modKPDGCSF.stan`) additionally includes:

- G-CSF subcutaneous absorption and elimination kinetics
- Stimulatory effects on both proliferation rate and maturation (transit) rate via saturable (Emax) functions

Both models are fitted using Hamiltonian Monte Carlo (NUTS) via Torsten Stan, with literature-informed lognormal population priors and inter-individual variability on key parameters. Generated quantities include pointwise log-likelihood for model comparison (LOO-CV / WAIC).

In the broader study, the G-CSF model demonstrated substantially better out-of-sample fit (ELPD LOO −4808 vs −5629), and model-informed ANC features improved early prediction of docetaxel-induced neutropenia when combined with machine learning methods.

---

## Associated Abstract

**KPD-ML prediction of docetaxel-induced neutropenia in large RWE datasets**

- **Conference:** PAGE 34 (2026) Abstr 12119
- **Category:** Poster — Clinical Applications
- **Authors:** Carlos Serra Traynor, Marija Kekic, David Boulton, Diansong Zhou
- **Affiliations:** Clinical Pharmacology & Quantitative Pharmacology and Predictive AI & Data, Clinical Pharmacology & Safety Sciences, R&D BioPharmaceuticals, AstraZeneca
- **Abstract:** [https://www.page-meeting.org/Abstracts/kpd-ml-prediction-of-docetaxel-induced-neutropenia-in-large-rwe-datasets/](https://www.page-meeting.org/Abstracts/kpd-ml-prediction-of-docetaxel-induced-neutropenia-in-large-rwe-datasets/)

---

## How to Cite

If you use or reference these models, please cite the associated abstract:

> Serra Traynor C, Kekic M, Boulton D, Zhou D. KPD-ML prediction of docetaxel-induced neutropenia in large RWE datasets. PAGE 34 (2026) Abstr 12119. [https://www.page-meeting.org/Abstracts/kpd-ml-prediction-of-docetaxel-induced-neutropenia-in-large-rwe-datasets/](https://www.page-meeting.org/Abstracts/kpd-ml-prediction-of-docetaxel-induced-neutropenia-in-large-rwe-datasets/)

A formal software DOI is not currently available for this repository.

---

## Usage

### Requirements

- [Stan](https://mc-stan.org/) (≥ 2.26 recommended)
- [Torsten](https://github.com/metrumresearchgroup/Torsten) — required for the `pmx_solve_group_rk45` ODE solver used in both models
- An interface to Stan such as [CmdStan](https://mc-stan.org/users/interfaces/cmdstan), [CmdStanR](https://mc-stan.org/cmdstanr/), or [CmdStanPy](https://cmdstanpy.readthedocs.io/)

### Notes

- The models use Torsten-specific ODE solver functions (`pmx_solve_group_rk45`) and require a Torsten-enabled Stan installation for compilation.
- Input data must be formatted according to the data blocks defined in each Stan file (NONMEM-style event records with `time`, `amt`, `evid`, `cmt`, `rate`, `ii`, `addl`).
- Population priors are supplied as data inputs, allowing flexible prior specification without recompilation.
- No runnable scripts or example data are currently included in this repository. Compilation and execution depend on the user's Stan/Torsten environment and appropriately formatted input data.

---

## Data Availability

Patient-level data used in the associated study are from the Flatiron Health database and are **not included** in this repository due to data use agreements and patient privacy considerations.

---

## Scope and Limitations

This repository currently contains the **core Stan model definitions** only. The following components described in the associated abstract are not included at present:

- Data extraction and preprocessing pipelines
- Model fitting and diagnostic scripts
- Prior-posterior updating workflow for validation
- Stochastic gates (STG) covariate screening
- XGBoost prediction pipeline
- Simulation and visual predictive check code

The models are provided for transparency and reproducibility of the KPD model structure. Users interested in the full analytical workflow should refer to the abstract and contact the authors.

---

## License

This project is licensed under the **Apache License 2.0** — see the [LICENSE](LICENSE) file for details.
