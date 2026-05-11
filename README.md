# Exact identification of unknown unitary processes

[![arXiv](https://img.shields.io/badge/arXiv-2605.04981-b31b1b.svg)](https://arxiv.org/abs/2605.04981)

Code and numerical calculations for the paper:\
**Exact identification of unknown unitary processes** – [arXiv:2605.04981](https://arxiv.org/abs/2605.04981)
## Overview

This repository contains the computational framework used to evaluate protocols for identifying anomalous unitary processes within the framework of quantum testers. In our paper, we address the problem of identifying $k$ malfunctioning devices (which apply an unknown unitary action) among a series of $n$ devices intended to apply a fixed, known unitary transformation.

We study the **zero-error protocol** for identifying these faulty devices under the assumption of complete ignorance regarding the anomalous unitary, considering the Haar average of the success probability as our figure of merit. The code here supports the theoretical derivations of the optimal success probability for the two-anomaly scenario, as well as numerical computations for small system sizes.

## Code

The problem of zero-error discrimination of quantum channels via quantum testers is structured as an optimization problem. In the single-anomaly case, local symmetries of the Choi operators allow the optimization to be restricted to a reduced family of testers, leading to a closed-form expression. 

However, moving to two or more anomalies introduces a substantial increase in complexity. To compute the optimal protocol and its associated success probability, we analyze the structural constraints of the symmetrized testers, and we rely on **Semidefinite Programming (SDP)**.

This repository provides the code to:
1. Construct the operators in the algebra of partially transposed permutations.
2. Run the Semidefinite Programs (SDPs) to find the maximal success probability $P_s$ for small system sizes.

## Files Included

- `Dual_operator.nb` – **Mathematica notebook** containing the construction of the operator in Eq. (29) in the manuscript. This file uses the following two files from the code associated with [arXiv:2310.02252](https://arxiv.org/abs/2310.02252v1), which can be found in [github:walledbrauer-gtbasis](https://github.com/dgrinko/walledbrauer-gtbasis):
    - `Mixed Schur Tools.wl` 
    - `walled_brauer_gt_basis.wl` 
- Two **MATLAB** scripts for computing the success probability of the parallel and sequential strategies, both of them using [QETLAB](https://qetlab.com)
    - `ZE_par.m` – for **parallel** strategies
    - `ZE_seq.m` – for **sequential** strategies uses the function in:
        - `odskok.m`

## Authors

- **Santiago Llorens** (UAB and UG)
- **Arnau Diebra** (UAB)
- **Michal Sedlák** (RCQI and MUNI)
- **Ramon Muñoz-Tapia** (UAB)

## Citation

If you find this code useful in your research, please consider citing our paper:\
**Exact identification of unknown unitary processes** – [arXiv:2605.04981](https://arxiv.org/abs/2605.04981)

