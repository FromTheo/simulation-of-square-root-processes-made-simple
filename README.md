### Simulation of square-root processes made simple: applications to the Heston model

This repository reproduces the results of [E. Abi Jaber (2025)](https://arxiv.org/pdf/2412.11264), [L. Andersen (2006)](https://www.ressources-actuarielles.net/EXT/ISFA/1226.nsf/0/1826b88b152e65a7c12574b000347c74/$FILE/LeifAndersenHeston.pdf) and [A. Alfonsi (2008)](https://hal.science/hal-00143723v5/document).

We implement: 
- Simulation of the iVi scheme (see `ivi.py`), along with other discretizations, including Andersen's Truncatured Gaussian (TG) and Quadratic Exponential (QE) schemes, and Alfonsi's second-order scheme. 
- Martingale corrections for the TG and QE schemes (see `martingale_correction.ipynb`). 
- Empirical comparison across several experiments: pricing of $q$-volatility swaps, European options pricing (with reference prices obtained via Fourier inversion), and implied volatility smiles under market revelant regimes (in particular, the Feller condition is violated) (see `demo.ipynb`). 
- Since the law of $V_t|V_s$ $s\leq t$ is proportional to a noncentral chi-square distribution, we also evaluate each scheme's ability to reproduce the conditional CDF then $s=0$.  

In addition, we provide a technical note (see `simulate_square_root_processes.pdf`) detailing the construction of the TG(-M), QE(-M), and Alfonsi second-order schemes. 

### Examples of illustrations 

![swaps](assets/swaps_case_1.png)

![prices](assets/prices_case_1.png)

![iv](assets/iv_case_1.png)

![distribution](assets/distribution_V.png)

### Disclaimer 
Source code is available upon request. Please contact me directly. 