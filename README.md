# CancerABC

## Tumor Subclonal Selection and Evolutionary Dynamics

This repository contains two Python codes for analyzing tumor growth and subclonal selection dynamics:

1. **TumorSimul3D_SMALT_simulate_monoclonal.py**: This code simulates neoplastic growth and subclonal selection in a 3D agent-based tumor model.
2. **CancerABC.py**: This code estimates evolutionary dynamics for monoclonal tumors using Approximate Bayesian Computation (ABC).

### TumorSimul3D_SMALT_simulate_monoclonal.py

This code simulates tumor growth and mutation accumulation after a polyclonal to monoclonal transition. It incorporates "deme" subpopulations to represent glandular structures in intestinal tumors. Each deme represents a cube in a 3D lattice with interacting neighboring subpopulations.

**Running the simulation:**

```
python TumorSimul3D_SMALT_simulate_monoclonal.py <replicate_name> <num_hotspots> <selection_coefficient>
```

* `<replicate_name>`: Name for the simulation replicate (integer)
* `<num_hotspots>`: Number of hotspot mutations to simulate (integer)
* `<selection_coefficient>`: Selection coefficient for advantageous mutations (float, e.g. 0 for neutral, 0.01, 0.02)

**Description:**

This code is an adaptation of a previously established 3D agent-based tumor modeling framework ([reference 2] in code). It simulates:

* Spatial tumor growth within a 3D lattice with interacting subpopulations (demes).
* Cell birth and death processes within each deme.
* Accumulation of mutations, including 3kb barcode mutations and cancer driver mutations.

### CancerABC.py

This code employs Approximate Bayesian Computation (ABC) to estimate the evolutionary dynamics of monoclonal tumors. It utilizes the `TumorSimul3D_SMALT_simulate_monoclonal.py` code for simulations.

**Running ABC inference:**

```
python CancerABC.py --inputRdata example.phy --nprocs 80 --nsim 500 --npopulation 10
```

* `--inputRdata`: Path to the input data file (e.g., example.phy).
* `--nprocs`: Number of processes for parallelization (e.g., 80).
* `--nsim`: Number of simulations per ABC iteration (e.g., 500).
* `--npopulation`: The number of parameter draws to accept in each SMC iteration (e.g., 10).

**Description:**

This code estimates evolutionary dynamics by:

* Defining key parameters:
    * `u`: Mutation rate per cell in the whole 3kb barcode.
    * `b`: Division probability per cell generation.
    * `w`: Selection coefficient benefit for advantageous mutations (converted to `s` later).

**References:**

* **Lu, Zhaolian, Shanlan Mo, Duo Xie, Xiangwei Zhai, Shanjun Deng, Kantian Zhou, Kun Wang, Xueling Kang, Hao Zhang, Juanzhen Tong, Liangzhen Hou, Huijuan Hu, Xuefei Li, Da Zhou, Leo Tsz On Lee, Li Liu, Yaxi Zhu, Jing Yu, Ping Lan, Jiguang Wang, Zhen He, Xionglei He and Zheng Hu. 2024. Polyclonal-to-Monoclonal Transition in Colorectal Precancerous Evolution. *Nature* 636(8041):233–40. doi: 10.1038/s41586-024-08133-1. [https://www.nature.com/articles/s41586-024-08133-1](https://www.nature.com/articles/s41586-024-08133-1)**

Contact: For questions or suggestions, please open an issue on this GitHub repository or reach out directly via email (duo.xie at siat.ac.cn).
