# Data-based Distributionally Robust Stochastic OPF Package

The distributionally robust stochastic optimal power flow (OPF) package is developed at the [Control, Optimization and Networks Laboratory](http://www.utdallas.edu/~tyler.summers/), [The University of Texas at Dallas](https://utdallas.edu/). This framework uses MATLAB to solve a multi-stage stochastic OPF problem based on limited information about forecast error distributions, which explicitly combines multi-stage feedback policies with any forecasting method and historical forecast error data. The objective is to determine power scheduling policies for controllable devices in a power network to balance operational cost and conditional value-at-risk (CVaR) of device and network constraint violations. These decisions include both nominal power schedules and reserve policies, which specify planned reactions to forecast errors in order to accommodate fluctuating renewable energy sources (RESs). Instead of assuming the uncertainties across the networks follow prescribed probability distributions, we consider ambiguity sets of distributions centered around a finite training dataset. By utilizing the Wasserstein metric to quantify differences between the empirical data-based distribution and the real unknown data-generating distribution, we formulate a multi-stage distributionally robust OPF problem to compute control policies that are robust to both forecast errors and sampling errors inherent in the dataset. This package includes two sub-packages, especially designed for distribution networks and transmission systems, respectively.

## Outline of README Documentation
* [Requirements](#requirements)
* [Package Organization](#package-organization)
* [Capabilities of framework](#capabilities-of-framework)
* [Getting Started](#getting-started)
* [Useful Tips](#useful-tips)
* [Useful Resources](#useful-resources)
* [Acknowledgement](#acknowledgement)
* [Notes](#notes)
* [License](#license)

## Requirements
* Operating System:
  * Windows:
    * The framework has been tested with Windows 7 and Windows 10.
  * macOS:
    * The framework has been tested with IOS 12.
* Software:
  * MATLAB:
    * This package has been tested with MATLAB R2017b and R2018a. The additional required tool boxes/solvers are:
      * CVX: Disciplined convex programming for MATLAB.
      * MATPOWER: A power system simulation and optimization tool for MATLAB
    * For more information on these packages, please visit [CVX Research](http://cvxr.com/cvx/), and [MATPOWER](http://www.pserc.cornell.edu/matpower/).

## Package Organization

The folders of the repository are:
* DistributionNetworks
  * Data-based distributionally robust stochastic OPF sub-package for distribution networks, which includes one main script and two additional functions.
    1. the main script of this sub-package is `main_distribution.m`.
    2. network admittance generation function is `form_admittance.m`.
    3. SDP relaxation solver for optimal power flow is `getVSDP.m`.
    4. Note: make sure ``CVX`` is installed properly.
* TransmissionSystems
  * Data-based distributionally robust stochastic OPF sub-package for transmission systems, which includes one main script.
    1. the main script of this sub-package is `main_transmission.m`.
    2. Note: make sure ``MATPOWER`` and ``CVX`` are installed properly.
* CVX
  * Convex optimization modeling language for MATLAB.
* MATPOWER
  * A power system simulation and optimization package for MATLAB.

## Capabilities of framework

### Overvoltage problem for distribution networks
This sub-package provides a computationally-affordable chance-contrained AC OPF framework to solve an overvoltage problem for IEEE 37-node distribution network. The proposed distributionally robust stochastic OPF methodologies mitigate overvoltages by controlling set points for renewable energy resources and energy storage devices. The set points of controllable devices are repeatedly optimized over a finite planning horizon within a MPC feedback scheme. The risk conservativeness of the voltage magnitude constraints and the out-of-sample performance robustness to sampling errors are explicitly adjustable by two scalar parameters (e.g., weight factor and Wasserstein distance).

### N-1 Security problem for transmission systems
This sub-package provides a computationally-affordable chance-contrained DC OPF framework to solve a N-1 security problem for IEEE 118-bus transmission system. The proposed distributionally robust stochastic OPF will determine scheduled power output adjustments and reserve policies of generators, which specify planned reactions to wind power forecast errors in order to accommodate fluctuating renewable energy sources. The risk conservativeness of the line flow constraints and the out-of-sample performance robustness to sampling errors are explicitly adjustable by two scalar parameters (e.g., weight factor and Wasserstein distance).

## Getting Started

The framework works in principle for any forecast method. The user needs to import forecast error data into this framework, (e.g., solar power forecast scenarios for distribution networks and wind power forecast error scenarios for transmission networks). In order to run this framework, the forecast data should be organized in the following format:

### For distribution networks
Two MATLAB data files `.mat` are required in this sub-package: solar power forecast scenarios dataset and electric power loads dataset. The solar power forecast scenarios dataset will be imported in `main_distribution.m` by the command
  ```
  load('solar.mat')
  ```
  * three matrices contained data for different usages.
  * first matrix, `PV_Real24` contains the nominal solar power forecast in kW, which defined in the dimension `PV_Real24 (Timestep,1)`.
  * second matrix, `PV_Errors_MPC_S` contains forecast error scenarios of solar power within finite time horizon in kW, which defined in the dimension `PV_Errors_MPC_S (Timestep, Horizon, DRScenarios)`.
  * third matrix, `PV_MC` contains scenarios of solar power output for Monte Carlo simulation in kW, which defined in the dimension `PV_MC (Timestep, MCScenarios)`.
  * Note that all three matrices represent aggregated power output of PV inverters one node, the framework will do the distribution over the feeder.
  
Electric power loads dataset will be imported in `main_distribution.m` by the commend
  ```
  load('load.mat')
  ```
  * two matrices contain active power and reactive power loads for every node.
  * `P_l` contains the deterministic active power loads in kW, which defined in the dimension `P_l (Nnode, Timestep)`.
  * `Q_l` contains the deterministic reactive power loads in kvar, which defined in the dimension `Q_l (Nnode, Timestep)`.

The `Timestep` indicates the number of total timesteps; `DRScenarios` is the number of total scenarios for distributionally robust optimization; `MCScenarios` is the number of total scenarios for Monte Carlo simulation; `Horizon` is the finite time horizon setting of model predictive control; `Nnode` indicates the number of nodes in the distribution network.

### For transmission systems
One MATLAB data file `.mat` is required in this sub-package: wind power forecast error scenarios dataset, which will be imported in `main_transmission.m` by the command
  ```
  load('wind.mat')
  ```
  * three vectors contain wind power forecast error scenarios (in MW) for each individual wind farm.
  * `G_error_1` contains the wind power forecast error scenarios for wind farm #1 (at bus 1), which defined in the dimension `G_error_1 (DRScenarios, 1)`.
  * `G_error_2` contains the wind power forecast error scenarios for wind farm #2 (at bus 9), which defined in the dimension `G_error_2 (DRScenarios, 1)`.
  * `G_error_3` contains the wind power forecast error scenarios for wind farm #3 (at bus 26), which defined in the dimension `G_error_3 (DRScenarios, 1)`.

The `DRScenarios` is the number of scenarios for distributionally robust optimization.

## Useful Tips

1. We recommend the `MOSEK` solver for distributionally robust optimizatoin, and `SDPT3` solver for SDP relaxation problem. Note that Implementing `MOSEK` solver may request CVX academic license, which you may possibly apply from [CVX Research](http://cvxr.com/cvx/).




## Useful Resources

Further details on the mathematical formulation and example numerical results can be found in the following papers.
1. Y. Guo, K. Baker, E. Dall'Anese, Z. Hu and T.H. Summers, "[Data-based distribubtionally robust stochastic optimal power flow, Part I: Methodologies](https://arxiv.org/abs/1804.06388)", IEEE Transactions on Power Systems, vol.34, no.2, pp.1483-1492, March 2019.
2. Y. Guo, K. Baker, E. Dall'Anese, Z. Hu and T.H. Summers, "[Data-based distribubtionally robust stochastic optimal power flow, Part II: Case studies](https://arxiv.org/abs/1804.06384)", IEEE Transactions on Power Systems, vol.34, no.2, pp.1493-1503, March 2019.
3. Y. Guo, K. Baker, E. Dall'Anese, Z. Hu and T.H. Summers, "[Stochastic optimal power flow based on data-driven distributionally robust optimization](https://arxiv.org/abs/1706.04267)", 2018 Annual American Control Conference (ACC), Milwaukee, WI, June 2018.

## Notes
The following links provide an access to publicly available data of renewable energy resources (RESs).
* [National solar radiation databased (NSRDB)](https://rredc.nrel.gov/solar/old_data/nsrdb/).
* [Global Energy Forecasting Competition 2012](https://ac.els-cdn.com/S0169207013000745/1-s2.0-S0169207013000745-main.pdf?_tid=8df2518c-5031-4865-bb7b-5c64a0d9c8f8&acdnat=1539208244_48be7ab05a631e416a83f1acb70f7d9a).
Note that the solar radiation data used in the papers above are different from the publicly available data, so the plots from this package using different forecase error data will differ from the figures presented in the above papers.

## Acknowledgement
This material is based on work supported by the National Science Foundation (NSF) under grant CNS-1566127.

## License
MIT License

Copyright (c) 2018, Yi Guo, Tyler Summers (PI),
Control, Optimization and Networks Laboratory,
Department of Mechanical Engineering,
The University of Texas at Dallas, Richardson, TX, USA.

Emails:   yi.guo2@utdallas.edu,
          tyler.summers@utdallas.edu.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
