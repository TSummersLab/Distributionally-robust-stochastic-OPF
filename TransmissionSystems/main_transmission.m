%% Data-based distributionally robust stochastic optimal power flow (OPF)
% N-1 Security problem for IEEE 118-bus transmission system
%
% Latest update: 2018-Oct-12
%
% Copyright (c) 2018-2019, Yi Guo, Tyler Summers (PI)
%
% Control, Optimization and Networks Laboratory
% Department of Mechanical Engineering
% The University of Texas at Dallas, Richardson, TX, USA.
%
% Emails:   yi.guo2@utdallas.edu
%           tyler.summers@utdallas.edu
%
% This sub-package provides a computationally-affordable chance-contrained 
% DC OPF framework to solve a N-1 security problem for IEEE 118-bus 
% transmission system. The proposed distributionally robust stochastic 
% OPF will determine scheduled power output adjustments and reserve 
% policies of generators, which specify planned reactions to wind power 
% forecast errors in order to accommodate fluctuating renewable energy 
% sources. The risk conservativeness of the line flow constraints and 
% the out-of-sample performance robustness to sampling errors are 
% explicitly adjustable by two scalar parameters (e.g., weight factor 
% and Wasserstein distance).
%
% In general, the problem formulated in this script could be treated as a
% block of a multi-stage distributionally robust stochastic OPF, which
% also works for any forecasting methods in any resolution format.
% please refer to "README.md" for how to import necessary data: wind 
% forecast errors.
%
% The mathematical formulation and results discussions can be found in the
% following two references.
%
% [1]Y.Guo, K.Baker, E.Dall'Anese, Z.Hu and T.Summers, "Data-based
% distribubtionally robust stochastic optimal power flow, Part
% I: Methodologies", available on arXiv.org.
%
% [2]Y.Guo, K.Baker, E.Dall'Anese, Z.Hu and T.Summers, "Data-based
% distribubtionally robust stochastic optimal power flow, Part
% II: Case studies", available on arXiv.org.
%
% This material is based on works supported by the National Science
% Foundation (NSF) under grant CNS-1566127. (PI: Dr. Tyler H. Summers)
%
% Questions and comments are welcome.
%
% MIT License
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.



clc; clear all; close all;
warning off


disp(['Running Data-based distributionally robust stochastic optimal power flow (OPF)'])
disp(['N-1 Security problem for IEEE 118 buses transmission system'])
disp(['Latest update: 2018-Oct-12'])

disp(['..................................................................................'])

disp(['Copyright (c) 2018-2019, Yi Guo, Tyler Summers (PI)'])
disp(['Control, Optimization and Network Laboratory'])
disp(['Department of Mechanical Engineering'])
disp(['The University of Texas at Dallas, Richardson, TX, USA.'])

%% load input data
disp(['..................................................................................'])
load('wind.mat')
disp(['wind power forecasting errors imported sucessfully ......'])
disp(['..................................................................................'])

tic
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% program and system setup  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setup system parameters
%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
    1	2	51	27	0	0	1	0.955	10.67	138	1	1.06	0.94;
    2	1	20	9	0	0	1	0.971	11.22	138	1	1.06	0.94;
    3	1	39	10	0	0	1	0.968	11.56	138	1	1.06	0.94;
    4	2	39	12	0	0	1	0.998	15.28	138	1	1.06	0.94;
    5	1	0	0	0	-40	1	1.002	15.73	138	1	1.06	0.94;
    6	2	52	22	0	0	1	0.99	13	    138	1	1.06	0.94;
    7	1	19	2	0	0	1	0.989	12.56	138	1	1.06	0.94;
    8	2	28	0	0	0	1	1.015	20.77	345	1	1.06	0.94;
    9	1	0	0	0	0	1	1.043	28.02	345	1	1.06	0.94;
    10	2	0	0	0	0	1	1.05	35.61	345	1	1.06	0.94;
    11	1	70	23	0	0	1	0.985	12.72	138	1	1.06	0.94;
    12	2	47	10	0	0	1	0.99	12.2	138	1	1.06	0.94;
    13	1	34	16	0	0	1	0.968	11.35	138	1	1.06	0.94;
    14	1	14	1	0	0	1	0.984	11.5	138	1	1.06	0.94;
    15	2	90	30	0	0	1	0.97	11.23	138	1	1.06	0.94;
    16	1	25	10	0	0	1	0.984	11.91	138	1	1.06	0.94;
    17	1	11	3	0	0	1	0.995	13.74	138	1	1.06	0.94;
    18	2	60	34	0	0	1	0.973	11.53	138	1	1.06	0.94;
    19	2	45	25	0	0	1	0.963	11.05	138	1	1.06	0.94;
    20	1	18	3	0	0	1	0.958	11.93	138	1	1.06	0.94;
    21	1	14	8	0	0	1	0.959	13.52	138	1	1.06	0.94;
    22	1	10	5	0	0	1	0.97	16.08	138	1	1.06	0.94;
    23	1	7	3	0	0	1	1	    21	    138	1	1.06	0.94;
    24	2	13	0	0	0	1	0.992	20.89	138	1	1.06	0.94;
    25	2	0	0	0	0	1	1.05	27.93	138	1	1.06	0.94;
    26	2	0	0	0	0	1	1.015	29.71	345	1	1.06	0.94;
    27	2	71	13	0	0	1	0.968	15.35	138	1	1.06	0.94;
    28	1	17	7	0	0	1	0.962	13.62	138	1	1.06	0.94;
    29	1	24	4	0	0	1	0.963	12.63	138	1	1.06	0.94;
    30	1	0	0	0	0	1	0.968	18.79	345	1	1.06	0.94;
    31	2	43	27	0	0	1	0.967	12.75	138	1	1.06	0.94;
    32	2	59	23	0	0	1	0.964	14.8	138	1	1.06	0.94;
    33	1	23	9	0	0	1	0.972	10.63	138	1	1.06	0.94;
    34	2	59	26	0	14	1	0.986	11.3	138	1	1.06	0.94;
    35	1	33	9	0	0	1	0.981	10.87	138	1	1.06	0.94;
    36	2	31	17	0	0	1	0.98	10.87	138	1	1.06	0.94;
    37	1	0	0	0	-25	1	0.992	11.77	138	1	1.06	0.94;
    38	1	0	0	0	0	1	0.962	16.91	345	1	1.06	0.94;
    39	1	27	11	0	0	1	0.97	8.41	138	1	1.06	0.94;
    40	2	66	23	0	0	1	0.97	7.35	138	1	1.06	0.94;
    41	1	37	10	0	0	1	0.967	6.92	138	1	1.06	0.94;
    42	2	96	23	0	0	1	0.985	8.53	138	1	1.06	0.94;
    43	1	18	7	0	0	1	0.978	11.28	138	1	1.06	0.94;
    44	1	16	8	0	10	1	0.985	13.82	138	1	1.06	0.94;
    45	1	53	22	0	10	1	0.987	15.67	138	1	1.06	0.94;
    46	2	28	10	0	10	1	1.005	18.49	138	1	1.06	0.94;
    47	1	34	0	0	0	1	1.017	20.73	138	1	1.06	0.94;
    48	1	20	11	0	15	1	1.021	19.93	138	1	1.06	0.94;
    49	2	87	30	0	0	1	1.025	20.94	138	1	1.06	0.94;
    50	1	17	4	0	0	1	1.001	18.9	138	1	1.06	0.94;
    51	1	17	8	0	0	1	0.967	16.28	138	1	1.06	0.94;
    52	1	18	5	0	0	1	0.957	15.32	138	1	1.06	0.94;
    53	1	23	11	0	0	1	0.946	14.35	138	1	1.06	0.94;
    54	2	113	32	0	0	1	0.955	15.26	138	1	1.06	0.94;
    55	2	63	22	0	0	1	0.952	14.97	138	1	1.06	0.94;
    56	2	84	18	0	0	1	0.954	15.16	138	1	1.06	0.94;
    57	1	12	3	0	0	1	0.971	16.36	138	1	1.06	0.94;
    58	1	12	3	0	0	1	0.959	15.51	138	1	1.06	0.94;
    59	2	277	113	0	0	1	0.985	19.37	138	1	1.06	0.94;
    60	1	78	3	0	0	1	0.993	23.15	138	1	1.06	0.94;
    61	2	0	0	0	0	1	0.995	24.04	138	1	1.06	0.94;
    62	2	77	14	0	0	1	0.998	23.43	138	1	1.06	0.94;
    63	1	0	0	0	0	1	0.969	22.75	345	1	1.06	0.94;
    64	1	0	0	0	0	1	0.984	24.52	345	1	1.06	0.94;
    65	2	0	0	0	0	1	1.005	27.65	345	1	1.06	0.94;
    66	2	39	18	0	0	1	1.05	27.48	138	1	1.06	0.94;
    67	1	28	7	0	0	1	1.02	24.84	138	1	1.06	0.94;
    68	1	0	0	0	0	1	1.003	27.55	345	1	1.06	0.94;
    69	3	0	0	0	0	1	1.035	30	    138	1	1.06	0.94;
    70	2	66	20	0	0	1	0.984	22.58	138	1	1.06	0.94;
    71	1	0	0	0	0	1	0.987	22.15	138	1	1.06	0.94;
    72	2	12	0	0	0	1	0.98	20.98	138	1	1.06	0.94;
    73	2	6	0	0	0	1	0.991	21.94	138	1	1.06	0.94;
    74	2	68	27	0	12	1	0.958	21.64	138	1	1.06	0.94;
    75	1	47	11	0	0	1	0.967	22.91	138	1	1.06	0.94;
    76	2	68	36	0	0	1	0.943	21.77	138	1	1.06	0.94;
    77	2	61	28	0	0	1	1.006	26.72	138	1	1.06	0.94;
    78	1	71	26	0	0	1	1.003	26.42	138	1	1.06	0.94;
    79	1	39	32	0	20	1	1.009	26.72	138	1	1.06	0.94;
    80	2	130	26	0	0	1	1.04	28.96	138	1	1.06	0.94;
    81	1	0	0	0	0	1	0.997	28.1	345	1	1.06	0.94;
    82	1	54	27	0	20	1	0.989	27.24	138	1	1.06	0.94;
    83	1	20	10	0	10	1	0.985	28.42	138	1	1.06	0.94;
    84	1	11	7	0	0	1	0.98	30.95	138	1	1.06	0.94;
    85	2	24	15	0	0	1	0.985	32.51	138	1	1.06	0.94;
    86	1	21	10	0	0	1	0.987	31.14	138	1	1.06	0.94;
    87	2	0	0	0	0	1	1.015	31.4	161	1	1.06	0.94;
    88	1	48	10	0	0	1	0.987	35.64	138	1	1.06	0.94;
    89	2	0	0	0	0	1	1.005	39.69	138	1	1.06	0.94;
    90	2	163	42	0	0	1	0.985	33.29	138	1	1.06	0.94;
    91	2	10	0	0	0	1	0.98	33.31	138	1	1.06	0.94;
    92	2	65	10	0	0	1	0.993	33.8	138	1	1.06	0.94;
    93	1	12	7	0	0	1	0.987	30.79	138	1	1.06	0.94;
    94	1	30	16	0	0	1	0.991	28.64	138	1	1.06	0.94;
    95	1	42	31	0	0	1	0.981	27.67	138	1	1.06	0.94;
    96	1	38	15	0	0	1	0.993	27.51	138	1	1.06	0.94;
    97	1	15	9	0	0	1	1.011	27.88	138	1	1.06	0.94;
    98	1	34	8	0	0	1	1.024	27.4	138	1	1.06	0.94;
    99	2	42	0	0	0	1	1.01	27.04	138	1	1.06	0.94;
    100	2	37	18	0	0	1	1.017	28.03	138	1	1.06	0.94;
    101	1	22	15	0	0	1	0.993	29.61	138	1	1.06	0.94;
    102	1	5	3	0	0	1	0.991	32.3	138	1	1.06	0.94;
    103	2	23	16	0	0	1	1.001	24.44	138	1	1.06	0.94;
    104	2	38	25	0	0	1	0.971	21.69	138	1	1.06	0.94;
    105	2	31	26	0	20	1	0.965	20.57	138	1	1.06	0.94;
    106	1	43	16	0	0	1	0.962	20.32	138	1	1.06	0.94;
    107	2	50	12	0	6	1	0.952	17.53	138	1	1.06	0.94;
    108	1	2	1	0	0	1	0.967	19.38	138	1	1.06	0.94;
    109	1	8	3	0	0	1	0.967	18.93	138	1	1.06	0.94;
    110	2	39	30	0	6	1	0.973	18.09	138	1	1.06	0.94;
    111	2	0	0	0	0	1	0.98	19.74	138	1	1.06	0.94;
    112	2	68	13	0	0	1	0.975	14.99	138	1	1.06	0.94;
    113	2	6	0	0	0	1	0.993	13.74	138	1	1.06	0.94;
    114	1	8	3	0	0	1	0.96	14.46	138	1	1.06	0.94;
    115	1	22	7	0	0	1	0.96	14.46	138	1	1.06	0.94;
    116	2	184	0	0	0	1	1.005	27.12	138	1	1.06	0.94;
    117	1	20	8	0	0	1	0.974	10.67	138	1	1.06	0.94;
    118	1	33	15	0	0	1	0.949	21.92	138	1	1.06	0.94;
    ];
%% generator data
% several generators are removed from orginal setting
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
    %     1	0	0	15	-5	0.955	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    4	0	0	300	-300	0.998	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    6	0	0	50	-13	0.99	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    8	0	0	300	-300	1.015	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    10	450	0	200	-147	1.05	100	1	550	0	0	0	0	0	0	0	0	0	0	0	0;
    12	85	0	120	-35	0.99	100	1	185	0	0	0	0	0	0	0	0	0	0	0	0;
    15	0	0	30	-10	0.97	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    18	0	0	50	-16	0.973	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    19	0	0	24	-8	0.962	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    24	0	0	300	-300	0.992	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    25	220	0	140	-47	1.05	100	1	320	0	0	0	0	0	0	0	0	0	0	0	0;
    %     26	314	0	1000	-1000	1.015	100	1	414	0	0	0	0	0	0	0	0	0	0	0	0;
    27	0	0	300	-300	0.968	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    31	7	0	300	-300	0.967	100	1	107	0	0	0	0	0	0	0	0	0	0	0	0;
    32	0	0	42	-14	0.963	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    34	0	0	24	-8	0.984	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    36	0	0	24	-8	0.98	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    40	0	0	300	-300	0.97	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    42	0	0	300	-300	0.985	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    46	19	0	100	-100	1.005	100	1	119	0	0	0	0	0	0	0	0	0	0	0	0;
    49	204	0	210	-85	1.025	100	1	304	0	0	0	0	0	0	0	0	0	0	0	0;
    54	48	0	300	-300	0.955	100	1	148	0	0	0	0	0	0	0	0	0	0	0	0;
    55	0	0	23	-8	0.952	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    56	0	0	15	-8	0.954	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    59	155	0	180	-60	0.985	100	1	255	0	0	0	0	0	0	0	0	0	0	0	0;
    61	160	0	300	-100	0.995	100	1	260	0	0	0	0	0	0	0	0	0	0	0	0;
    62	0	0	20	-20	0.998	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    65	391	0	200	-67	1.005	100	1	491	0	0	0	0	0	0	0	0	0	0	0	0;
    66	392	0	200	-67	1.05	100	1	492	0	0	0	0	0	0	0	0	0	0	0	0;
    69	516.4	0	300	-300	1.035	100	1	805.2	0	0	0	0	0	0	0	0	0	0	0	0;
    70	0	0	32	-10	0.984	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    72	0	0	100	-100	0.98	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    73	0	0	100	-100	0.991	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    74	0	0	9	-6	0.958	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    76	0	0	23	-8	0.943	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    77	0	0	70	-20	1.006	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    80	477	0	280	-165	1.04	100	1	577	0	0	0	0	0	0	0	0	0	0	0	0;
    85	0	0	23	-8	0.985	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    87	4	0	1000	-100	1.015	100	1	104	0	0	0	0	0	0	0	0	0	0	0	0;
    89	607	0	300	-210	1.005	100	1	707	0	0	0	0	0	0	0	0	0	0	0	0;
    90	0	0	300	-300	0.985	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    91	0	0	100	-100	0.98	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    92	0	0	9	-3	0.99	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    99	0	0	100	-100	1.01	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    100	252	0	155	-50	1.017	100	1	352	0	0	0	0	0	0	0	0	0	0	0	0;
    103	40	0	40	-15	1.01	100	1	140	0	0	0	0	0	0	0	0	0	0	0	0;
    104	0	0	23	-8	0.971	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    105	0	0	23	-8	0.965	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    107	0	0	200	-200	0.952	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    110	0	0	23	-8	0.973	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    111	36	0	1000	-100	0.98	100	1	136	0	0	0	0	0	0	0	0	0	0	0	0;
    112	0	0	1000	-100	0.975	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    113	0	0	200	-100	0.993	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    116	0	0	1000	-1000	1.005	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    ];
%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
    1	2	0.0303	0.0999	0.0254	0	0	0	0	0	1	-360	360;
    1	3	0.0129	0.0424	0.01082	0	0	0	0	0	1	-360	360;
    4	5	0.00176	0.00798	0.0021	0	0	0	0	0	1	-360	360;
    3	5	0.0241	0.108	0.0284	0	0	0	0	0	1	-360	360;
    5	6	0.0119	0.054	0.01426	0	0	0	0	0	1	-360	360;
    6	7	0.00459	0.0208	0.0055	0	0	0	0	0	1	-360	360;
    8	9	0.00244	0.0305	1.162	0	0	0	0	0	1	-360	360;
    8	5	0	0.0267	0	0	0	0	0.985	0	1	-360	360;
    9	10	0.00258	0.0322	1.23	0	0	0	0	0	1	-360	360;
    4	11	0.0209	0.0688	0.01748	0	0	0	0	0	1	-360	360;
    5	11	0.0203	0.0682	0.01738	0	0	0	0	0	1	-360	360;
    11	12	0.00595	0.0196	0.00502	0	0	0	0	0	1	-360	360;
    2	12	0.0187	0.0616	0.01572	0	0	0	0	0	1	-360	360;
    3	12	0.0484	0.16	0.0406	0	0	0	0	0	1	-360	360;
    7	12	0.00862	0.034	0.00874	0	0	0	0	0	1	-360	360;
    11	13	0.02225	0.0731	0.01876	0	0	0	0	0	1	-360	360;
    12	14	0.0215	0.0707	0.01816	0	0	0	0	0	1	-360	360;
    13	15	0.0744	0.2444	0.06268	0	0	0	0	0	1	-360	360;
    14	15	0.0595	0.195	0.0502	0	0	0	0	0	1	-360	360;
    12	16	0.0212	0.0834	0.0214	0	0	0	0	0	1	-360	360;
    15	17	0.0132	0.0437	0.0444	0	0	0	0	0	1	-360	360;
    16	17	0.0454	0.1801	0.0466	0	0	0	0	0	1	-360	360;
    17	18	0.0123	0.0505	0.01298	0	0	0	0	0	1	-360	360;
    18	19	0.01119	0.0493	0.01142	0	0	0	0	0	1	-360	360;
    19	20	0.0252	0.117	0.0298	0	0	0	0	0	1	-360	360;
    15	19	0.012	0.0394	0.0101	0	0	0	0	0	1	-360	360;
    20	21	0.0183	0.0849	0.0216	0	0	0	0	0	1	-360	360;
    21	22	0.0209	0.097	0.0246	0	0	0	0	0	1	-360	360;
    22	23	0.0342	0.159	0.0404	0	0	0	0	0	1	-360	360;
    23	24	0.0135	0.0492	0.0498	0	0	0	0	0	1	-360	360;
    23	25	0.0156	0.08	0.0864	0	0	0	0	0	1	-360	360;
    26	25	0	0.0382	0	0	0	0	0.96	0	1	-360	360;
    25	27	0.0318	0.163	0.1764	0	0	0	0	0	1	-360	360;
    27	28	0.01913	0.0855	0.0216	0	0	0	0	0	1	-360	360;
    28	29	0.0237	0.0943	0.0238	0	0	0	0	0	1	-360	360;
    30	17	0	0.0388	0	0	0	0	0.96	0	1	-360	360;
    8	30	0.00431	0.0504	0.514	0	0	0	0	0	1	-360	360;
    26	30	0.00799	0.086	0.908	0	0	0	0	0	1	-360	360;
    17	31	0.0474	0.1563	0.0399	0	0	0	0	0	1	-360	360;
    29	31	0.0108	0.0331	0.0083	0	0	0	0	0	1	-360	360;
    23	32	0.0317	0.1153	0.1173	0	0	0	0	0	1	-360	360;
    31	32	0.0298	0.0985	0.0251	0	0	0	0	0	1	-360	360;
    27	32	0.0229	0.0755	0.01926	0	0	0	0	0	1	-360	360;
    15	33	0.038	0.1244	0.03194	0	0	0	0	0	1	-360	360;
    19	34	0.0752	0.247	0.0632	0	0	0	0	0	1	-360	360;
    35	36	0.00224	0.0102	0.00268	0	0	0	0	0	1	-360	360;
    35	37	0.011	0.0497	0.01318	0	0	0	0	0	1	-360	360;
    33	37	0.0415	0.142	0.0366	0	0	0	0	0	1	-360	360;
    34	36	0.00871	0.0268	0.00568	0	0	0	0	0	1	-360	360;
    34	37	0.00256	0.0094	0.00984	0	0	0	0	0	1	-360	360;
    38	37	0	0.0375	0	0	0	0	0.935	0	1	-360	360;
    37	39	0.0321	0.106	0.027	0	0	0	0	0	1	-360	360;
    37	40	0.0593	0.168	0.042	0	0	0	0	0	1	-360	360;
    30	38	0.00464	0.054	0.422	0	0	0	0	0	1	-360	360;
    39	40	0.0184	0.0605	0.01552	0	0	0	0	0	1	-360	360;
    40	41	0.0145	0.0487	0.01222	0	0	0	0	0	1	-360	360;
    40	42	0.0555	0.183	0.0466	0	0	0	0	0	1	-360	360;
    41	42	0.041	0.135	0.0344	0	0	0	0	0	1	-360	360;
    43	44	0.0608	0.2454	0.06068	0	0	0	0	0	1	-360	360;
    34	43	0.0413	0.1681	0.04226	0	0	0	0	0	1	-360	360;
    44	45	0.0224	0.0901	0.0224	0	0	0	0	0	1	-360	360;
    45	46	0.04	0.1356	0.0332	0	0	0	0	0	1	-360	360;
    46	47	0.038	0.127	0.0316	0	0	0	0	0	1	-360	360;
    46	48	0.0601	0.189	0.0472	0	0	0	0	0	1	-360	360;
    47	49	0.0191	0.0625	0.01604	0	0	0	0	0	1	-360	360;
    42	49	0.0715	0.323	0.086	0	0	0	0	0	1	-360	360;
    42	49	0.0715	0.323	0.086	0	0	0	0	0	1	-360	360;
    45	49	0.0684	0.186	0.0444	0	0	0	0	0	1	-360	360;
    48	49	0.0179	0.0505	0.01258	0	0	0	0	0	1	-360	360;
    49	50	0.0267	0.0752	0.01874	0	0	0	0	0	1	-360	360;
    49	51	0.0486	0.137	0.0342	0	0	0	0	0	1	-360	360;
    51	52	0.0203	0.0588	0.01396	0	0	0	0	0	1	-360	360;
    52	53	0.0405	0.1635	0.04058	0	0	0	0	0	1	-360	360;
    53	54	0.0263	0.122	0.031	0	0	0	0	0	1	-360	360;
    49	54	0.073	0.289	0.0738	0	0	0	0	0	1	-360	360;
    49	54	0.0869	0.291	0.073	0	0	0	0	0	1	-360	360;
    54	55	0.0169	0.0707	0.0202	0	0	0	0	0	1	-360	360;
    54	56	0.00275	0.00955	0.00732	0	0	0	0	0	1	-360	360;
    55	56	0.00488	0.0151	0.00374	0	0	0	0	0	1	-360	360;
    56	57	0.0343	0.0966	0.0242	0	0	0	0	0	1	-360	360;
    50	57	0.0474	0.134	0.0332	0	0	0	0	0	1	-360	360;
    56	58	0.0343	0.0966	0.0242	0	0	0	0	0	1	-360	360;
    51	58	0.0255	0.0719	0.01788	0	0	0	0	0	1	-360	360;
    54	59	0.0503	0.2293	0.0598	0	0	0	0	0	1	-360	360;
    56	59	0.0825	0.251	0.0569	0	0	0	0	0	1	-360	360;
    56	59	0.0803	0.239	0.0536	0	0	0	0	0	1	-360	360;
    55	59	0.04739	0.2158	0.05646	0	0	0	0	0	1	-360	360;
    59	60	0.0317	0.145	0.0376	0	0	0	0	0	1	-360	360;
    59	61	0.0328	0.15	0.0388	0	0	0	0	0	1	-360	360;
    60	61	0.00264	0.0135	0.01456	0	0	0	0	0	1	-360	360;
    60	62	0.0123	0.0561	0.01468	0	0	0	0	0	1	-360	360;
    61	62	0.00824	0.0376	0.0098	0	0	0	0	0	1	-360	360;
    63	59	0	0.0386	0	0	0	0	0.96	0	1	-360	360;
    63	64	0.00172	0.02	0.216	0	0	0	0	0	1	-360	360;
    64	61	0	0.0268	0	0	0	0	0.985	0	1	-360	360;
    38	65	0.00901	0.0986	1.046	0	0	0	0	0	1	-360	360;
    64	65	0.00269	0.0302	0.38	0	0	0	0	0	1	-360	360;
    49	66	0.018	0.0919	0.0248	0	0	0	0	0	1	-360	360;
    49	66	0.018	0.0919	0.0248	0	0	0	0	0	1	-360	360;
    62	66	0.0482	0.218	0.0578	0	0	0	0	0	1	-360	360;
    62	67	0.0258	0.117	0.031	0	0	0	0	0	1	-360	360;
    65	66	0	0.037	0	0	0	0	0.935	0	1	-360	360;
    66	67	0.0224	0.1015	0.02682	0	0	0	0	0	1	-360	360;
    65	68	0.00138	0.016	0.638	0	0	0	0	0	1	-360	360;
    47	69	0.0844	0.2778	0.07092	0	0	0	0	0	1	-360	360;
    49	69	0.0985	0.324	0.0828	0	0	0	0	0	1	-360	360;
    68	69	0	0.037	0	0	0	0	0.935	0	1	-360	360;
    69	70	0.03	0.127	0.122	0	0	0	0	0	1	-360	360;
    24	70	0.00221	0.4115	0.10198	0	0	0	0	0	1	-360	360;
    70	71	0.00882	0.0355	0.00878	0	0	0	0	0	1	-360	360;
    24	72	0.0488	0.196	0.0488	0	0	0	0	0	1	-360	360;
    71	72	0.0446	0.18	0.04444	0	0	0	0	0	1	-360	360;
    71	73	0.00866	0.0454	0.01178	0	0	0	0	0	1	-360	360;
    70	74	0.0401	0.1323	0.03368	0	0	0	0	0	1	-360	360;
    70	75	0.0428	0.141	0.036	0	0	0	0	0	1	-360	360;
    69	75	0.0405	0.122	0.124	0	0	0	0	0	1	-360	360;
    74	75	0.0123	0.0406	0.01034	0	0	0	0	0	1	-360	360;
    76	77	0.0444	0.148	0.0368	0	0	0	0	0	1	-360	360;
    69	77	0.0309	0.101	0.1038	0	0	0	0	0	1	-360	360;
    75	77	0.0601	0.1999	0.04978	0	0	0	0	0	1	-360	360;
    77	78	0.00376	0.0124	0.01264	0	0	0	0	0	1	-360	360;
    78	79	0.00546	0.0244	0.00648	0	0	0	0	0	1	-360	360;
    77	80	0.017	0.0485	0.0472	0	0	0	0	0	1	-360	360;
    77	80	0.0294	0.105	0.0228	0	0	0	0	0	1	-360	360;
    79	80	0.0156	0.0704	0.0187	0	0	0	0	0	1	-360	360;
    68	81	0.00175	0.0202	0.808	0	0	0	0	0	1	-360	360;
    81	80	0	0.037	0	0	0	0	0.935	0	1	-360	360;
    77	82	0.0298	0.0853	0.08174	0	0	0	0	0	1	-360	360;
    82	83	0.0112	0.03665	0.03796	0	0	0	0	0	1	-360	360;
    83	84	0.0625	0.132	0.0258	0	0	0	0	0	1	-360	360;
    83	85	0.043	0.148	0.0348	0	0	0	0	0	1	-360	360;
    84	85	0.0302	0.0641	0.01234	0	0	0	0	0	1	-360	360;
    85	86	0.035	0.123	0.0276	0	0	0	0	0	1	-360	360;
    86	87	0.02828	0.2074	0.0445	0	0	0	0	0	1	-360	360;
    85	88	0.02	0.102	0.0276	0	0	0	0	0	1	-360	360;
    85	89	0.0239	0.173	0.047	0	0	0	0	0	1	-360	360;
    88	89	0.0139	0.0712	0.01934	0	0	0	0	0	1	-360	360;
    89	90	0.0518	0.188	0.0528	0	0	0	0	0	1	-360	360;
    89	90	0.0238	0.0997	0.106	0	0	0	0	0	1	-360	360;
    90	91	0.0254	0.0836	0.0214	0	0	0	0	0	1	-360	360;
    89	92	0.0099	0.0505	0.0548	0	0	0	0	0	1	-360	360;
    89	92	0.0393	0.1581	0.0414	0	0	0	0	0	1	-360	360;
    91	92	0.0387	0.1272	0.03268	0	0	0	0	0	1	-360	360;
    92	93	0.0258	0.0848	0.0218	0	0	0	0	0	1	-360	360;
    92	94	0.0481	0.158	0.0406	0	0	0	0	0	1	-360	360;
    93	94	0.0223	0.0732	0.01876	0	0	0	0	0	1	-360	360;
    94	95	0.0132	0.0434	0.0111	0	0	0	0	0	1	-360	360;
    80	96	0.0356	0.182	0.0494	0	0	0	0	0	1	-360	360;
    82	96	0.0162	0.053	0.0544	0	0	0	0	0	1	-360	360;
    94	96	0.0269	0.0869	0.023	0	0	0	0	0	1	-360	360;
    80	97	0.0183	0.0934	0.0254	0	0	0	0	0	1	-360	360;
    80	98	0.0238	0.108	0.0286	0	0	0	0	0	1	-360	360;
    80	99	0.0454	0.206	0.0546	0	0	0	0	0	1	-360	360;
    92	100	0.0648	0.295	0.0472	0	0	0	0	0	1	-360	360;
    94	100	0.0178	0.058	0.0604	0	0	0	0	0	1	-360	360;
    95	96	0.0171	0.0547	0.01474	0	0	0	0	0	1	-360	360;
    96	97	0.0173	0.0885	0.024	0	0	0	0	0	1	-360	360;
    98	100	0.0397	0.179	0.0476	0	0	0	0	0	1	-360	360;
    99	100	0.018	0.0813	0.0216	0	0	0	0	0	1	-360	360;
    100	101	0.0277	0.1262	0.0328	0	0	0	0	0	1	-360	360;
    92	102	0.0123	0.0559	0.01464	0	0	0	0	0	1	-360	360;
    101	102	0.0246	0.112	0.0294	0	0	0	0	0	1	-360	360;
    100	103	0.016	0.0525	0.0536	0	0	0	0	0	1	-360	360;
    100	104	0.0451	0.204	0.0541	0	0	0	0	0	1	-360	360;
    103	104	0.0466	0.1584	0.0407	0	0	0	0	0	1	-360	360;
    103	105	0.0535	0.1625	0.0408	0	0	0	0	0	1	-360	360;
    100	106	0.0605	0.229	0.062	0	0	0	0	0	1	-360	360;
    104	105	0.00994	0.0378	0.00986	0	0	0	0	0	1	-360	360;
    105	106	0.014	0.0547	0.01434	0	0	0	0	0	1	-360	360;
    105	107	0.053	0.183	0.0472	0	0	0	0	0	1	-360	360;
    105	108	0.0261	0.0703	0.01844	0	0	0	0	0	1	-360	360;
    106	107	0.053	0.183	0.0472	0	0	0	0	0	1	-360	360;
    108	109	0.0105	0.0288	0.0076	0	0	0	0	0	1	-360	360;
    103	110	0.03906	0.1813	0.0461	0	0	0	0	0	1	-360	360;
    109	110	0.0278	0.0762	0.0202	0	0	0	0	0	1	-360	360;
    110	111	0.022	0.0755	0.02	0	0	0	0	0	1	-360	360;
    110	112	0.0247	0.064	0.062	0	0	0	0	0	1	-360	360;
    17	113	0.00913	0.0301	0.00768	0	0	0	0	0	1	-360	360;
    32	113	0.0615	0.203	0.0518	0	0	0	0	0	1	-360	360;
    32	114	0.0135	0.0612	0.01628	0	0	0	0	0	1	-360	360;
    27	115	0.0164	0.0741	0.01972	0	0	0	0	0	1	-360	360;
    114	115	0.0023	0.0104	0.00276	0	0	0	0	0	1	-360	360;
    68	116	0.00034	0.00405	0.164	0	0	0	0	0	1	-360	360;
    12	117	0.0329	0.14	0.0358	0	0	0	0	0	1	-360	360;
    75	118	0.0145	0.0481	0.01198	0	0	0	0	0	1	-360	360;
    76	118	0.0164	0.0544	0.01356	0	0	0	0	0	1	-360	360;
    ];
%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
    %     2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.0222222222	20	0;
    2	0	0	3	0.117647059	20	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.0454545455	20	0;
    %     2	0	0	3	0.0318471338	20	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	1.42857143	20	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.526315789	20	0;
    2	0	0	3	0.0490196078	20	0;
    2	0	0	3	0.208333333	20	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.064516129	20	0;
    2	0	0	3	0.0625	20	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.0255754476	20	0;
    2	0	0	3	0.0255102041	20	0;
    2	0	0	3	0.0193648335	20	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.0209643606	20	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	2.5	20	0;
    2	0	0	3	0.0164744646	20	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.0396825397	20	0;
    2	0	0	3	0.25	20	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.277777778	20	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    2	0	0	3	0.01	40	0;
    ];

Fmax = [500 % changed
    500 ; % changed
    500 ;
    175 ;
    175 ;
    175 ;
    500 ;
    500 ;
    500 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    500 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    500 ;
    500 ;
    500 ;
    175 ;
    175 ;
    500 ;
    500 ; % changed to 500 (initial 175)
    500 ;
    175 ;
    175 ;
    140 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    175 ;
    500 ;
    500 ;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    500;
    175;
    175;
    500;
    500;
    500;
    500;
    500;
    500;
    500;
    175;
    175;
    500;
    175;
    500;
    175;
    175;
    500;
    500;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    500;
    175;
    175;
    175;
    175;
    175;
    175;
    500;
    500;
    175;
    500;
    500;
    200;
    200;
    175;
    175;
    175;
    500;
    500;
    175;
    175;
    500;
    500;
    500;
    175;
    500;
    500;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    200;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    500;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    175;
    500;
    175;
    175;
    175;
    500;
    175;
    175;
    175;
    ];

mpc.baseMVA = 100; % in MVA
%% Input Data Section
% System MVA base
baseMVA = mpc.baseMVA;
% Number of Nodes
bus = mpc.bus;
branch = mpc.branch;
gen = mpc.gen;
gencost = mpc.gencost;

Nbus = size(bus,1); % number of buses
Nlines = size(branch,1); % number of lines
NGen = size(gen,1); % number of generators

%% for each branch, cause conducting the DC Power Flow here, the element
%% of B matrix is computed here.
%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

stat = branch(:, BR_STATUS);                    %% ones at in-service branches
b = stat ./ branch(:, BR_X);                    %% series susceptance
tap = ones(Nlines, 1);                              %% default tap ratio = 1
i = find(branch(:, TAP));                       %% indices of non-zero tap ratios
tap(i) = branch(i, TAP);                        %% assign non-zero tap ratios
b = b ./ tap;

%% build connection matrix Cft = Cf - Ct for line and from - to buses
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses
i = [(1:Nlines)'; (1:Nlines)'];                         %% double set of row indices
Cft = sparse(i, [f;t], [ones(Nlines, 1); -ones(Nlines, 1)], Nlines, Nbus);    %% connection matrix

%% build Bf such that Bf * Va is the vector of real branch powers injected
%% at each branch's "from" bus
Bf = sparse(i, [f; t], [b; -b], Nlines, Nbus);    % = spdiags(b, 0, nl, nl) * Cft;

%% build Bbus
Bbus = Cft' * Bf;

%% load wind energy forecast errors data
% Before import the wind power forecast errors, Please shift the forecast 
% errors data to make sure zero mean.
load('wind.mat')
% load three vector for each individual wind farm
% G_error_1;
% G_error_2;
% G_error_3;


% generator located matrix
Cg = zeros(Nbus,NGen);
for i = 1:NGen
    Cg(gen(i,1),i) = 1;
end

% Got the load demand at each bus
Pd = bus(:,3);

% Define a Gamma Matrix (power projection matrix),  nl x nb for the line flow calculation
Gamma = Bf * inv(Bbus);
% set a wind energy injection vector in dimension Nbus x 1
Pwinj = zeros(Nbus,1);

% The nominal wind power injection at bus 1 is 500 MW
Pwinj(1) = 500; % in MW
% The nominal wind power injection at bus 9 is 500 MW
Pwinj(9) = 500; % in MW
% The nominal wind power injection at bus 26 is 800 MW
Pwinj(26) = 800; % in MW

% set line flow limitation in MW
Fmax(7) = 600;
Fmax(37) = 500;  
Fmax(38) = 500;
Fmax(54) = 500; 
Fmax(96) = 300; 

% save the power project matrix for N-1 security constraints
GammaN1 = Gamma;

%% Get the number of loads
loadlocation = find(bus(:,3)>0);
load = mpc.bus(find(bus(:,3)>0),[1,3]);
Nload = size(load,1);

Ns = 30;% Number of scenarios for distributionally robust optimization
alpha = 0.01; % tolerant for line flow constraint CVaR
rho_Matrix = [1,10,30:30:900];

%% parameter check before start
disp(['program and system setup done !']);
disp(['Please check the following settings before start optimization']);
disp(['number of buses:    ', num2str(Nbus)]);
disp(['number of lines:   ', num2str(Nlines)]);
disp(['number of generators:    ', num2str(NGen)]);
disp(['number of scenarios for distributionally robust stochastic OPF:  ', num2str(Ns)]);
disp(['range of weight factors:', 'start from ', num2str(min(rho_Matrix)),', up to ', num2str(max(rho_Matrix))]);
disp(['wind injection locations: ', '  bus 1,  ','    bus 9,  ','    bus 26']);
disp(['wind nominal injection:   ', '  500 (MW), ','  500 (MW), ','  800 (MW)']);
disp(['..................................................................................'])

disp(['press any keys to continue .......']);
pause

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  data-based distributionally robust stochastic OPF   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Distribution Robust Optimization (wasserstein distance = 0.00)
Varepsilon = 0.000; % wasserstein distance
Num = 1;% location index of rho matrix
for  i = 1:length(rho_Matrix)
    rho = rho_Matrix(i);
    disp(['..................................................................................'])
    disp(['Wasserstein distance:    ', num2str(Varepsilon)]);
    disp(['Current weight factor:   ', num2str(rho)]);
    
    cvx_begin
    cvx_solver MOSEK
    
    variable D1(NGen) % Reserve policy to forecast errors of wind farm #1
    variable D2(NGen) % Reserve policy to forecast errors of wind farm #2
    variable D3(NGen) % Reserve policy to forecast errors of wind farm #3
    variable e(NGen) % scheduled power adjustment of generators
    variable d1(NGen, Nload) % reserve policy respond to the losing a load
    variable eN1(NGen, NGen) % scheduled power adjustment of generators (N-1 security)
    variable ea(NGen) % respond to line 9 trip, node 10 to node 9, lost generator 4
    variable eb(NGen) % respond to line 7 trip, node 9  to node 8, lost generator 4 and wind at bus 9
    variable ec(NGen) % respond to line 184 trip, node 12 to node 117, lost load at bus 117
    variable ed(NGen) % respond to line 176 trip, node 111 to node 110, lost generator 51
    variable ee(NGen) % respond to line 177 trip, node 112 to node 110, lost generator 52
    variable ef(NGen) % respond to line 134 trip, node 86 to node 87, lost generator 39
    variable eg(NGen) % respond to line 133 trip, node 86 to node 85, lost generator 39
    variable eh(NGen) % respond to line 113 trip, node 73 to node 71, lost generator 32
    
    % CVaR auxillary variables
    variable t11
    variable t21
    variable t31
    variable t41
    variable t51
     
    % distributionally robust auxillary variables
    variable lambda11
    variable lambda21
    variable lambda31
    variable lambda41
    variable lambda51
    
    % distributionally robust auxillary variables
    variable s11(Ns)
    variable s21(Ns)
    variable s31(Ns)
    variable s41(Ns)
    variable s51(Ns)

    F = 0;
    % opteration costs
    for i = 1:NGen
        F = F +  sum(gencost(i,5)*(D1(i).*G_error_1(:,1) + D2(i).*G_error_2(:,1) + D3(i).*G_error_3(:,1) + e(i)).^2 + (gencost(i,6))*(D1(i).*G_error_1(:,1) + D2(i).*G_error_2(:,1) + D3(i).*G_error_3(:,1) + e(i)) + gencost(i,7));
    end
    F = F./Ns;
    
    % distributoinally robust optimization (constraint violation)
    F2 = sum(s11(:)+s21(:)+s31(:)+s41(:)+s51(:))/Ns + (lambda11+lambda21+lambda31+lambda41+lambda51)*Varepsilon;% + sum(s1(:))/Ns + lambda1*Varepsilon;
    F = F + F2;
    
    minimize F
    subject to
    
    % calculate total load demand over system
    Pdemand = sum(mpc.bus(:,3));
    
    % sum-up reserve policies inequality constraints
    sum(D1(:)) + 1 == 0;
    sum(D2(:)) + 1 == 0;
    sum(D3(:)) + 1 == 0;
    Pwinj_error = zeros(Nbus,1);
    
    for i = 1:Ns
        % wind power injection with forecast errors
        Pwinj_error(9)  = Pwinj(9) + G_error_1(i,1);
        Pwinj_error(26) = Pwinj(26) + G_error_2(i,1);
        Pwinj_error(1)  = Pwinj(1) + G_error_3(i,1);
        
        rho*(-Gamma(7,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(7) + t11 - t11*alpha) <= s11(i);
        rho*(-t11*alpha) <= s11(i);
     
        % Distribution Robust Optimization on line 37
        rho*(Gamma(37,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(37) + t21 - t21*alpha) <= s21(i);
        rho*(-t21*alpha) <= s21(i);

        % Distribution Robust Optimization on line 38
        rho*(Gamma(38,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(38) + t31 - t31*alpha) <= s31(i);
        rho*(-t31*alpha) <= s31(i);
        
        % Distribution Robust Optimization on line 54
        rho*(Gamma(54,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(54) + t41 - t41*alpha) <= s41(i);
        rho*(-t41*alpha) <= s41(i);

        % Distribution Robust Optimization on line 96
        rho*(Gamma(96,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(96) + t51 - t51*alpha) <= s51(i);
        rho*(-t51*alpha) <= s51(i);
    end

    % Distribution Robust Optimization on Line 7
    norm([rho*Gamma(7,:)*Cg.*D1',rho*Gamma(7,:)*Cg.*D2',rho*Gamma(7,:)*Cg.*D3', rho*Gamma(7,:)*Cg.*e'],Inf) <= lambda11;
    norm(0,Inf) <= lambda11;
    
    % Distribution Robust Optimization on Line 37
    norm([-rho*Gamma(37,:)*Cg.*D1',-rho*Gamma(37,:)*Cg.*D2',-rho*Gamma(37,:)*Cg.*D3', -rho*Gamma(37,:)*Cg.*e'],Inf) <= lambda21;
    norm(0,Inf) <= lambda21;
    
    % Distribution Robust Optimization on Line 38
    norm([-rho*Gamma(38,:)*Cg.*D1',-rho*Gamma(38,:)*Cg.*D2',-rho*Gamma(38,:)*Cg.*D3', -rho*Gamma(38,:)*Cg.*e'],Inf) <= lambda31;
    norm(0,Inf) <= lambda31;
    
    % Distribution Robust Optimization on Line 54
    norm([-rho*Gamma(54,:)*Cg.*D1',-rho*Gamma(54,:)*Cg.*D2',-rho*Gamma(54,:)*Cg.*D3', -rho*Gamma(54,:)*Cg.*e'],Inf) <= lambda41;
    norm(0,Inf) <= lambda41;
    
    % Distribution Robust Optimization on Line 96
    norm([-rho*Gamma(96,:)*Cg.*D1',-rho*Gamma(96,:)*Cg.*D2',-rho*Gamma(96,:)*Cg.*D3', -rho*Gamma(96,:)*Cg.*e'],Inf) <= lambda51;
    norm(0,Inf) <= lambda51;
    
    % schedule power output of generators
    e>=0;
    %% Added the N-1 Security Constraints
    % N-1 security  constraints can also be formulated by distributionally
    % robust optimization. For simplcity, we formulate N-1 security
    % constraints here by sample-average approxiamtion
    
    % average all forecasting errors
    G_error_ave_1 = mean(G_error_1);
    G_error_ave_2 = mean(G_error_2);
    G_error_ave_3 = mean(G_error_3);
    
    
    % Case 1: if we lost one load
    for i = 1:Nload
        Pdismatch = load(i);
        
        Plimit = Fmax;
        GammaN1 = Gamma;

        PdN1 = Pd;
        PdN1(loadlocation(i)) = 0;

        GammaN1*(-PdN1 + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + d1(:,i).*Pdismatch) + Pwinj) <= Plimit;
        -GammaN1*(-PdN1 + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + d1(:,i).*Pdismatch) + Pwinj) <= Plimit;
        
        sum(d1(:,i)) == 1;
        sum(Pd) == sum(e(:)) + sum(Pwinj(:));
    end
    
    % case 2: if we lost one generator
    for i = 1:NGen
        
        CgN1 = zeros(Nbus,NGen);
        Plimit = Fmax;
        GammaN1 = Gamma;

        for j = 1:NGen
            CgN1(gen(j,1),j) = 1;
        end
        CgN1(gen(i,1),i) = 0;
        
        GammaN1*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eN1(:,i)) + Pwinj) <= Plimit;
        -GammaN1*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eN1(:,i)) + Pwinj) <= Plimit;
        eN1(i,i) == 0;
        sum(eN1(:,i)) == e(i);
    end
    
    %   case 3: if we lost one line
    for i = 1:Nlines
        branchN1 = branch;
        branchN1(i) = [];
        
        % re-formulate the Bf and Bbus matrix
        Nbus = size(bus,1); 
        NlinesN1 = size(branchN1,1);
        NGen = size(gen,1); 
        
        %% for each branch, cause conducting the DC Power Flow here, the element
        %% of B matrix is computed here.
        %% define named indices into bus, branch matrices
        [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
            VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
        [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
            TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
            ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
        
        stat = branchN1(:, BR_STATUS);                    %% ones at in-service branches
        b = stat ./ branchN1(:, BR_X);                    %% series susceptance
        tap = ones(NlinesN1, 1);                              %% default tap ratio = 1
        i = find(branchN1(:, TAP));                       %% indices of non-zero tap ratios
        tap(i) = branchN1(i, TAP);                        %% assign non-zero tap ratios
        b = b ./ tap;
        
        %% build connection matrix Cft = Cf - Ct for line and from - to buses
        f = branchN1(:, F_BUS);                           %% list of "from" buses
        t = branchN1(:, T_BUS);                           %% list of "to" buses
        i = [(1:NlinesN1)'; (1:NlinesN1)'];                         %% double set of row indices
        Cft = sparse(i, [f;t], [ones(NlinesN1, 1); -ones(NlinesN1, 1)], NlinesN1, Nbus);    %% connection matrix
        
        %% build Bf such that Bf * Va is the vector of real branch powers injected
        %% at each branch's "from" bus
        Bf = sparse(i, [f; t], [b; -b], NlinesN1, Nbus);    % = spdiags(b, 0, nl, nl) * Cft;
        
        %% build Bbus
        Bbus = Cft' * Bf;
        
        GammaN11 = Bf * inv(Bbus);
        Plimit = Fmax;
        Plimit(i) = [];
        
        if i == 9 % if trip line 9, lost generator 4
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(4,1),4) = 0;
            GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ea) + Pwinj) <= Plimit;
            -GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ea) + Pwinj) <= Plimit;
            ea(4) == 0;
            sum(ea) == e(4);
            
        elseif i == 7 % if trip line 7, lost generator 4, and wind injection at bus 9
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(4,1),4) = 0;
            PwinjN1 = Pwinj;
            PwinjN1(9) = 0;
            
            GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eb) + PwinjN1) <= Plimit;
            -GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eb) + PwinjN1) <= Plimit;
            eb(4) == 0;
            sum(eb) == e(4) + Pwinj(9);
            
        elseif i == 184 % if trip line 184, lost load at bus 117
            PdN1 = Pd;
            PdN1(117)= 0;
            GammaN11*(-PdN1 + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ec) + PwinjN1) <= Plimit;
            -GammaN11*(-PdN1+ Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ec) + PwinjN1) <= Plimit;
            sum(ec) == Pd(117);
            
        elseif i == 176 % if trip line 176, lost generator 51
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(51,1),51) = 0;
            GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ed) + Pwinj) <= Plimit;
            -GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ed) + Pwinj) <= Plimit;
            ed(51) == 0;
            sum(ed) == e(51);    
            
        elseif i == 177 % if trip line 177, lost generator 52, lost load at bus 112
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            PdN1 = Pd;
            PdN1(112)= 0;
            CgN1(gen(52,1),52) = 0;
            GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ee) + Pwinj) <= Plimit;
            -GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ee) + Pwinj) <= Plimit;
            
            ee(52) == 0;
            sum(ee) == e(52)+Pd(112); 
            
       elseif i == 134 % if trip line 134, lost generator 39
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(39,1),39) = 0;
            GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ef) + Pwinj) <= Plimit;
            -GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ef) + Pwinj) <= Plimit;
            ef(39) == 0;
            sum(ef) == e(39);
            
      elseif i == 133 % if trip line 133, lost generator 39, lost load at 86
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            PdN1 = Pd;
            PdN1(86)= 0;
            
            CgN1(gen(39,1),39) = 0;
            GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eg) + Pwinj) <= Plimit;
            -GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eg) + Pwinj) <= Plimit;
            eg(39) == 0;
            sum(eg) == e(39) + Pd(86);
            
      elseif i == 113 % if trip line 113, lost generator 32, lost load at 73
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(32,1),32) = 0;
            PdN1 = Pd;
            PdN1(73)= 0;
            
            GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eh) + Pwinj) <= Plimit;
            -GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eh) + Pwinj) <= Plimit;
            eh(39) == 0;
            sum(eh) == e(32) + Pd(73);
            
        else
            GammaN11*(-Pd + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e) + Pwinj) <= Plimit;
            -GammaN11*(-Pd + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e) + Pwinj) <= Plimit;
        end
    end
    
    cvx_end
    
    %% Data save and plot
    % save total operational cost and sum-up CVaR values of all line
    % constraints violation
    Cost_Data_0(Num) = cvx_optval - sum(s11(:)+s21(:)+s31(:)+s41(:)+s51(:))/Ns + (lambda11+lambda21+lambda31+lambda41+lambda51)*Varepsilon;    
    CVaR_Data_0(Num) = (sum(s11(:)+s21(:)+s31(:)+s41(:)+s51(:))/Ns + (lambda11+lambda21+lambda31+lambda41+lambda51)*Varepsilon)/rho;% + sum(s1(:))/Ns + lambda1*Varepsilon;
   
    % save CVaR value of each line constraint violation
    CVaR_Line1_0(Num) = (sum(s11(:))/(Ns) + lambda11*Varepsilon)/rho;
    CVaR_Line2_0(Num) = (sum(s21(:))/(Ns) + lambda21*Varepsilon)/rho;
    CVaR_Line3_0(Num) = (sum(s31(:))/(Ns) + lambda31*Varepsilon)/rho;
    CVaR_Line4_0(Num) = (sum(s41(:))/(Ns) + lambda41*Varepsilon)/rho;
    CVaR_Line5_0(Num) = (sum(s51(:))/(Ns) + lambda51*Varepsilon)/rho;
    
    % save reserve policies
    D1_Data_0(Num,:) = D1;
    D2_Data_0(Num,:) = D2;
    D3_Data_0(Num,:) = D3;
    
    % save schedule power output
    e_Data_0(Num,:) = e;
    
    % save CVaR auxillary variables
    t11_Data_0(Num) = t11;
    t21_Data_0(Num) = t21;
    t31_Data_0(Num) = t31;
    t41_Data_0(Num) = t41;   
    t51_Data_0(Num) = t51;
    
    % save optimal status for each weight factor \rho
    OptimalStatus_0{Num} = cvx_status;
    
    % average actual line flow (mean of forecast errors is zero)
    LineFlow_0(Num,:) = Gamma*(- Pd + Cg*e+ Pwinj);
    
%     cvx_status
%     cvx_optval
%     Cost_Data_0
%     CVaR_Data_0
    Num = Num + 1; % move to next weight factor \rho
    toc
end


%% Distribution Robust Optimization (wasserstein distance = 0.02)
Varepsilon = 0.02; % wasserstein distance
Num = 1;% location index of rho matrix
for  i = 1:length(rho_Matrix)
    rho = rho_Matrix(i);
    disp(['..................................................................................'])
    disp(['Wasserstein distance:    ', num2str(Varepsilon)]);
    disp(['Current weight factor:   ', num2str(rho)]);
    
    cvx_begin
    cvx_solver MOSEK
    
    variable D1(NGen) % Reserve policy to forecast errors of wind injection 1
    variable D2(NGen) % Reserve policy to forecast errors of wind injection 2
    variable D3(NGen) % Reserve policy to forecast errors of wind injection 3
    variable e(NGen) % power schedule of generators
    variable d1(NGen, Nload) % reserve policy respond to the losing a load
    variable eN1(NGen, NGen) % power schedule of generators (N-1 security)
    variable ea(NGen) % respond to the line 9 trip, node 10 to node 9, loss generator 4
    variable eb(NGen) % respond to the line 7 trip, node 9  to node 8, loss generator 4 and wind at bus 9
    variable ec(NGen) % respond to the line 184 trip, node 12 to node 117, loss load at bus 117
    variable ed(NGen) % respond to the line 176 trip, node 111 to node 110, loss generator 51
    variable ee(NGen) % respond to the line 177 trip, node 112 to node 110, loss generator 52
    variable ef(NGen) % respond to the line 134 trip, node 86 to node 87, loss generator 39
    variable eg(NGen) % respond to the line 133 trip, node 86 to node 85, loss generator 39
    variable eh(NGen) % respond to the line 113 trip, node 73 to node 71, loss generator 32
    
    % CVaR auxillary variables
    variable t11
    variable t21
    variable t31
    variable t41
    variable t51
     
    % distributionally robust auxillary variables
    variable lambda11
    variable lambda21
    variable lambda31
    variable lambda41
    variable lambda51
    
    % distributionally robust auxillary variables
    variable s11(Ns)
    variable s21(Ns)
    variable s31(Ns)
    variable s41(Ns)
    variable s51(Ns)

    F = 0;
    % opteration costs
    for i = 1:NGen
        F = F +  sum(gencost(i,5)*(D1(i).*G_error_1(:,1) + D2(i).*G_error_2(:,1) + D3(i).*G_error_3(:,1) + e(i)).^2 + (gencost(i,6))*(D1(i).*G_error_1(:,1) + D2(i).*G_error_2(:,1) + D3(i).*G_error_3(:,1) + e(i)) + gencost(i,7));
    end
    F = F./Ns;
    
    % distributoinally robust optimization (constraint violation)
    F2 = sum(s11(:)+s21(:)+s31(:)+s41(:)+s51(:))/Ns + (lambda11+lambda21+lambda31+lambda41+lambda51)*Varepsilon;% + sum(s1(:))/Ns + lambda1*Varepsilon;
    F = F + F2;
    
    minimize F
    subject to
    
    % calculate total load demand over system
    Pdemand = sum(mpc.bus(:,3));
    
    % sum-up reserve policies inequality constraints
    sum(D1(:)) + 1 == 0;
    sum(D2(:)) + 1 == 0;
    sum(D3(:)) + 1 == 0;
    Pwinj_error = zeros(Nbus,1);
    
    for i = 1:Ns
        % wind power injection with forecast errors
        Pwinj_error(9)  = Pwinj(9) + G_error_1(i,1);
        Pwinj_error(26) = Pwinj(26) + G_error_2(i,1);
        Pwinj_error(1)  = Pwinj(1) + G_error_3(i,1);
        
        rho*(-Gamma(7,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(7) + t11 - t11*alpha) <= s11(i);
        rho*(-t11*alpha) <= s11(i);
     
        % Distribution Robust Optimization on line 37
        rho*(Gamma(37,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(37) + t21 - t21*alpha) <= s21(i);
        rho*(-t21*alpha) <= s21(i);

        % Distribution Robust Optimization on line 38
        rho*(Gamma(38,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(38) + t31 - t31*alpha) <= s31(i);
        rho*(-t31*alpha) <= s31(i);
        
        % Distribution Robust Optimization on line 54
        rho*(Gamma(54,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(54) + t41 - t41*alpha) <= s41(i);
        rho*(-t41*alpha) <= s41(i);

        % Distribution Robust Optimization on line 96
        rho*(Gamma(96,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(96) + t51 - t51*alpha) <= s51(i);
        rho*(-t51*alpha) <= s51(i);
    end

    % Distribution Robust Optimization on Line 7
    norm([rho*Gamma(7,:)*Cg.*D1',rho*Gamma(7,:)*Cg.*D2',rho*Gamma(7,:)*Cg.*D3', rho*Gamma(7,:)*Cg.*e'],Inf) <= lambda11;
    norm(0,Inf) <= lambda11;
    
    % Distribution Robust Optimization on Line 37
    norm([-rho*Gamma(37,:)*Cg.*D1',-rho*Gamma(37,:)*Cg.*D2',-rho*Gamma(37,:)*Cg.*D3', -rho*Gamma(37,:)*Cg.*e'],Inf) <= lambda21;
    norm(0,Inf) <= lambda21;
    
    % Distribution Robust Optimization on Line 38
    norm([-rho*Gamma(38,:)*Cg.*D1',-rho*Gamma(38,:)*Cg.*D2',-rho*Gamma(38,:)*Cg.*D3', -rho*Gamma(38,:)*Cg.*e'],Inf) <= lambda31;
    norm(0,Inf) <= lambda31;
    
    % Distribution Robust Optimization on Line 54
    norm([-rho*Gamma(54,:)*Cg.*D1',-rho*Gamma(54,:)*Cg.*D2',-rho*Gamma(54,:)*Cg.*D3', -rho*Gamma(54,:)*Cg.*e'],Inf) <= lambda41;
    norm(0,Inf) <= lambda41;
    
    % Distribution Robust Optimization on Line 96
    norm([-rho*Gamma(96,:)*Cg.*D1',-rho*Gamma(96,:)*Cg.*D2',-rho*Gamma(96,:)*Cg.*D3', -rho*Gamma(96,:)*Cg.*e'],Inf) <= lambda51;
    norm(0,Inf) <= lambda51;
    
    % schedule power output of generators
    e>=0;
    %% Added the N-1 Security Constraints
    % N-1 security  constraints can also be formulated by distributionally
    % robust optimization. For simplcity, we formulate N-1 security
    % constraints here by sample-average approxiamtion
    
    % average all forecasting errors
    G_error_ave_1 = mean(G_error_1);
    G_error_ave_2 = mean(G_error_2);
    G_error_ave_3 = mean(G_error_3);
    
    
    % Case 1: if we lost one load
    for i = 1:Nload
        Pdismatch = load(i);
        
        Plimit = Fmax;
        GammaN1 = Gamma;

        PdN1 = Pd;
        PdN1(loadlocation(i)) = 0;

        GammaN1*(-PdN1 + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + d1(:,i).*Pdismatch) + Pwinj) <= Plimit;
        -GammaN1*(-PdN1 + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + d1(:,i).*Pdismatch) + Pwinj) <= Plimit;
        
        sum(d1(:,i)) == 1;
        sum(Pd) == sum(e(:)) + sum(Pwinj(:));
    end
    
    % case 2: if we lost one generator
    for i = 1:NGen
        
        CgN1 = zeros(Nbus,NGen);
        Plimit = Fmax;
        GammaN1 = Gamma;

        for j = 1:NGen
            CgN1(gen(j,1),j) = 1;
        end
        CgN1(gen(i,1),i) = 0;
        
        GammaN1*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eN1(:,i)) + Pwinj) <= Plimit;
        -GammaN1*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eN1(:,i)) + Pwinj) <= Plimit;
        eN1(i,i) == 0;
        sum(eN1(:,i)) == e(i);
    end
    
    %   case 3: if we lost one line
    for i = 1:Nlines
        branchN1 = branch;
        branchN1(i) = [];
        
        % re-formulate the Bf and Bbus matrix
        Nbus = size(bus,1); 
        NlinesN1 = size(branchN1,1);
        NGen = size(gen,1); 
        
        %% for each branch, cause conducting the DC Power Flow here, the element
        %% of B matrix is computed here.
        %% define named indices into bus, branch matrices
        [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
            VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
        [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
            TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
            ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
        
        stat = branchN1(:, BR_STATUS);                    %% ones at in-service branches
        b = stat ./ branchN1(:, BR_X);                    %% series susceptance
        tap = ones(NlinesN1, 1);                              %% default tap ratio = 1
        i = find(branchN1(:, TAP));                       %% indices of non-zero tap ratios
        tap(i) = branchN1(i, TAP);                        %% assign non-zero tap ratios
        b = b ./ tap;
        
        %% build connection matrix Cft = Cf - Ct for line and from - to buses
        f = branchN1(:, F_BUS);                           %% list of "from" buses
        t = branchN1(:, T_BUS);                           %% list of "to" buses
        i = [(1:NlinesN1)'; (1:NlinesN1)'];                         %% double set of row indices
        Cft = sparse(i, [f;t], [ones(NlinesN1, 1); -ones(NlinesN1, 1)], NlinesN1, Nbus);    %% connection matrix
        
        %% build Bf such that Bf * Va is the vector of real branch powers injected
        %% at each branch's "from" bus
        Bf = sparse(i, [f; t], [b; -b], NlinesN1, Nbus);    % = spdiags(b, 0, nl, nl) * Cft;
        
        %% build Bbus
        Bbus = Cft' * Bf;
        
        GammaN11 = Bf * inv(Bbus);
        Plimit = Fmax;
        Plimit(i) = [];
        
        if i == 9 % if trip line 9, lost generator 4
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(4,1),4) = 0;
            GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ea) + Pwinj) <= Plimit;
            -GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ea) + Pwinj) <= Plimit;
            ea(4) == 0;
            sum(ea) == e(4);
            
        elseif i == 7 % if trip line 7, lost generator 4, and wind injection at bus 9
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(4,1),4) = 0;
            PwinjN1 = Pwinj;
            PwinjN1(9) = 0;
            
            GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eb) + PwinjN1) <= Plimit;
            -GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eb) + PwinjN1) <= Plimit;
            eb(4) == 0;
            sum(eb) == e(4) + Pwinj(9);
            
        elseif i == 184 % if trip line 184, lost load at bus 117
            PdN1 = Pd;
            PdN1(117)= 0;
            GammaN11*(-PdN1 + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ec) + PwinjN1) <= Plimit;
            -GammaN11*(-PdN1+ Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ec) + PwinjN1) <= Plimit;
            sum(ec) == Pd(117);
            
        elseif i == 176 % if trip line 176, lost generator 51
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(51,1),51) = 0;
            GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ed) + Pwinj) <= Plimit;
            -GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ed) + Pwinj) <= Plimit;
            ed(51) == 0;
            sum(ed) == e(51);    
            
        elseif i == 177 % if trip line 177, lost generator 52, lost load at bus 112
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            PdN1 = Pd;
            PdN1(112)= 0;
            CgN1(gen(52,1),52) = 0;
            GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ee) + Pwinj) <= Plimit;
            -GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ee) + Pwinj) <= Plimit;
            
            ee(52) == 0;
            sum(ee) == e(52)+Pd(112); 
            
       elseif i == 134 % if trip line 134, lost generator 39
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(39,1),39) = 0;
            GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ef) + Pwinj) <= Plimit;
            -GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ef) + Pwinj) <= Plimit;
            ef(39) == 0;
            sum(ef) == e(39);
            
      elseif i == 133 % if trip line 133, lost generator 39, lost load at 86
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            PdN1 = Pd;
            PdN1(86)= 0;
            
            CgN1(gen(39,1),39) = 0;
            GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eg) + Pwinj) <= Plimit;
            -GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eg) + Pwinj) <= Plimit;
            eg(39) == 0;
            sum(eg) == e(39) + Pd(86);
            
      elseif i == 113 % if trip line 113, lost generator 32, lost load at 73
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(32,1),32) = 0;
            PdN1 = Pd;
            PdN1(73)= 0;
            
            GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eh) + Pwinj) <= Plimit;
            -GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eh) + Pwinj) <= Plimit;
            eh(39) == 0;
            sum(eh) == e(32) + Pd(73);
            
        else
            GammaN11*(-Pd + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e) + Pwinj) <= Plimit;
            -GammaN11*(-Pd + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e) + Pwinj) <= Plimit;
        end
    end
    
    cvx_end
    
    %% Data save and plot
    % save total operational cost and sum-up CVaR values of all line
    % constraints violation
    Cost_Data_1(Num) = cvx_optval - sum(s11(:)+s21(:)+s31(:)+s41(:)+s51(:))/Ns + (lambda11+lambda21+lambda31+lambda41+lambda51)*Varepsilon;    
    CVaR_Data_1(Num) = (sum(s11(:)+s21(:)+s31(:)+s41(:)+s51(:))/Ns + (lambda11+lambda21+lambda31+lambda41+lambda51)*Varepsilon)/rho;% + sum(s1(:))/Ns + lambda1*Varepsilon;
   
    % save CVaR value of each line constraint violation
    CVaR_Line1_1(Num) = (sum(s11(:))/(Ns) + lambda11*Varepsilon)/rho;
    CVaR_Line2_1(Num) = (sum(s21(:))/(Ns) + lambda21*Varepsilon)/rho;
    CVaR_Line3_1(Num) = (sum(s31(:))/(Ns) + lambda31*Varepsilon)/rho;
    CVaR_Line4_1(Num) = (sum(s41(:))/(Ns) + lambda41*Varepsilon)/rho;
    CVaR_Line5_1(Num) = (sum(s51(:))/(Ns) + lambda51*Varepsilon)/rho;
    
    % save reserve policies
    D1_Data_1(Num,:) = D1;
    D2_Data_1(Num,:) = D2;
    D3_Data_1(Num,:) = D3;
    
    % save schedule power output
    e_Data_1(Num,:) = e;
    
    % save CVaR auxillary variables
    t11_Data_1(Num) = t11;
    t21_Data_1(Num) = t21;
    t31_Data_1(Num) = t31;
    t41_Data_1(Num) = t41;   
    t51_Data_1(Num) = t51;
    
    % save optimal status for each weight factor \rho
    OptimalStatus_1{Num} = cvx_status;
    
    % average actual line flow (mean of forecast errors is zero)
    LineFlow_1(Num,:) = Gamma*(- Pd + Cg*e+ Pwinj);
    
%     cvx_status
%     cvx_optval
%     Cost_Data_0
%     CVaR_Data_0
    Num = Num + 1; % move to next weight factor \rho
    toc
end

%% Distribution Robust Optimization (wasserstein distance = 0.04)
Varepsilon = 0.04; % wasserstein distance
Num = 1;% location index of rho matrix
for  i = 1:length(rho_Matrix)
    rho = rho_Matrix(i);
    disp(['..................................................................................'])
    disp(['Wasserstein distance:    ', num2str(Varepsilon)]);
    disp(['Current weight factor:   ', num2str(rho)]);
    
    cvx_begin
    cvx_solver MOSEK
    
    variable D1(NGen) % Reserve policy to forecast errors of wind injection 1
    variable D2(NGen) % Reserve policy to forecast errors of wind injection 2
    variable D3(NGen) % Reserve policy to forecast errors of wind injection 3
    variable e(NGen) % power schedule of generators
    variable d1(NGen, Nload) % reserve policy respond to the losing a load
    variable eN1(NGen, NGen) % power schedule of generators (N-1 security)
    variable ea(NGen) % respond to the line 9 trip, node 10 to node 9, loss generator 4
    variable eb(NGen) % respond to the line 7 trip, node 9  to node 8, loss generator 4 and wind at bus 9
    variable ec(NGen) % respond to the line 184 trip, node 12 to node 117, loss load at bus 117
    variable ed(NGen) % respond to the line 176 trip, node 111 to node 110, loss generator 51
    variable ee(NGen) % respond to the line 177 trip, node 112 to node 110, loss generator 52
    variable ef(NGen) % respond to the line 134 trip, node 86 to node 87, loss generator 39
    variable eg(NGen) % respond to the line 133 trip, node 86 to node 85, loss generator 39
    variable eh(NGen) % respond to the line 113 trip, node 73 to node 71, loss generator 32
    
    % CVaR auxillary variables
    variable t11
    variable t21
    variable t31
    variable t41
    variable t51
     
    % distributionally robust auxillary variables
    variable lambda11
    variable lambda21
    variable lambda31
    variable lambda41
    variable lambda51
    
    % distributionally robust auxillary variables
    variable s11(Ns)
    variable s21(Ns)
    variable s31(Ns)
    variable s41(Ns)
    variable s51(Ns)

    F = 0;
    % opteration costs
    for i = 1:NGen
        F = F +  sum(gencost(i,5)*(D1(i).*G_error_1(:,1) + D2(i).*G_error_2(:,1) + D3(i).*G_error_3(:,1) + e(i)).^2 + (gencost(i,6))*(D1(i).*G_error_1(:,1) + D2(i).*G_error_2(:,1) + D3(i).*G_error_3(:,1) + e(i)) + gencost(i,7));
    end
    F = F./Ns;
    
    % distributoinally robust optimization (constraint violation)
    F2 = sum(s11(:)+s21(:)+s31(:)+s41(:)+s51(:))/Ns + (lambda11+lambda21+lambda31+lambda41+lambda51)*Varepsilon;
    F = F + F2;
    
    minimize F
    subject to
    
    % calculate total load demand over system
    Pdemand = sum(mpc.bus(:,3));
    
    % sum-up reserve policies inequality constraints
    sum(D1(:)) + 1 == 0;
    sum(D2(:)) + 1 == 0;
    sum(D3(:)) + 1 == 0;
    Pwinj_error = zeros(Nbus,1);
    
    for i = 1:Ns
        % wind power injection with forecast errors
        Pwinj_error(9)  = Pwinj(9) + G_error_1(i,1);
        Pwinj_error(26) = Pwinj(26) + G_error_2(i,1);
        Pwinj_error(1)  = Pwinj(1) + G_error_3(i,1);
        
        rho*(-Gamma(7,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(7) + t11 - t11*alpha) <= s11(i);
        rho*(-t11*alpha) <= s11(i);
     
        % Distribution Robust Optimization on line 37
        rho*(Gamma(37,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(37) + t21 - t21*alpha) <= s21(i);
        rho*(-t21*alpha) <= s21(i);

        % Distribution Robust Optimization on line 38
        rho*(Gamma(38,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(38) + t31 - t31*alpha) <= s31(i);
        rho*(-t31*alpha) <= s31(i);
        
        % Distribution Robust Optimization on line 54
        rho*(Gamma(54,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(54) + t41 - t41*alpha) <= s41(i);
        rho*(-t41*alpha) <= s41(i);

        % Distribution Robust Optimization on line 96
        rho*(Gamma(96,:)*( - Pd + Cg*(D1.*G_error_1(i,1) + D2.*G_error_2(i,1) + D3.*G_error_3(i,1) + e) + Pwinj_error) - Fmax(96) + t51 - t51*alpha) <= s51(i);
        rho*(-t51*alpha) <= s51(i);
    end

    % Distribution Robust Optimization on Line 7
    norm([rho*Gamma(7,:)*Cg.*D1',rho*Gamma(7,:)*Cg.*D2',rho*Gamma(7,:)*Cg.*D3', rho*Gamma(7,:)*Cg.*e'],Inf) <= lambda11;
    norm(0,Inf) <= lambda11;
    
    % Distribution Robust Optimization on Line 37
    norm([-rho*Gamma(37,:)*Cg.*D1',-rho*Gamma(37,:)*Cg.*D2',-rho*Gamma(37,:)*Cg.*D3', -rho*Gamma(37,:)*Cg.*e'],Inf) <= lambda21;
    norm(0,Inf) <= lambda21;
    
    % Distribution Robust Optimization on Line 38
    norm([-rho*Gamma(38,:)*Cg.*D1',-rho*Gamma(38,:)*Cg.*D2',-rho*Gamma(38,:)*Cg.*D3', -rho*Gamma(38,:)*Cg.*e'],Inf) <= lambda31;
    norm(0,Inf) <= lambda31;
    
    % Distribution Robust Optimization on Line 54
    norm([-rho*Gamma(54,:)*Cg.*D1',-rho*Gamma(54,:)*Cg.*D2',-rho*Gamma(54,:)*Cg.*D3', -rho*Gamma(54,:)*Cg.*e'],Inf) <= lambda41;
    norm(0,Inf) <= lambda41;
    
    % Distribution Robust Optimization on Line 96
    norm([-rho*Gamma(96,:)*Cg.*D1',-rho*Gamma(96,:)*Cg.*D2',-rho*Gamma(96,:)*Cg.*D3', -rho*Gamma(96,:)*Cg.*e'],Inf) <= lambda51;
    norm(0,Inf) <= lambda51;
    
    % schedule power output of generators
    e>=0;
    %% Added the N-1 Security Constraints
    % N-1 security  constraints can also be formulated by distributionally
    % robust optimization. For simplcity, we formulate N-1 security
    % constraints here by sample-average approxiamtion
    
    % average all forecasting errors
    G_error_ave_1 = mean(G_error_1);
    G_error_ave_2 = mean(G_error_2);
    G_error_ave_3 = mean(G_error_3);
    
    % Case 1: if we lost one load
    for i = 1:Nload
        Pdismatch = load(i);
        
        Plimit = Fmax;
        GammaN1 = Gamma;

        PdN1 = Pd;
        PdN1(loadlocation(i)) = 0;

        GammaN1*(-PdN1 + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + d1(:,i).*Pdismatch) + Pwinj) <= Plimit;
        -GammaN1*(-PdN1 + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + d1(:,i).*Pdismatch) + Pwinj) <= Plimit;
        
        sum(d1(:,i)) == 1;
        sum(Pd) == sum(e(:)) + sum(Pwinj(:));
    end
    
    % case 2: if we lost one generator
    for i = 1:NGen
        
        CgN1 = zeros(Nbus,NGen);
        Plimit = Fmax;
        GammaN1 = Gamma;

        for j = 1:NGen
            CgN1(gen(j,1),j) = 1;
        end
        CgN1(gen(i,1),i) = 0;
        
        GammaN1*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eN1(:,i)) + Pwinj) <= Plimit;
        -GammaN1*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eN1(:,i)) + Pwinj) <= Plimit;
        eN1(i,i) == 0;
        sum(eN1(:,i)) == e(i);
    end
    
    %   case 3: if we lost one line
    for i = 1:Nlines
        branchN1 = branch;
        branchN1(i) = [];
        
        % re-formulate the Bf and Bbus matrix
        Nbus = size(bus,1); 
        NlinesN1 = size(branchN1,1);
        NGen = size(gen,1); 
        
        %% for each branch, cause conducting the DC Power Flow here, the element
        %% of B matrix is computed here.
        %% define named indices into bus, branch matrices
        [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
            VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
        [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
            TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
            ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
        
        stat = branchN1(:, BR_STATUS);                    %% ones at in-service branches
        b = stat ./ branchN1(:, BR_X);                    %% series susceptance
        tap = ones(NlinesN1, 1);                              %% default tap ratio = 1
        i = find(branchN1(:, TAP));                       %% indices of non-zero tap ratios
        tap(i) = branchN1(i, TAP);                        %% assign non-zero tap ratios
        b = b ./ tap;
        
        %% build connection matrix Cft = Cf - Ct for line and from - to buses
        f = branchN1(:, F_BUS);                           %% list of "from" buses
        t = branchN1(:, T_BUS);                           %% list of "to" buses
        i = [(1:NlinesN1)'; (1:NlinesN1)'];                         %% double set of row indices
        Cft = sparse(i, [f;t], [ones(NlinesN1, 1); -ones(NlinesN1, 1)], NlinesN1, Nbus);    %% connection matrix
        
        %% build Bf such that Bf * Va is the vector of real branch powers injected
        %% at each branch's "from" bus
        Bf = sparse(i, [f; t], [b; -b], NlinesN1, Nbus);    % = spdiags(b, 0, nl, nl) * Cft;
        
        %% build Bbus
        Bbus = Cft' * Bf;
        
        GammaN11 = Bf * inv(Bbus);
        Plimit = Fmax;
        Plimit(i) = [];
        
        if i == 9 % if trip line 9, lost generator 4
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(4,1),4) = 0;
            GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ea) + Pwinj) <= Plimit;
            -GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ea) + Pwinj) <= Plimit;
            ea(4) == 0;
            sum(ea) == e(4);
            
        elseif i == 7 % if trip line 7, lost generator 4, and wind injection at bus 9
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(4,1),4) = 0;
            PwinjN1 = Pwinj;
            PwinjN1(9) = 0;
            
            GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eb) + PwinjN1) <= Plimit;
            -GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eb) + PwinjN1) <= Plimit;
            eb(4) == 0;
            sum(eb) == e(4) + Pwinj(9);
            
        elseif i == 184 % if trip line 184, lost load at bus 117
            PdN1 = Pd;
            PdN1(117)= 0;
            GammaN11*(-PdN1 + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ec) + PwinjN1) <= Plimit;
            -GammaN11*(-PdN1+ Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ec) + PwinjN1) <= Plimit;
            sum(ec) == Pd(117);
            
        elseif i == 176 % if trip line 176, lost generator 51
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(51,1),51) = 0;
            GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ed) + Pwinj) <= Plimit;
            -GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ed) + Pwinj) <= Plimit;
            ed(51) == 0;
            sum(ed) == e(51);    
            
        elseif i == 177 % if trip line 177, lost generator 52, lost load at bus 112
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            PdN1 = Pd;
            PdN1(112)= 0;
            CgN1(gen(52,1),52) = 0;
            GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ee) + Pwinj) <= Plimit;
            -GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ee) + Pwinj) <= Plimit;
            
            ee(52) == 0;
            sum(ee) == e(52)+Pd(112); 
            
       elseif i == 134 % if trip line 134, lost generator 39
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(39,1),39) = 0;
            GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ef) + Pwinj) <= Plimit;
            -GammaN11*(-Pd + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + ef) + Pwinj) <= Plimit;
            ef(39) == 0;
            sum(ef) == e(39);
            
      elseif i == 133 % if trip line 133, lost generator 39, lost load at 86
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            PdN1 = Pd;
            PdN1(86)= 0;
            
            CgN1(gen(39,1),39) = 0;
            GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eg) + Pwinj) <= Plimit;
            -GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eg) + Pwinj) <= Plimit;
            eg(39) == 0;
            sum(eg) == e(39) + Pd(86);
            
      elseif i == 113 % if trip line 113, lost generator 32, lost load at 73
            CgN1 = zeros(Nbus,NGen);
            for j = 1:NGen
                CgN1(gen(j,1),j) = 1;
            end
            CgN1(gen(32,1),32) = 0;
            PdN1 = Pd;
            PdN1(73)= 0;
            
            GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eh) + Pwinj) <= Plimit;
            -GammaN11*(-PdN1 + CgN1*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e + eh) + Pwinj) <= Plimit;
            eh(39) == 0;
            sum(eh) == e(32) + Pd(73);
            
        else
            GammaN11*(-Pd + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e) + Pwinj) <= Plimit;
            -GammaN11*(-Pd + Cg*(D1.*G_error_ave_1  + D2.*G_error_ave_2 + D3.*G_error_ave_3 + e) + Pwinj) <= Plimit;
        end
    end
    
    cvx_end
    
    %% Data save and plot
    % save total operational cost and sum-up CVaR values of all line
    % constraints violation
    Cost_Data_2(Num) = cvx_optval - sum(s11(:)+s21(:)+s31(:)+s41(:)+s51(:))/Ns + (lambda11+lambda21+lambda31+lambda41+lambda51)*Varepsilon;    
    CVaR_Data_2(Num) = (sum(s11(:)+s21(:)+s31(:)+s41(:)+s51(:))/Ns + (lambda11+lambda21+lambda31+lambda41+lambda51)*Varepsilon)/rho;% + sum(s1(:))/Ns + lambda1*Varepsilon;
   
    % save CVaR value of each line constraint violation
    CVaR_Line1_2(Num) = (sum(s11(:))/(Ns) + lambda11*Varepsilon)/rho;
    CVaR_Line2_2(Num) = (sum(s21(:))/(Ns) + lambda21*Varepsilon)/rho;
    CVaR_Line3_2(Num) = (sum(s31(:))/(Ns) + lambda31*Varepsilon)/rho;
    CVaR_Line4_2(Num) = (sum(s41(:))/(Ns) + lambda41*Varepsilon)/rho;
    CVaR_Line5_2(Num) = (sum(s51(:))/(Ns) + lambda51*Varepsilon)/rho;
    
    % save reserve policies
    D1_Data_2(Num,:) = D1;
    D2_Data_2(Num,:) = D2;
    D3_Data_2(Num,:) = D3;
    
    % save schedule power output
    e_Data_2(Num,:) = e;
    
    % save CVaR auxillary variables
    t11_Data_2(Num) = t11;
    t21_Data_2(Num) = t21;
    t31_Data_2(Num) = t31;
    t41_Data_2(Num) = t41;   
    t51_Data_2(Num) = t51;
    
    % save optimal status for each weight factor \rho
    OptimalStatus_2{Num} = cvx_status;
    
    % average actual line flow (mean of forecast errors is zero)
    LineFlow_2(Num,:) = Gamma*(- Pd + Cg*e+ Pwinj);
    
%     cvx_status
%     cvx_optval
%     Cost_Data_0
%     CVaR_Data_0
    Num = Num + 1; % move to next weight factor \rho
    toc
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  results plot section   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1-st plot
% This plot visualizes the predicted tradeoff between operational cost and
% CVaR of line constraint violation

figure;
% Line 7
subplot(3,2,1)
plot(Cost_Data_0, CVaR_Line1_0,'o-','LineWidth',1.5)
hold on
plot(Cost_Data_1, CVaR_Line1_1,'*-','LineWidth',1.5)
hold on
plot(Cost_Data_2, CVaR_Line1_2,'s-','LineWidth',1.5)
hold on
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('Operation Cost in Dollars')
ylabel('CVaR')
title('Line 7')
set(gca,'FontSize',13);
grid

% Line 37
subplot(3,2,2)
plot(Cost_Data_0, CVaR_Line2_0,'o-','LineWidth',1.5)
hold on
plot(Cost_Data_1, CVaR_Line2_1,'*-','LineWidth',1.5)
hold on
plot(Cost_Data_2, CVaR_Line2_2,'s-','LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('Operation Cost in Dollars')
ylabel('CVaR')
title('Line 37')
set(gca,'FontSize',13);
grid

% Line 38
subplot(3,2,3)
plot(Cost_Data_0, CVaR_Line3_0,'o-','LineWidth',1.5)
hold on
plot(Cost_Data_1, CVaR_Line3_1,'*-','LineWidth',1.5)
hold on
plot(Cost_Data_2, CVaR_Line3_2,'s-','LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('Operation Cost in Dollars')
ylabel('CVaR')
title('Line 38')
set(gca,'FontSize',13);
grid

% Line 54
subplot(3,2,4)
plot(Cost_Data_0, CVaR_Line4_0,'o-','LineWidth',1.5)
hold on
plot(Cost_Data_1, CVaR_Line4_1,'*-','LineWidth',1.5)
hold on
plot(Cost_Data_2, CVaR_Line4_2,'s-','LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('Operation Cost in Dollars')
ylabel('CVaR')
title('Line 54')
set(gca,'FontSize',13);
grid

% Line 96
subplot(3,2,5)
plot(Cost_Data_0, CVaR_Line5_0,'o-','LineWidth',1.5)
hold on
plot(Cost_Data_1, CVaR_Line5_1,'*-','LineWidth',1.5)
hold on
plot(Cost_Data_2, CVaR_Line5_2,'s-','LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('Operation Cost in Dollars')
ylabel('CVaR')
title('Line 96')
set(gca,'FontSize',13);
grid

%% 2-nd plot
% This plot visualizes out-of-sample performance with controllable level of
% conservativeness
figure;
% Line 7
subplot(3,2,1)
plot(rho_Matrix_0, LineFlow_0(:,7),'o-','LineWidth',1.5)
hold on
plot(rho_Matrix_1, LineFlow_1(:,7), '*-','LineWidth',1.5)
hold on
plot(rho_Matrix_2, LineFlow_2(:,7), 's-','LineWidth', 2)
hold on 
plot(rho_Matrix_0,-600.*ones(32,1),'k--','LineWidth',2,'color',[0.5 0.5 0.5]) % nominal power limitations
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'], 'nominal flow limit')
xlabel('weight factor \rho')
ylabel('Line flow (MW)')
title('Line 7')
set(gca,'FontSize',13);
grid

% Line 37
subplot(3,2,2)
plot(rho_Matrix_0, LineFlow_0(:,37),'o-','LineWidth',1.5)
hold on
plot(rho_Matrix_1, LineFlow_1(:,37), '*-','LineWidth',1.5)
hold on
plot(rho_Matrix_2, LineFlow_2(:,37), 's-','LineWidth', 2)
hold on 
plot(rho_Matrix_0, 500.*ones(32,1),'k--','LineWidth',2,'color',[0.5 0.5 0.5]) % nominal power limitations
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'], 'nominal flow limit')
xlabel('weight factor \rho')
ylabel('Line flow (MW)')
title('Line 37')
set(gca,'FontSize',13);
grid

% Line 38
subplot(3,2,3)
plot(rho_Matrix_0, LineFlow_0(:,38),'o-','LineWidth',1.5)
hold on
plot(rho_Matrix_1, LineFlow_1(:,38), '*-','LineWidth',1.5)
hold on
plot(rho_Matrix_2, LineFlow_2(:,38), 's-','LineWidth', 2)
hold on
plot(rho_Matrix_0, 500.*ones(32,1),'k--','LineWidth',2,'color',[0.5 0.5 0.5]) % nominal power limitations
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'], 'nominal flow limit')
xlabel('weight factor \rho')
ylabel('Line flow (MW)')
title('Line 38')
set(gca,'FontSize',13);
axis([0 1000 350 550])
grid


% Line 54
subplot(3,2,4)
plot(rho_Matrix_0, LineFlow_0(:,54),'o-','LineWidth',1.5)
hold on
plot(rho_Matrix_1, LineFlow_1(:,54), '*-','LineWidth',1.5)
hold on
plot(rho_Matrix_2, LineFlow_2(:,54), 's-','LineWidth', 2)
hold on
plot(rho_Matrix_0, 500.*ones(32,1),'k--','LineWidth',2,'color',[0.5 0.5 0.5]) % nominal power limitations
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'], 'nominal flow limit')
xlabel('weight factor \rho')
ylabel('Line flow (MW)')
title('Line 54')
set(gca,'FontSize',13);
axis([0 1000, 450 800])
grid 

% Line 96
subplot(3,2,5)
plot(rho_Matrix_0, LineFlow_0(:,96),'o-','LineWidth',1.5)
hold on
plot(rho_Matrix_1, LineFlow_1(:,96), '*-','LineWidth',1.5)
hold on
plot(rho_Matrix_2, LineFlow_2(:,96), 's-','LineWidth', 2)
hold on
plot(rho_Matrix_0, 300.*ones(32,1),'k--','LineWidth',2,'color',[0.5 0.5 0.5]) % nominal power limitations
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'], 'nominal flow limit' )
xlabel('weight factor \rho')
ylabel('Line flow (MW)')
title('Line 96')
set(gca,'FontSize',13);
grid

%% 3-rd plot
% This plot visualizes the policies and output powers of selected
% generators for various values of risk aversion and Wasserstein metric

figure;
% select the index of Generator you want to show the policy and power outputs;
Gen1_Index = 4;
Gen2_Index = 5;
Gen3_Index = 35;
Gen4_Index = 49;

% The first column,  power output adjustment e
subplot(4,4,1)
plot(rho_Matrix_0,e_Data_0(:,Gen1_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,e_Data_1(:,Gen1_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,e_Data_2(:,Gen1_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('power output (MW)')
title(['#',num2str(Gen1_Index),' Generator at bus ', num2str(mpc.gen(Gen1_Index,1))])
grid

subplot(4,4,5)
plot(rho_Matrix_0,e_Data_0(:,Gen2_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,e_Data_1(:,Gen2_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,e_Data_2(:,Gen2_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('power output (MW)')
title(['#',num2str(Gen2_Index),' Generator at bus ', num2str(mpc.gen(Gen2_Index,1))])
grid

subplot(4,4,9)
plot(rho_Matrix_0,e_Data_0(:,Gen3_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,e_Data_1(:,Gen3_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,e_Data_2(:,Gen3_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('power output (MW)')
title(['#',num2str(Gen3_Index),' Generator at bus ', num2str(mpc.gen(Gen3_Index,1))])
grid

subplot(4,4,13)
plot(rho_Matrix_0,e_Data_0(:,Gen4_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,e_Data_1(:,Gen4_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,e_Data_2(:,Gen4_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('power output (MW)')
title(['#',num2str(Gen4_Index),' Generator at bus ', num2str(mpc.gen(Gen4_Index,1))])
grid

% The second column, the reserve policy D1 
subplot(4,4,2)
plot(rho_Matrix_0,D1_Data_0(:,Gen1_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,D1_Data_1(:,Gen1_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,D1_Data_2(:,Gen1_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('D1')
title(['#',num2str(Gen1_Index),' Generator at bus ', num2str(mpc.gen(Gen1_Index,1))])
grid

subplot(4,4,6)
plot(rho_Matrix_0,D1_Data_0(:,Gen2_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,D1_Data_1(:,Gen2_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,D1_Data_2(:,Gen2_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('D1')
title(['#',num2str(Gen2_Index),' Generator at bus ', num2str(mpc.gen(Gen2_Index,1))])
grid

subplot(4,4,10)
plot(rho_Matrix_0,D1_Data_0(:,Gen3_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,D1_Data_1(:,Gen3_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,D1_Data_2(:,Gen3_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('D1')
title(['#',num2str(Gen3_Index),' Generator at bus ', num2str(mpc.gen(Gen3_Index,1))])
grid

subplot(4,4,14)
plot(rho_Matrix_0,D1_Data_0(:,Gen4_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,D1_Data_1(:,Gen4_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,D1_Data_2(:,Gen4_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('D1')
title(['#',num2str(Gen4_Index),' Generator at bus ', num2str(mpc.gen(Gen4_Index,1))])
grid

% The third column, the reserve policy D2 
subplot(4,4,3)
plot(rho_Matrix_0,D2_Data_0(:,Gen1_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,D2_Data_1(:,Gen1_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,D2_Data_2(:,Gen1_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('D2')
title(['#',num2str(Gen1_Index),' Generator at bus ', num2str(mpc.gen(Gen1_Index,1))])
grid

subplot(4,4,7)
plot(rho_Matrix_0,D2_Data_0(:,Gen2_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,D2_Data_1(:,Gen2_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,D2_Data_2(:,Gen2_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('D2')
title(['#',num2str(Gen2_Index),' Generator at bus ', num2str(mpc.gen(Gen2_Index,1))])
grid

subplot(4,4,11)
plot(rho_Matrix_0,D2_Data_0(:,Gen3_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,D2_Data_1(:,Gen3_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,D2_Data_2(:,Gen3_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('D2')
title(['#',num2str(Gen3_Index),' Generator at bus ', num2str(mpc.gen(Gen3_Index,1))])
grid

subplot(4,4,15)
plot(rho_Matrix_0,D2_Data_0(:,Gen4_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,D2_Data_1(:,Gen4_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,D2_Data_2(:,Gen4_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('D2')
title(['#',num2str(Gen4_Index),' Generator at bus ', num2str(mpc.gen(Gen4_Index,1))])
grid

% The fourth column, the reserve policy D3 
subplot(4,4,4)
plot(rho_Matrix_0,D3_Data_0(:,Gen1_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,D3_Data_1(:,Gen1_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,D3_Data_2(:,Gen1_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('D3')
title(['#',num2str(Gen1_Index),' Generator at bus ', num2str(mpc.gen(Gen1_Index,1))])
grid

subplot(4,4,8)
plot(rho_Matrix_0,D3_Data_0(:,Gen2_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,D3_Data_1(:,Gen2_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,D3_Data_2(:,Gen2_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('D3')
title(['#',num2str(Gen2_Index),' Generator at bus ', num2str(mpc.gen(Gen2_Index,1))])
grid

subplot(4,4,12)
plot(rho_Matrix_0,D3_Data_0(:,Gen3_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,D3_Data_1(:,Gen3_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,D3_Data_2(:,Gen3_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('D3')
title(['#',num2str(Gen3_Index),' Generator at bus ', num2str(mpc.gen(Gen3_Index,1))])
grid

subplot(4,4,16)
plot(rho_Matrix_0,D3_Data_0(:,Gen4_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_1,D3_Data_1(:,Gen4_Index),'LineWidth',1.5)
hold on
plot(rho_Matrix_2,D3_Data_2(:,Gen4_Index),'LineWidth',1.5)
legend([char(949),'  = 0.00'],[char(949), '  = 0.02'], [char(949), '  = 0.04'])
xlabel('weight factor \rho')
ylabel('D3')
title(['#',num2str(Gen4_Index),' Generator at bus ', num2str(mpc.gen(Gen4_Index,1))])
grid

disp('Finished!')

toc







