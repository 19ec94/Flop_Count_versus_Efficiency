# Flop_Count_versus_Efficiency

This repository contains programming solutions developed as part of the *Simulation Science Seminar* coursework at SiSc, RWTH Aachen University. The primary focus is on analyzing the correlation between floating-point operations (FLOPs) and computational efficiency across various numerical methods.

## 📚 Overview

The project aims to provide insights into how the number of FLOPs influences the performance of different numerical algorithms. By examining metrics such as FLOP counts and efficiency ratios, the repository serves as a resource for anyone interested in computational performance analysis.

## ⚙️ Features

- **FLOP Count Analysis** – Detailed computation of floating-point operations involved in various numerical methods.
- **Efficiency Metrics** – Evaluation of computational efficiency based on FLOP counts and execution time.
- **Visualization** – Graphs and charts illustrating the relationship between FLOP counts and efficiency.
- **Comparative Studies** – Side-by-side comparisons of different algorithms in terms of performance.

## 🔧 Installation

Clone the repository to your local machine:

```bash
git clone https://github.com/19ec94/Flop_Count_versus_Efficiency.git
cd Flop_Count_versus_Efficiency/src
```
## Usage
Tha main program file contains 5 matrices multiplication. It uses OpenBlas for matrix computation.

For ease of compilation a makefile is provided. The given program uses OpenBLAS libraray. 
Before executing the makefile make sure to link the OpenBLAS library path in your system.
For example-see makefile- (in my case) 
```bash
gcc -o output main.o multi.o -I /opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas -lpthread -lgfortran   
```
make sure to provide your system path
```bash
gcc -o output main.o multi.o -I /<your system path>/OpenBLAS/include/ -L/<your system path>/OpenBLAS/lib -lopenblas -lpthread -lgfortran
```

if it is done, then run the follwing commands

```bash
>>> make
>>>./output arg1 arg2 arg3 arg4 arg5 arg6
```

you should see some output on your screen at this point :).

