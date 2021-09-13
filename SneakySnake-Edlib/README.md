# SneakySnake-Edlib :  Integrating SneakySnake and Edlib
This directory contains the source codes and results of integrating SneakySnake and Edlib. By integrating SneakySnake and Edlib, it is aimed to achieve such that sequence alignment leverages the data structure calculated by SneakySnake for fast and efficient read mapping.

#### :zap: Edlib is a open-source C/C++ library for sequence alignment using edit distance. For more information, please visit link [here.](https://github.com/Martinsos/edlib)

#### :zap: If you want to learn more about SneakySnake, please go back to [root directory.](https://github.com/ilknurbas/SneakySnake)

## Table of Contents
- [Algorithm](#logic)
- [Directory Structure](#directory)
- [Getting Started](#install)


##  <a name="logic"></a>Algorithm
Explain algo

##  <a name="directory"></a>Directory Structure:
```
SneakySnake-Edlib
├───1. Evaluation Results
    ├───1.1 
    ├───1.2 
├───2.  
├───3. kthread.c
├───4. kthread.h
├───5. main.c
├───6. Makefile
├───7. README.md
├───8. SneakySnake.c
├───9. SneakySnake.h
```      
:information_source: #1 is the directory that contains all evaluations and outputs using datasets from "../Datasets". <br />
:information_source: #2 is the directory that contains ??  <br />
:information_source: #5 is the C file where the logic of integrating SneakySnake and Edlib is implemented. <br />
:information_source: #8 and #9 are the C files that contains the original SneakySnake algorithm with modifications.


## <a name="install"></a>Getting Started

### <a name="run"></a>Installation
```sh
git clone https://github.com/ilknurbas/SneakySnake
cd SneakySnake-Edlib && make

./main [DebugMode] [KmerSize] [ReadLength] [IterationNo] [ReadRefFile] [# of reads] [# of threads] [EditThreshold]
# Short sequences
./main 0 100 100 100 ../Datasets/ERR240727_1_E2_30000Pairs.txt 30000 10 10
# Long sequences
./main 0 20500 100000 20500 ../Datasets/LongSequences_100K_PBSIM_10Pairs.txt 10 40 20000
```
#### :exclamation: updating...