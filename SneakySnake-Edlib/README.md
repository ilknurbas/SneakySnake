# SneakySnake-Edlib :  Integrating SneakySnake and Edlib
This directory contains the source codes and results of integrating SneakySnake and Edlib. By integrating SneakySnake and Edlib, it is aimed to achieve such that sequence alignment leverages the data structure calculated by SneakySnake for fast and efficient read mapping.

#### :zap: Edlib is a open-source C/C++ library for sequence alignment using edit distance. For more information, please visit link [here.](https://github.com/Martinsos/edlib)

#### :zap: If you want to learn more about SneakySnake, please go back to [root directory.](https://github.com/ilknurbas/SneakySnake)

## Table of Contents
- [Algorithm](#logic)
- [Directory Structure](#directory)
- [Get Started](#run)

##  <a name="logic"></a>Algorithm 
SneakySnake checks all rows and finds the row with the largest number of zeros. Once the checkpoint is located, the cigar string coming from SneakySnake is also generated. SneakySnake is done at this point. After that, the Edlib function is called and takes a sequence of length 11 **(refer to line 154 in main.c)** as an input having the character located in the checkpoint is in the middle. Then, Edlib returns a cigar string. The program excludes the region that is sent to Edlib from the cigar string that SneakySnake generates and merges this cigar string with the one coming from Edlib. The steps above are repeated until the end of the sequence is reached.  


##  <a name="directory"></a>Directory Structure:
```
SneakySnake-Edlib
├───1. Evaluation Results
    ├───1.1 
    ├───1.2   
├───2. Makefile
├───3. README.md
├───4. SneakySnake.c
├───5. SneakySnake.h
├───6. kthread.c
├───7. kthread.h
├───8. main.c

```      
:information_source: #1 is the directory that contains all evaluations and outputs using datasets from "../Datasets". <br />
:information_source: #4 and #5 are the C files that contains the original SneakySnake algorithm with some modifications. <br />
:information_source: #8 is the C file where the logic of integrating SneakySnake and Edlib is implemented.

## <a name="run"></a> Get Started
```sh
git clone https://github.com/ilknurbas/SneakySnake
cd SneakySnake-Edlib && make

./main [DebugMode] [KmerSize] [ReadLength] [IterationNo] [ReadRefFile] [# of reads] [# of threads] [EditThreshold]
# Short sequences
./main 0 100 100 100 ../Datasets/ERR240727_1_E2_30000Pairs.txt 30000 10 10
# Long sequences
./main 0 20500 100000 20500 ../Datasets/LongSequences_100K_PBSIM_10Pairs.txt 10 40 20000
```
#### :exclamation: Set "DebugMode" paramater to 0 all the time.
