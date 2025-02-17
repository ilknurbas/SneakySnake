
## <a name="test"></a>Evaluation Results
* The Excel document contains all evaluation results using the datasets in **"../../Datasets"** directory. </br>
* **"Files"** directory contains all .txt files that are generated for each evaluation case in the Excel document.

## <a name="output"></a>Output
Run 
```sh
./main 0 100 100 100 ../Datasets/ERR240727_1_E2_30000Pairs.txt 100 4 20
```
and as a output, **File.txt** will be generated inside SneakySnake-Edlib directory.
 
```sh
Seq No      X      Y    X-Y
     0     -1     -1     -1
     1     -1     -1     -1
     2     -1     -1     -1
     3     -1     -1     -1
     4     -1     -1     -1
     5     -1     -1     -1
     6     -1     -1     -1
     7     -1     -1     -1
     8     -1     -1     -1
     9     -1     -1     -1
    10     -1     -1     -1
    11     -1     19     20
    12     -1     -1     -1
    13     -1     -1     -1
    14     -1     -1     -1
    15      0     -1     -1
    16     -1     -1     -1
    17     18     -1    -19
    18     -1     -1     -1
    19     20     -1    -21
    20     -1     -1     -1
    21      3      2      1
    22     -1     -1     -1
    23     15     -1    -16
    24     11      9      2
    25     17     -1    -18
    26     -1     -1     -1
    27     -1     -1     -1
    28     -1     -1     -1
    29     15     19     -4
    30     -1     -1     -1
    31     17     -1    -18
    32     20     -1    -21
    33     18     -1    -19
    34     20     -1    -21
    35     18     -1    -19
    36     19     16      3
    37     19     18      1
    38     -1     -1     -1
    39     -1     -1     -1
    40     19     -1    -20
    41     -1     -1     -1
    42     -1     -1     -1
    43     -1     -1     -1
    44     -1     -1     -1
    45      2      1      1
    46      1      1      0
    47     16     -1    -17
    48     10     12     -2
    49     14      8      6
    50     13     10      3
    51     20     -1    -21
    52     11     10      1
    53     19     -1    -20
    54     10     15     -5
    55     18     -1    -19
    56     18     -1    -19
    57     18     17      1
    58     12     12      0
    59     13      8      5
    60     14     -1    -15
    61     16     -1    -17
    62      8      1      7
    63      1      2     -1
    64      9      4      5
    65     -1     -1     -1
    66     -1     -1     -1
    67     16     18     -2
    68     -1     -1     -1
    69     18     -1    -19
    70     -1     -1     -1
    71     -1     -1     -1
    72     -1     -1     -1
    73     -1      2      3
    74     -1     -1     -1
    75     -1     -1     -1
    76     -1     -1     -1
    77     -1     -1     -1
    78     -1     -1     -1
    79     14     15     -1
    80     18     -1    -19
    81     20     -1    -21
    82      0     -1     -1
    83     -1     -1     -1
    84     -1     -1     -1
    85     20     -1    -21
    86     -1     -1     -1
    87     -1     -1     -1
    88     -1     -1     -1
    89     -1     -1     -1
    90      5      3      2
    91     -1     19     20
    92     12     15     -3
    93     12      9      3
    94     12      7      5
    95     -1     -1     -1
    96     14     10      4
    97      9      8      1
    98     13     11      2
    99     17     15      2
For ../Datasets/ERR240727_1_E2_30000Pairs.txt
 #positive 22
 #neutral 2
 #negative 76

```
Also, on console the following information is printed.
```sh
This system has 4 processors configured and 4 processors available.
        Sneaky Snake and Edlib Total Time: 7 milliseconds
Seq No       X      Y    X-Y
For ../Datasets/ERR240727_1_E2_30000Pairs.txt
 #positive 21
 #neutral 2
 #negative 77
```
