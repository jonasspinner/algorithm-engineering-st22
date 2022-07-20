# Algorithm Engineering 2022
This repository contains the framework for the exercises regarding the Algorithm Engineering 2022 lecture.

There are three parts to this repository:
<ul>
    <li><code>benchmark</code> provides code for benchmarking your implementations. 
        <code>mst_construction_benchmark.cpp</code> and <code>mst_verification_benchmark.cpp</code> contain 
        the program entries for benchmarking MST construction and MST verification algorithms, respectively.</li>
    <li><code>src</code> is for your code. Any file you place here will be picked up by the CMake build 
        system automatically, although you still need to register your implementation in the corresponding 
        <code>*_parameters.hpp</code> file. The given example contenders indicate how to register your 
        implementation. </li>
    <li><code>binaries</code> contains a library as a binary file <code>libbinaries-*.a</code>/<code>binaries-*.lib</code> which provides methods 
        used by the framework as defined in <code>includes/binary_includes.hpp</code>. You are not allowed 
        to call these methods in your implementations.</li>
</ul>


## Usage
You need cmake (version 3.16+) and a C++17-capable compiler to build this project.
We strongly encourage using a Linux system to develop your code.

(If you are on Windows, you can also look into Linux virtual machines.)
Our binary library currently only supports GCC (g++) version 10+, Clang version 10+ 
and MSVC version 14+. If you use Windows, make sure you use the latest version 
of VS 2022 Community Edition. If you absolutely need to use a different (publicly 
available) compiler, we may be able to provide another library binary file to 
use with that compiler. In that case, please contact us with the compiler 
you wish to use.

CMake helps build and package your code.
If you're unfamiliar, this might help you get started (on UNIX):

``` shell
mkdir build && cd build
cmake ..
make
```


Add your own code in <code>src</code> and register your implementation 
in [`mst_construction_parameters.hpp`](src/mst_construction_parameters.hpp) 
or [`mst_verification_parameters.hpp`](src/mst_verification_parameters.hpp)
, as appropriate.
You're free to remove the example implementations or adjust any of the 
other parameters of course. In particular, you may want to adjust the 
size and density of the generated random input graph for your tests. 
(It is advised to disable the naive example implementations, e.g. 
NaiveKruskal, for larger graphs.)

You shouldn't have to touch any code in <code>benchmark</code>. 
If you want to adapt the provided implementations, copy them to <code>src</code> 
and place them into a custom namespace (see below).
This way, you'll still have the baseline implementation as a reference for your evaluation.
If you feel you need to make changes to the benchmarking or reference code, let us know.

If you work on "Aufgabe 4: Linear-Time MST Algorithmus mit Blackbox", you may use 
<code>algen::EdgeClassifier</code> as the black box for determining edges that are light wrt. to a sample MST.
Please acquire the singleton instance using <code>algen::getEdgeClassifier()</code> as the framework uses it to measure
the time spent in the black box. The framework subtracts this time from the total running time of your algorithm.
The class may not be used for any task other than Aufgabe 4.

## Requirements
In order to work with the included benchmarks, your code will need to provide the expected API.
You can find the specifics in the <code>*_parameters.hpp</code> files.

You should place all your code in a namespace. You can use your own name, a name describing 
your implementation or any other name you can reasonably expect to be unique. 
This is so we can merge your code with that of other students without name clashes.
(It is strongly forbidden to use top-level <code>using namespace</code> in your headers as that defeats the point
of namespaces all together.)

## Plots
You find some suggestions/help for plots in <code>plots</code>. You can execute the plotting scripts as follows:
``` shell
Rscript ./evaluation_mst_[construction,verification].r --infile <path to .csv file> --outfile <name of pdf file with plots>
```

## Submission of your code
To submit the final version of your project (deadline will be at the beginning of september) you can simply send us a link to the repository (on github, scc-gitlab, etc.) into which you have forked this framework + the final commit hash/release via email.

