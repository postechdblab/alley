# Alley

The source code of "Combining Sampling and Synopses with Worst-Case Optimal Runtime and Quality Guarantees for Graph Pattern Cardinality Estimation" (SIGMOD 2021, [link](https://dl.acm.org/doi/10.1145/3448016.3457246)) on top of G-CARE ([link](https://dl.acm.org/doi/10.1145/3318464.3389702)).

Technical Paper ([link](https://arxiv.org/abs/2103.13681)).

Please see the detailed instructions of running G-CARE first ([link](https://github.com/yspark-dblab/gcare)).


## A short README

### 1. Installation

```
  mkdir build
  cd build
  cmake ..
  make
```

### 2. Build mode (convert data graph in txt file to binary & build synopses)

```
  <bin> -b -m <method> -i <input_data_graph> -d <output_data_header> -o <output_file>
  
  <bin>: binary to execute (gcare_graph, gcare_relation, gcare_inter)
  <method>: method (alley, alleyTPI, wj, ibjs, cs, sumrdf, cset, impr, bsk, jsub)
  <input_data_graph>: input data graph txt file
  <output_data_header>: output data header, files named <output_data_header>.* will be created
  <output_file>: log reporting the elapsed time
```

### 3. Query mode (run estimation)
  
``` 
  <bin> -q -m <method> -i <input_query_dir> -d <data_header> -p <sampling_ratio> -n <num_repeat> -o <output_file>
  
  <input_query_dir>: input query directory, all .txt files are searched recursively & regared as query graphs
  <data_header>: corresponds to <output_data_header> in build mode
  <sampling_ratio>: sampling ratio (e.g., 0.01 denotes 1%, 0.001 denotes 0.1% (default)) 
  <num_repeat>: number of repetition for each query
  <output_file>: log reporting individual estimates and average elapsed time (per query)
```
