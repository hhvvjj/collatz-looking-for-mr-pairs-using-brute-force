# Looking for mr Pairs Using Brute Force

A high-performance parallel C program for discovering unique `mr` values in Collatz sequences using the tuple-based transform.

## Features

This program analyzes Collatz sequences to find the first repeated value in the transformed sequence of `m` parameters. For each positive integer `n`, it computes the Collatz sequence and applies the transformation `m = (c - p) / 2` where `p = 2` if `c` is even, `p = 1` if `c` is odd. The first repeated `m` value in this sequence is called `mr`.

- **High Performance**: Optimized parallel implementation using OpenMP
- **Real-time Discovery**: Shows unique `mr` values as they are discovered
- **Adaptive Scheduling**: Automatically chooses optimal work distribution strategy
- **Memory Efficient**: Dynamic memory allocation with intelligent growth
- **Large Scale**: Supports analysis up to 2^64 numbers
- **Thread Safe**: Lock-optimized concurrent access to shared data structures

## Dependencies

### Ubuntu/Debian
```bash
sudo apt-get install gcc libomp-dev
```

### CentOS/RHEL
```bash
sudo yum install gcc libgomp-devel
```

## Installation

```bash
# Clone the repository
git clone https://github.com/hhvvjj/collatz-looking-for-mr-pairs-using-brute-force.git
cd collatz-looking-for-mr-pairs-using-brute-force

# Compile 
gcc -fopenmp -O3 -march=native looking_for_mr_pairs_using_brute_force.c -lgomp -lpthread -o looking_for_mr_pairs_using_brute_force
```

## Usage

```bash
./looking_for_mr_pairs_using_brute_force <exponent>
```

### Examples

```bash
# Quick test - analyze numbers 1 to 2^20 (1M numbers)
./looking_for_mr_pairs_using_brute_force 20

# Standard analysis - numbers 1 to 2^25 (33M numbers)  
./looking_for_mr_pairs_using_brute_force 25

# Large scale - numbers 1 to 2^30 (1B numbers)
./looking_for_mr_pairs_using_brute_force 30
```

## Performance Guide

| Exponent | Range | Numbers | Est. Time* | Memory Usage |      Use Cases     |
|----------|-------|---------|------------|--------------|--------------------|
|    15    |  2^15 |   32K   |    < 1 sec |    < 1 MB    | Testing            |
|    20    |  2^20 |    1M   |     ~5 sec |    ~10 MB    | Quick analysis     |
|    25    |  2^25 |   33M   |   ~2-5 min |    ~50 MB    | Standard research  |
|    30    |  2^30 |    1B   | ~1-3 hours |   ~200 MB    | Intensive analysis |
|    35    |  2^35 |   34B   |  ~1-2 days |     ~1 GB    | Extreme analysis   |

*Times are approximate and depend on hardware (8-core CPU assumed)

## Output

### Real-time Progress
```
Using scheduling strategy: guided (medium range)
 [!] NEW UNIQUE mr = 0 found at n = 1 (total unique: 1)
 [!] NEW UNIQUE mr = 2 found at n = 3 (total unique: 2)
 [!] NEW UNIQUE mr = 6 found at n = 7 (total unique: 3)
 ...
Progress: (31.2%) | Processed: 163354 | Unique values of mr found: 42 | 429496.7 nums/sec | ETA: 0.0 min
```

### Final Results
```
Final Summary:

 [*] Exponent: 19
 [*] Total time: 1.10 seconds
 [*] Speed: 476191.9 numbers/second
 [*] Numbers processed: 2^19 = 524287
 [*] Numbers with mr found: 524287
 [*] Percentage with mr: 100.00%
 [*] Unique mr values found: 42

Final list of mr values:
0, 1, 2, 3, 6, 7, 8, 9, 12, 16, 19, 25, 45, 53, 60, 79, 91, 121, 125, 141, 166, 188, 205, 243, 250, 324, 333, 432, 444, 487, 576, 592, 649, 667, 683, 865, 889, 1153, 1214, 1821, 2428, 3643

Results saved to `mr_results_2_19.txt`:
```

## Research Background

By analyzing the m parameters derived from each step using the tuple-based transform on a Collatz sequence, we can identify patterns in how numbers "approach" the final cycle. The mr values represent the first repeated elements in these sequences, providing insight into the internal structure of Collatz dynamics.
Only a finite number of unique mr values exist, with all discoveries occurring within relatively small bounds (n ≤ 7287). This finding has implications for understanding the global behavior of Collatz sequences and provides pathways towards convergence patterns.
This brute-force approach serves as both a verification tool for theoretical bounds and a discovery mechanism for identifying the complete set of possible mr values across different ranges of starting values of n.

**Reference:** The complete theoretical framework is detailed in my research paper: http://dx.doi.org/10.5281/zenodo.15546925

## Algorithm Details

### Core Algorithm
1. **Sequence Generation**: Standard 3n+1 Collatz sequence
2. **Transformation**: Apply m = (c-p)/2 for each step
3. **Repetition Detection**: Track all m values until first repetition
4. **Optimization**: Early termination for long sequences likely going to trivial cycle

### Parallel Strategy
- **Work Distribution**: Adaptive scheduling (static/guided/dynamic)
- **Thread Safety**: Lock-optimized unique value detection
- **Memory Management**: Dynamic allocation with intelligent growth
- **Progress Tracking**: Real-time discovery reporting

## Contributing

Contributions are welcome! Areas of interest:

  - Extended validation testing using large exponents
  - Documentation improvements

## Academic Citation

For academic use, please cite both the original theoretical work and this implementation.

## License

CC-BY-NC-SA 4 License - see LICENSE file for details.

## Contact

For questions about the algorithm implementation or mathematical details drop me an email.
