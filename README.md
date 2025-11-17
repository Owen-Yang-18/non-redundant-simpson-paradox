# Finding Non-redundant Simpson's Paradox

A Rust implementation for detecting non-redundant Simpson's Paradox in multi-dimensional datasets.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Input Format](#input-format)
- [Output Format](#output-format)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Understanding Results](#understanding-results)
- [Citation](#citation)


## Prerequisites

- Rust 1.70 or later (with Cargo)
- Dependencies (specified in `Cargo.toml`):
  - `rayon` - for parallel computation

## Installation

### Step 1: Install Rust

If you don't have Rust installed, get it from [rustup.rs](https://rustup.rs/):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### Step 2: Project Directory

```bash
# Navigate to code directory
cd simpson-rust
```

### Step 3: Build

```bash
# Build in release mode for best performance
cargo build --release
```

The executable will be in `target/release/simpson-rust`.

## Usage

### Basic Command

```bash
cargo run --release <SIZE> <NUM_LABELS> <FILENAME> [THRESHOLD]
```

### Parameters

- **SIZE**: Number of categorical attributes
- **NUM_LABELS**: Number of label attributes
- **FILENAME**: Path to the input CSV file
- **THRESHOLD**: (Optional) Minimum size threshold for population coverage (default: 0)

### Example Commands

```bash
# Basic usage with adult dataset
cargo run --release 8 1 adult.csv

# With threshold to filter population with limited coverage
cargo run --release 8 1 adult.csv 100

# With a different dataset
cargo run --release 12 1 loan.csv 100
```

## Input Format

The input must be a CSV file with the following structure:

### Format Requirements

1. **Header Row**: First row contains categorical attribute and label names
2. **Categorical Attribute Columns**: First `SIZE` columns (from the left) are categorical attributes
3. **Label Columns**: Next `NUM_LABELS` columns are binary labels

### Example: Adult Income Dataset

| workclass | education | m-status | occupation | relationship | race | sex | native-country | income |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| State-gov | Bachelors | Never-married | Adm-clerical | Not-in-family | White | Male | United-States | 0 |
| Self-emp-not-inc | Bachelors | Married-civ-spouse | Exec-managerial | Husband | White | Male | United-States | 0 |
| Private | HS-grad | Divorced | Handlers-cleaners | Not-in-family | White | Male | United-States | 0 |
| Private | 11th | Married-civ-spouse | Handlers-cleaners | Husband | Black | Male | United-States | 0 |


Where:
- First 8 columns: Categorical Attributes (demographic and work attributes)
- Last column: Income label (1 for >50K, 0 for <=50K)

## Output Format

The program generates two types of output:

### Console Output

```
Dataset: adult.csv
Number of records: 32561
Agg size: 8945
Total time for materialization: 5.6789 seconds
Total number of non-empty populations: 15678
Total number of coverage equivalent subsets: 892
Total runtime: 8.4532 seconds
Total number of SP: 234
Total number of redundant relations: 52
Results written to adult_infos.txt
```

### Output File (`<dataset>_infos.txt`)

The program creates a detailed text file with all detected paradoxes:

```
Total redundant paradox groups: 52
Showing top 52 groups (by size)
======================================================================================

Group 1 (Signature: 1234567890, Count: 15)
--------------------------------------------------------------------------------------
  Paradox 1: s1=(relationship=Not-in-family, workclass=*, education=Bachelors, ...), 
             s2=(relationship=Husband, workclass=*, education=Bachelors, ...), 
             separator=sex, label=income
  Paradox 2: s1=(relationship=Not-in-family, workclass=Private, education=*, ...), 
             s2=(relationship=Husband, workclass=Private, education=*, ...), 
             separator=sex, label=income
...
```

### Understanding the Output

Each paradox entry shows:
- **s1, s2**: Two sibling populations being compared
- **separator**: The separator attribute for population partitioning
- **label**: The label attribute for measuring frequency statistics

## Examples

### Example 2: Adult Income Dataset

Predict whether income exceeds $50K/year based on census data.

**Dataset**: `adult.csv` (available from [UCI Machine Learning Repository](https://archive.ics.uci.edu/ml/datasets/adult))

```bash
# Download dataset
wget https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data
mv adult.data adult.csv

# Run detection
# Please remove continuous attributes first
cargo run --release 8 1 adult.csv 100
```

This uses a threshold of 100, meaning only populations covering at least 100 records are materialized.

## Troubleshooting

### Common Issues

#### 1. Out of Memory
```
error: memory allocation failed
```

**Solution**: Increase threshold or reduce dataset size:
```bash
# Use higher threshold
cargo run --release 10 1 data.csv 100

# Or sample the data first
head -n 10000 data.csv > data_sample.csv
```

#### 2. Compilation Errors
```
error: could not compile `simpson_paradox_detector`
```

**Solution**:
```bash
# Update Rust
rustup update

# Clean and rebuild
cargo clean
cargo build --release
```

#### 3. Slow Performance

**Optimization Tips**:
```bash
# Ensure release mode
cargo build --release  # NOT cargo build

# Increase threshold
cargo run --release 10 1 data.csv 100

# Check CPU usage
top  # Should see near 100% * num_cores
```

## Understanding Results

### Interpreting Paradoxes

When you find a paradox like:

```
s1=(age=young, gender=M, treatment=*, outcome=*), 
s2=(age=old, gender=M, treatment=*, outcome=*), 
separator=smoking, label=recovery
```

This means:

1. **Overall**: Young males and old males have similar recovery rates
2. **When partitioned by smoking**:
   - Among smokers: Young males have better recovery than old males
   - Among non-smokers: Young males have better recovery than old males
3. **Paradox**: The relationship appears different when you consider smoking status

### Understanding Redundancy

The tool identifies and groups **redundant paradoxes**: multiple paradoxes that describe the same underlying phenomenon. Redundancy arises in three ways:

#### 1. Sibling Child Equivalence

This occurs when different populations cover the same set of records. For example:

```
Paradox 1: s1=(age=young, workclass=*, education=Bachelors, race=*)
           s2=(age=old, workclass=*, education=Bachelors, race=*)
           
Paradox 2: s1=(age=young, hours=40, education=Bachelors, race=*)
           s2=(age=old, hours=40, education=Bachelors, race=*)
```

If `workclass=*` and `hours=40` cover exactly the same records, these paradoxes are redundant; they're comparing the same groups of people, just described differently.

**Why it happens**: In sparse datasets, different attribute combinations can accidentally cover identical record sets.

#### 2. Separator Equivalence

This occurs when different separator attributes partition the data identically. For example:

```
Paradox 1: separator=marital_status (partitions: married, single, divorced)
Paradox 2: separator=spouse_present (partitions: yes, no, N/A)
```

If there's a one-to-one correspondence (married↔yes, single↔no, divorced↔N/A) and each pair covers the same records, these separators reveal the same structure.

**Why it happens**: Correlated or dependent attributes in the dataset create redundant ways to partition the data.

#### 3. Statistic Equivalence

This occurs when different label attributes exhibit the same statistical behavior. For example, if two health outcomes are perfectly correlated in your dataset, any paradox found with one label will also exist with the other.

**Why it happens**: Dependent or highly correlated outcome variables in the data.

### Redundancy Groups

The tool groups redundant paradoxes together and reports them as a single **redundant paradox group**. Each group is characterized by:

- **E1 × E2**: Sets of equivalent sibling populations being compared
- **X**: Set of equivalent separator attributes
- **Y**: Set of equivalent label attributes

For instance, a group with 8 paradoxes might represent:
- 2 equivalent sibling populations (s1, s2)
- 2 equivalent separator attributes
- 2 equivalent label attributes

Total: 2 × 2 × 2 = 8 redundant paradoxes, all describing the same underlying phenomenon.


## Citation

If you use find this project interesting or helpful for your research, please consider citing our paper:

```bibtex
@article{yang2025finding,
  title={Finding Non-Redundant Simpson's Paradox from Multidimensional Data},
  author={Yang, Yi and Pei, Jian and Yang, Jun and Xie, Jichun},
  journal={arXiv preprint arXiv:2511.00748},
  year={2025}
}
```