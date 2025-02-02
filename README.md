# Gene Motif Analysis

This project analyzes gene promoter regions and classifies motifs based on their length. It processes CSV files containing motif data and outputs aggregated results.

## 📂 Project Structure
```
project_root/
│── data/                          # Contains input CSV files
│   ├── promoter_gene_ranges_within_2000_genes_lab_new.csv
│   ├── promoter_gene_ranges_within_2000_genes_random.csv
│── src/                           # Contains Python scripts
│   ├── finder_list.py              # Main script for processing data
│── results/                        # Stores processed output
│── README.md                       # Project documentation
│── environment.yml               # Dependencies
│── .gitignore                      # Ignore unnecessary files
```

## 📦 Installation

Ensure you have Python installed. Install required dependencies:
```bash
conda env create -f environment.yml
```

## 🚀 Usage

To process a CSV file, run:
```bash
python Project_root/finder_list.py data/input.csv results/output.csv
```

## 📊 Expected Output
The script processes motif data and outputs an aggregated CSV file containing summed nucleotide counts (A, T, C, G).

## 📝 Author
Yi-Shan Lan
