# Data Analysis

This folder contains a Jupyter notebook for analyzing search algorithm performance data.

## Prerequisites

You'll need [uv](https://docs.astral.sh/uv/) installed to manage Python dependencies and run the notebook.

## Running the Notebook

1. Navigate to this directory
2. Start Jupyter:
   ```bash
   uv run jupyter notebook analytics.ipynb
   ```

This will automatically install all required dependencies (pandas, numpy, jupyter) and launch the notebook in your browser.

## What's Inside

The notebook analyzes search algorithm results from CSV files, computing metrics like:
- Execution time
- Nodes expanded
- Solution length
- Heuristic values
- Success rates

Results are compared across different algorithms (BFS, IDFS, GBFS, A*, IDA*) and puzzle sizes (8-puzzle and 15-puzzle).