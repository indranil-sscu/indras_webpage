---
layout: page
title: "Batch Plotting in Python"
---

## Batch Plotting multiple datasets in Python

In many experimental workflows, you may encounter situations where you need to process multiple data files from different trials and overlay them on a single plot. This script streamlines what would otherwise be a tedious, manual process. A basic understanding of Bash or PowerShell is helpful. Since I primarily use a MacBook for data processing, this tutorial will focus on macOS-specific instructions. The same set of instructions should work seamlessly on Windows as well.

This code is based on an original script developed by Dr. Adam Luke Baskerville (https://adambaskerville.github.io/posts/PythonGUIPlotter/)[https://adambaskerville.github.io/posts/PythonGUIPlotter/] for a Python GUI plotting tool. While my version does not include the GUI functionality, it incorporates improved error handling and enhanced formatting for better usability.

```python
"""
## Python Batch Plotting Code
## author: indranil
## email: indranilroy@iisc.ac.in
## forked from: https://adambaskerville.github.io/posts/PythonGUIPlotter/
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import sys
import itertools

sns.set_style("ticks", {
    'axes.facecolor': 'beige',       # background color inside plot
    'grid.color': '0.8',             # light gray gridlines
    'grid.linestyle': '--',          # dashed gridlines
    'axes.edgecolor': 'black',       # color of border lines
    'axes.linewidth': 1.0,           # thickness of border lines
    'xtick.direction': 'in',         # ticks point inward
    'ytick.direction': 'in'
})

if len(sys.argv) < 3:
    print("\nNo data files provided as arguments. Usage: python script.py path data1.txt data2.txt ...")
    sys.exit(1)

data_dir = sys.argv[1]
if not os.path.isdir(data_dir):
    print(f"\nError: The directory '{data_dir}' does not exist.")
    sys.exit(1)
os.chdir(sys.argv[1])

# === Get axis labels ===
xAxisLabel = input("\nEnter x-axis label: ").strip()
yAxisLabel = input("\nEnter y-axis label: ").strip()

# === Setup ===
ListOfDataSets = []
LegendLabels = []
xcols = []
ycols = []
cols_to_use = []
plot_type = []

print("\nEnter the indices of the columns you want to plot. Eg. 1 for column 1, 2 for column 2, ...", end="")

# === Iterate over all files ===
for arg in sys.argv[2:]:
    if not os.path.exists(arg):
        print(f"Error: File '{arg}' not found.")
        sys.exit(1)

    ListOfDataSets.append(arg)

    try:
        xcol_idx = int(input(f"\nEnter x col idx for data set '{arg}': ")) - 1
        ycol_idx = int(input(f"Enter y col idx for data set '{arg}': ")) - 1
    except ValueError:
        print("Invalid column index. Please enter numeric values only.")
        sys.exit(1)

    if xcol_idx < 0 or ycol_idx < 0:
        print("Column indices must be >= 1.")
        sys.exit(1)

    xcols.append(xcol_idx)
    ycols.append(ycol_idx)
    cols_to_use.append([xcol_idx, ycol_idx])

    label = input(f"Enter legend label for data set '{arg}': ").strip()
    LegendLabels.append(label if label else f"{arg}")

    ptype = input("Plot Needed [Scatter/Line]: ").strip().capitalize()
    if ptype not in ["Scatter", "Line"]:
        print("Invalid plot type specified. Defaulting to 'Line'.")
        ptype = "Line"
    plot_type.append(ptype)

# === Plotting ===
WIDTH, HEIGHT, DPI = 500, 500, 100
fig, ax = plt.subplots(figsize=(WIDTH/DPI, HEIGHT/DPI))

N = len(ListOfDataSets)
palette = sns.color_palette("husl", N)  # or "tab20", "Set2", "Paired", etc.
colors = itertools.cycle(palette)  # infinite generator of colors

for i in range(len(ListOfDataSets)):
    path = ListOfDataSets[i]
    try:
        df = pd.read_csv(path, usecols=cols_to_use[i],
                         sep=r"\s+|\t+|\s+\t+|\t+\s+|,\s+|\s+,",
                         header=None, engine='python')
    except Exception as e:
        print(f"Failed to read file '{path}': {e}")
        sys.exit(1)

    color = next(colors)

    x_data = df[xcols[i]]
    y_data = df[ycols[i]]

    if plot_type[i] == "Scatter":
        ax.scatter(x_data, y_data, color=color, s=10, label=LegendLabels[i])
    else:
        ax.plot(x_data, y_data, color=color, label=LegendLabels[i], marker='o')

    ax.set_xlim(min(ax.get_xlim()[0], x_data.min()) - 0.1, max(ax.get_xlim()[1], x_data.max()) + 0.1)
    ax.set_ylim(min(ax.get_ylim()[0], y_data.min()) - 1, max(ax.get_ylim()[1], y_data.max()) + 1)

ax.set_xlabel(xAxisLabel, fontsize=14)
ax.set_ylabel(yAxisLabel, fontsize=14)
ax.legend()

plt.tight_layout()
plt.show()
```

To run the script, ensure that Seaborn, Pandas, and Matplotlib are installed in your Python environment. Open your terminal and use the following shell command to execute the script. It’s recommended to use the absolute path to the script and change the working directory (cd) to the location of your data files before running the code. There's no need to worry about delimiters—this script is designed to handle most commonly used separators automatically.
```bash
cd ~\path
python3 script_name.py . file_name1 file_name2 ...
```
If your data files are named sequentially or follow a consistent naming pattern, you can simplify the command using a wildcard. For example, the following will process all .txt files in the directory:
```bash
cd ~\path
python3 script_name.py . *.txt
```
Simply follow the on-screen prompts—and voilà! Your plot will be generated in no time.



