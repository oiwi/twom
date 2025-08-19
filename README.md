# TWOM - Threshold-Weighted Overlap Mapping

Vectorised, command-line Python implementation of TWOM for characterising individual variation in neuroimaging. Quantifies spatial overlap across multiple statistical maps. 

## Quick Setup

```bash
git clone https://github.com/yourusername/twom.git
cd twom
chmod +x create_env.sh
./create_env.sh
conda activate twom
```

## Usage

### 1. Get Help
```bash
python twom.py --help
```

### 2. Basic Analysis

Create file list (one NIfTI file per line):
```
subject1_tstat.nii.gz
subject2_tstat.nii.gz
subject3_tstat.nii.gz
```

Run standard analysis:
```bash
# Group-consistent voxels only (default)
python twom.py file_list.txt

# With custom input directory
python twom.py file_list.txt --input_dir /path/to/data/

# Custom sphere radius
python twom.py file_list.txt --radius 2 --radius_units voxels
```

### 3. Comprehensive Analysis
```bash
# Analyse all voxels (not just group-consistent ones)
python twom.py file_list.txt --analyse_all_voxels --radius 2 --radius_units voxels
```

### 4. Visualise Individual Voxel Effects
```bash
# Plot effects at specific voxel coordinates (x y z)
python twom.py file_list.txt --plotvoxel 50 50 50

# With custom parameters
python twom.py file_list.txt --plotvoxel 50 50 50 --radius 2 --radius_units voxels
```

This creates a bar plot showing effect sizes across subjects at the specified voxel. If an existing TWOM map is found, it will also display the TWOM value at that location.

## Key Options

- `--analyse_all_voxels` - Analyse all brain voxels (default: relevant voxels only)
- `--radius 2` - Sphere radius in mm or voxels (default: 1)
- `--radius_units voxels` - Units for radius: 'mm' or 'voxels' (default: mm)
- `--input_dir /path/` - Directory prefix for input files
- `--prob_min 0.001` - Minimum p-value threshold (default: 0.01)
- `--prob_max 0.05` - Maximum p-value threshold (default: 0.05)
- `--weight_type quadratic` - Weighting function: 'linear', 'quadratic', 'none' (default: linear)
- `--maskfile mask.nii.gz` - Optional brain mask
- `--plotvoxel x y z` - Plot effects at voxel coordinates (skips TWOM processing)
- `--min_subjects 3` - Minimum subjects for relevant voxels (default: 2)
- `--min_cluster_size 20` - Minimum cluster size for relevant voxels (default: 10)

## Analysis Modes

**Relevant Voxels (Default):**
- Focuses on group-consistent activations
- Requires activation in ≥2 subjects with clusters ≥10 voxels
- Faster processing, emphasises reproducible effects

**All Voxels:**
- Comprehensive analysis including single-subject effects
- Use `--analyse_all_voxels` flag
- Slower but more complete coverage

## Example Workflows

### Standard Group Analysis
```bash
python twom.py group_data.txt --input_dir /data/group1/ --radius 2 --radius_units voxels
```

### Individual Differences Analysis
```bash
python twom.py individual_data.txt --analyse_all_voxels --prob_min 0.001
```

### Explore Specific Region
```bash
# First run analysis
python twom.py data.txt --radius 2 --radius_units voxels

# Then plot specific voxel
python twom.py data.txt --plotvoxel 45 30 25
```

## Output

Results saved as `TWOM_*.nii.gz` in `outputs/` directory. Values range 0-1 (higher = more overlap).

Plotting generates matplotlib figures showing individual subject effects at specified voxels.

## Requirements

Python 3.8+, installed automatically with `create_env.sh`

## Citation

If you use this tool, please cite both:

**Original TWOM method paper (includes Matlab code):**
Seghier ML, Price CJ. Visualising inter-subject variability in fMRI using threshold-weighted overlap maps. *Sci Rep*. 2016;6(1):20170. doi:10.1038/srep20170

**This Python implementation:**
Voets, NL, Parker Jones, O, Seghier, ML, Stacey, R, Plaha, P. (2025) Separating the forest from the palm trees: Individual variation in a presurgical language mapping task. *Manuscript in preparation*
