#!/bin/bash

echo "=== TWOM Voxel Plotting Example ==="
echo "Plotting individual subject effects at example coordinates..."
echo ""

# Plot effects at an example voxel
echo "Command: python ../twom.py nifti_files.txt --input_dir data/ --plotvoxel 45 54 18"
echo ""
python ../twom.py nifti_files.txt --input_dir data/ --plotvoxel 45 54 18

echo ""
echo "✓ Plot generated! This shows individual subject activation at voxel [45, 54, 18]"
echo "✓ If ../outputs/demo_twom_output.nii.gz exists, the overlap value is also displayed"
echo ""
echo "Try different coordinates:"
echo "  python ../twom.py nifti_files.txt --input_dir data/ --plotvoxel X Y Z"