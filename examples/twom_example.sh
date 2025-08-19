# examples/twom_example.sh
#!/bin/bash

echo "=== TWOM Example Analysis ==="
echo "Running TWOM on example subjects (language mapping data)..."
echo ""

echo "Key analysis modes:"
echo "  1. Relevant voxels (default): Group-consistent activations only"
echo "  2. All voxels: Includes single-subject effects"
echo ""

# Show both filename options
echo "Filename options:"
echo "  # Auto-generated: python ../twom.py nifti_files.txt --input_dir data/ --radius 2 --radius_units voxels"
echo "  # Custom name:    python ../twom.py nifti_files.txt --input_dir data/ --radius 2 --radius_units voxels --output_filename demo_twom_output.nii.gz"
echo ""

echo "Running: Relevant voxels analysis with custom filename..."
python ../twom.py nifti_files.txt --input_dir data/ --radius 2 --radius_units voxels --output_filename demo_twom_output.nii.gz

echo ""
echo "✓ Analysis complete! Results in ../outputs/"
echo "✓ TWOM map saved as: demo_relevant_voxels.nii.gz"
echo ""
echo "Next steps:"
echo "  - View results in your favourite neuroimaging viewer"
echo "  - Try: ./plot_example.sh to visualise individual voxel effects"
echo "  - Try: python ../twom.py nifti_files.txt --input_dir data/ --analyse_all_voxels --radius 2 --radius_units voxels --output_filename demo_all_voxels.nii.gz"
echo "  - The last option reruns twom.py but with explicit all_voxels option rather than default (relevant_voxels)"
