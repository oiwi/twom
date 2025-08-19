#!/bin/bash

# Exit on any error
set -e

echo "=== Creating twom environment ==="

# Remove existing environment if it exists
if conda env list | grep -q "twom"; then
    echo "Removing existing 'twom' environment..."
    conda env remove --name twom -y
fi

# Create environment
echo "Creating new environment..."
conda create --name twom python=3.8 -y

# Deactivate any current environment
conda deactivate 2>/dev/null || true

# Activate the new environment
echo "Activating environment..."
eval "$(conda shell.bash hook)"  # Ensure conda works in scripts
conda activate twom

# Check if requirements.txt exists
if [ ! -f "requirements.txt" ]; then
    echo "ERROR: requirements.txt not found!"
    exit 1
fi

echo "Installing packages from requirements.txt..."
echo "Contents of requirements.txt:"
cat requirements.txt

# Try bulk install first, then individual if it fails
echo "Attempting bulk install..."
if pip install -r requirements.txt -v; then
    echo "✓ Successfully installed packages from requirements.txt"
else
    echo "❌ Bulk install failed, trying individual installations..."
    
    # Fallback: install each package individually
    while IFS= read -r package; do
        # Skip empty lines and comments
        if [[ ! -z "$package" && ! "$package" =~ ^# ]]; then
            echo "Installing $package..."
            if pip install "$package"; then
                echo "✓ $package installed successfully"
            else
                echo "❌ Failed to install $package"
                exit 1
            fi
        fi
    done < requirements.txt
fi

# Verify installation
echo "=== Verifying installation ==="
echo "Installed packages:"
pip list

echo "Testing imports..."
python -c "
import sys
packages = ['numpy', 'nibabel', 'scipy', 'matplotlib', 'tqdm']
failed = []

for pkg in packages:
    try:
        __import__(pkg)
        print(f'✓ {pkg} imported successfully')
    except ImportError as e:
        print(f'❌ {pkg} failed: {e}')
        failed.append(pkg)

if failed:
    print(f'Failed packages: {failed}')
    sys.exit(1)
else:
    print('🎉 All packages working!')
"

echo "=== Environment ready! ==="
echo "To use: conda activate twom"