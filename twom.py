import argparse
import os
import numpy as np
import nibabel as nib
from scipy.stats import norm, t
from datetime import datetime
from time import time
import matplotlib.pyplot as plt
from tqdm import tqdm

def parse_arguments():
    parser = argparse.ArgumentParser(description="Threshold-Weighted Overlap Mapping (TWOM) script - Vectorized v1.")
    parser.add_argument('file_list', type=str, help='Path to the text file containing the list of NIfTI files.')
    parser.add_argument('--input_dir', type=str, default='', help='Input directory path to be prefixed to each item in the file list.')
    parser.add_argument('--output_filename', type=str, default=None, help='Output filename for the saved NIfTI map.')
    parser.add_argument('--radius', type=float, default=1, help='Radius in mm or voxels for the spherical region (default=1).')
    parser.add_argument('--radius_units', type=str, default='mm', choices=['mm', 'voxels'], help='Units of radius, either "mm" or "voxels" (default=mm).')
    parser.add_argument('--maskfile', type=str, default=None, help='Optional mask file (e.g. mask.nii.gz).')
    parser.add_argument('--analyse_all_voxels', action='store_true',
                       help='Analyse all voxels instead of relevant voxels (default: use relevant voxels like MATLAB).')
    parser.add_argument('--min_subjects', type=int, default=2,
                       help='Minimum subjects required for relevant voxels auto-mask (default=2).')
    parser.add_argument('--min_cluster_size', type=int, default=10,
                       help='Minimum cluster size for relevant voxels auto-mask (default=10 voxels).')
    parser.add_argument('--prob_min', type=float, default=0.01, help='Minimum probability threshold (default=0.01).')
    parser.add_argument('--prob_max', type=float, default=0.05, help='Maximum probability threshold (default=0.05).')
    parser.add_argument('--weight_type', type=str, default='linear', choices=['linear', 'quadratic', 'none'], help='Type of weighting function (default=linear).')
    parser.add_argument('--plotvoxel', type=int, nargs=3, help='Voxel coordinates to plot effects (x y z) - skips TWOM processing.')
    parser.add_argument('--suffix_to_drop', type=str, default='.nii.gz', help='Suffix to drop in file names when plotting effects (optional) (default=.nii.gz).')
    return parser.parse_args()

def read_file_list(file_list_path, input_dir):
    with open(file_list_path, 'r') as file:
        files = [os.path.join(input_dir, line.strip()) for line in file]
    return files

def get_data(file, verbose=False, return_img=False):
    if file is None:
        raise ValueError("File path is None")
    if not os.path.exists(file):
        raise FileNotFoundError(f"File does not exist: {file}")
    
    img = nib.load(file)
    data = img.get_fdata()
    if verbose:
        print(f'Loading: {file}')
        print(f'-- Dimensions: {img.shape}')
        print(f'-- Voxel sizes: {img.header.get_zooms()}')
    if return_img:
        return data, img
    return data

def check_dimensions(subject_files, mask_file_path=None):
    import sys
    first_file_dim = None
    for file_path in subject_files:
        subj_nii = nib.load(file_path)
        subj_dim = subj_nii.shape
        if first_file_dim is None:
            first_file_dim = subj_dim
        if subj_dim != first_file_dim:
            print(f"Error: Dimensions mismatch in subject files. {file_path} has dimensions {subj_dim}, expected {first_file_dim}.")
            sys.exit(1)
    if mask_file_path is not None:
        mask_nii = nib.load(mask_file_path)
        mask_dim = mask_nii.shape
        if mask_dim != first_file_dim:
            print(f"Error: Mask dimensions {mask_dim} do not match subject file dimensions {first_file_dim}.")
            sys.exit(1)
    print("All files are consistent in dimensions.")

def timeit(enabled=True):
    def decorator(func):
        if not enabled:
            return func
        def wrapper(*args, **kwargs):
            start_time = time()
            result = func(*args, **kwargs)
            end_time = time()
            print(f"Function '{func.__name__}' executed in {end_time - start_time:.4f} seconds")
            return result
        return wrapper
    return decorator

def extract_degrees_of_freedom(img_files):
    """Extract degrees of freedom from NIfTI headers"""
    df_values = []
    
    for file_path in img_files:
        try:
            img = nib.load(file_path)
            descrip = img.header.get('descrip', b'').decode('utf-8')
            
            import re
            df_match = re.search(r'\[(\d+)\]', descrip)
            if df_match:
                df_values.append(int(df_match.group(1)))
            else:
                df_match = re.search(r'df=(\d+)', descrip, re.IGNORECASE)
                if df_match:
                    df_values.append(int(df_match.group(1)))
        except:
            continue
    
    if df_values:
        unique_df = list(set(df_values))
        if len(unique_df) == 1:
            return unique_df[0]
        else:
            print(f"Variable d.f. values found: {unique_df}")
            print("Taking average d.f. value")
            return int(np.mean(df_values))
    else:
        print("Warning: Could not extract degrees of freedom from headers")
        print("Using large df approximation (df=1000) - equivalent to normal distribution")
        return 1000

def calculate_thresholds(prob_min, prob_max, num_thresholds=100, tails=2, df=None):
    """Calculate t-thresholds using t-distribution or normal distribution"""
    p_values = np.logspace(np.log10(prob_max), np.log10(prob_min), num_thresholds) / tails
    
    if df is not None and df < 1000:
        t_thresholds = t.ppf(1 - p_values, df)
        print(f"Using t-distribution with df={df}")
    else:
        t_thresholds = norm.ppf(1 - p_values)
        if df is not None:
            print(f"Using normal distribution (df={df} is large)")
        else:
            print("Using normal distribution (no df information)")
    
    if prob_min == 0.5:
        t_thresholds[0] = 0
    return t_thresholds

def linear_weighting_function(midpoints, t_min, t_max):
    return 2 * (midpoints - t_min) / (t_max - t_min)

def quadratic_weighting_function(midpoints, t_min, t_max):
    return 3 * ((midpoints - t_min) ** 2) / ((t_max - t_min) ** 2)

def calculate_thresholds_and_weights(prob_min, prob_max, num_thresholds=100, weight_type='linear', df=None):
    t_thresholds = calculate_thresholds(prob_min, prob_max, num_thresholds, df=df)
    t_min, t_max = np.min(t_thresholds), np.max(t_thresholds)
    if weight_type == 'linear':
        weights = linear_weighting_function(t_thresholds, t_min, t_max)
    elif weight_type == 'quadratic':
        weights = quadratic_weighting_function(t_thresholds, t_min, t_max)
    else:
        weights = np.ones_like(t_thresholds)
    return t_thresholds, weights

def convert_radius_to_voxels(radius_mm, voxel_sizes):
    return np.array(radius_mm) / np.array(voxel_sizes)

def get_sphere_offsets(radius):
    """Generate 3D offsets for spherical neighborhood"""
    r = int(radius)
    offsets = []
    for x in range(-r, r+1):
        for y in range(-r, r+1):
            for z in range(-r, r+1):
                if x**2 + y**2 + z**2 <= r**2:
                    offsets.append([x, y, z])
    return np.array(offsets)

def create_relevant_voxels_mask(all_data, threshold, min_subjects=2, min_cluster_size=10):
    """
    Create relevant voxels mask
    Excludes voxels activated in <min_subjects and removes tiny clusters
    """
    n_subjects, ny, nx, nz = all_data.shape
    
    activation_count = np.sum(all_data > threshold, axis=0)
    mask = activation_count >= min_subjects
    
    print(f"Initial relevant voxels: {np.sum(mask):,} voxels (≥{min_subjects} subjects above threshold {threshold:.3f})")
    
    try:
        from scipy import ndimage
        
        structure = ndimage.generate_binary_structure(3, 2)
        labeled_mask, num_clusters = ndimage.label(mask, structure=structure)
        
        print(f"Found {num_clusters} clusters before size filtering")
        
        for cluster_id in range(1, num_clusters + 1):
            cluster_voxels = np.sum(labeled_mask == cluster_id)
            if cluster_voxels < min_cluster_size:
                mask[labeled_mask == cluster_id] = False
        
        final_voxels = np.sum(mask)
        removed_voxels = np.sum(activation_count >= min_subjects) - final_voxels
        print(f"Removed {removed_voxels:,} voxels in clusters <{min_cluster_size} voxels")
        print(f"Final relevant voxels: {final_voxels:,} voxels")
        
    except ImportError:
        print("Warning: scipy.ndimage not available, skipping cluster size filtering")
        print("Install scipy for full relevant voxels functionality")
    
    return mask

def simple_sphere_extraction(all_data, center, sphere_offsets):
    """Extract data from spherical neighborhood around center voxel"""
    n_subjects = all_data.shape[0]
    ny, nx, nz = all_data.shape[1:]
    
    sphere_values = []
    
    for offset in sphere_offsets:
        neighbor = center + offset
        
        if (neighbor[0] >= 0 and neighbor[0] < ny and
            neighbor[1] >= 0 and neighbor[1] < nx and
            neighbor[2] >= 0 and neighbor[2] < nz):
            
            neighbor_data = all_data[:, neighbor[0], neighbor[1], neighbor[2]]
            sphere_values.append(neighbor_data)
    
    if len(sphere_values) > 0:
        return np.array(sphere_values).T
    else:
        return np.array([]).reshape(n_subjects, 0)

def calculate_voxel_twom_vectorized(sphere_data, t_thresholds, weights, t_min, t_max):
    """
    VECTORIZED: Calculate TWOM value for a single voxel using vectorized operations
    This eliminates the inner threshold loop for significant speedup
    """
    if sphere_data.shape[1] == 0:
        return 0.0
    
    # Vectorized threshold exceedance calculation
    # Shape: (n_subjects, n_neighbors, n_thresholds)
    exceedances = sphere_data[:, :, np.newaxis] > t_thresholds[np.newaxis, np.newaxis, :]
    
    # Calculate proportions across subjects for each neighbor and threshold
    # Shape: (n_neighbors, n_thresholds)
    proportions = exceedances.mean(axis=0)
    
    # Take maximum across sphere (like MATLAB) for each threshold
    # Shape: (n_thresholds,)
    max_proportions = proportions.max(axis=0)
    
    # Apply weights and integrate
    weighted_proportions = max_proportions * weights
    
    # Calculate AUC (area under curve)
    auc_value = np.trapz(weighted_proportions, t_thresholds) / (t_max - t_min)
    
    return auc_value

@timeit(enabled=True)
def process_volume(subject_files, mask_file, prob_min, prob_max, weight_type, radius, radius_units='mm', 
                  analyse_all_voxels=False, min_subjects=2, min_cluster_size=10):
    """Process volume using threshold-weighted overlap mapping - VECTORIZED v1"""
    example_img_data, example_img = get_data(subject_files[0], return_img=True)
    voxel_sizes = example_img.header.get_zooms()[:3]
    ny, nx, nz = example_img.shape
    n_subjects = len(subject_files)

    if radius_units == 'mm':
        radius_voxels = np.round(convert_radius_to_voxels(radius, voxel_sizes)).astype(int).max()
    elif radius_units == 'voxels':
        radius_voxels = int(radius)
    else:
        raise ValueError("radius_units must be either 'mm' or 'voxels'.")

    print(f"Using sphere radius of {radius_voxels} voxels")

    cst = np.zeros((ny, nx, nz))

    df = extract_degrees_of_freedom(subject_files)

    print("Loading subject data...")
    all_data = np.array([get_data(file) for file in subject_files])
    print(f"Loaded data shape: {all_data.shape}")

    if mask_file is not None:
        mask_img_data, mask_img = get_data(mask_file, return_img=True)
        mask = mask_img_data.astype(bool)
        print("Using user-provided mask")
    elif analyse_all_voxels:
        mask = np.ones((ny, nx, nz), dtype=bool)
        print("Using all voxels")
    else:
        print("Using MATLAB default: auto-generating relevant voxels mask")
        print("Focus: Group-level consistency (excludes subject-specific activations)")
        print("To analyse all voxels instead, use --analyse_all_voxels flag")
        
        t_thresholds_temp, _ = calculate_thresholds_and_weights(prob_min, prob_max, df=df)
        t_min = np.min(t_thresholds_temp)
        
        mask = create_relevant_voxels_mask(
            all_data, t_min, 
            min_subjects=min_subjects,
            min_cluster_size=min_cluster_size
        )

    t_thresholds, weights = calculate_thresholds_and_weights(
        prob_min, prob_max, num_thresholds=100, weight_type=weight_type, df=df
    )
    t_min, t_max = np.min(t_thresholds), np.max(t_thresholds)

    sphere_offsets = get_sphere_offsets(radius_voxels)
    print(f"Sphere contains {len(sphere_offsets)} voxels")

    mask_coords = np.array(np.where(mask)).T
    n_masked_voxels = len(mask_coords)
    print(f"Processing {n_masked_voxels} masked voxels")

    with tqdm(total=n_masked_voxels, desc="Processing Voxels", unit="voxel") as pbar:
        for i, coord in enumerate(mask_coords):
            
            # Extract sphere data using proven method (unchanged)
            sphere_data = simple_sphere_extraction(all_data, coord, sphere_offsets)
            
            # VECTORIZED: Calculate TWOM value using vectorized threshold processing
            auc_value = calculate_voxel_twom_vectorized(sphere_data, t_thresholds, weights, t_min, t_max)
            
            # Store result
            cst[coord[0], coord[1], coord[2]] = auc_value
            pbar.update(1)

    cst = np.clip(cst, 0, 1)
    print(f"Result range: [{cst.min():.6f}, {cst.max():.6f}]")
    
    return cst

def save_results(cst, affine, header, n_subj, input_list_name, radius, radius_units, weight_type, savefile=None, directory='outputs', method_tag='relevantvoxels'):
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    list_basename = os.path.splitext(os.path.basename(input_list_name))[0]
    
    if not savefile:
        savefile = f'TWOM_{list_basename}_N{n_subj}_r{radius}{radius_units}_{weight_type}_{method_tag}_{datetime.now().strftime("%Y%m%d_%H%M%S")}.nii.gz'
    
    save_path = os.path.join(directory, savefile)
    
    comment = f'TWOM {{N={n_subj} from {list_basename}}} radius={radius}{radius_units} weight={weight_type} method={method_tag} created={datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'
    
    new_header = header.copy()
    if len(comment) > 79:
        comment = comment[:76] + '...'
    new_header['descrip'] = comment.encode('utf-8')
    
    nib.save(nib.Nifti1Image(cst, affine, new_header), save_path)
    print(f"TWOM map has been saved to {save_path}")
    print(f"Header description: {comment}")
    
    return save_path

def plot_effects(subj_effects, all_subj, voxel_coordinates, suffix_to_drop="_zmap3_MNI2.nii.gz"):
    all_subj = [os.path.basename(subj).replace(suffix_to_drop, '') for subj in all_subj]
    plt.figure(figsize=(10,6))
    plt.bar(all_subj, subj_effects, color='skyblue')
    plt.title(f'Summary: Effect by Subject, Voxel: {voxel_coordinates}')
    plt.xlabel('Subjects', fontsize=14)
    plt.ylabel('Effect Sizes', fontsize=14)
    plt.grid(axis='y')
    plt.xticks(rotation=45, ha='right', fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.show()

def get_voxel_effects(all_data, voxel):
    return all_data[:, voxel[0], voxel[1], voxel[2]]

def view_voxel_effects(voxel, all_data, all_subj, suffix_to_drop):
    subj_effects = get_voxel_effects(all_data, voxel)
    plot_effects(subj_effects, all_subj, voxel, suffix_to_drop)

def find_existing_twom_file(file_list_name):
    """Try to find existing TWOM output file"""
    list_basename = os.path.splitext(os.path.basename(file_list_name))[0]
    output_dir = 'outputs'
    
    if os.path.exists(output_dir):
        for filename in os.listdir(output_dir):
            if filename.startswith(f'TWOM_{list_basename}') and filename.endswith('.nii.gz'):
                return os.path.join(output_dir, filename)
    return None

def main():
    args = parse_arguments()

    files = read_file_list(args.file_list, args.input_dir)
    n_subj = len(files)

    check_dimensions(files, args.maskfile)

    all_data = np.array([get_data(file) for file in files])
    print(f"Data value range: min={all_data.min()}, max={all_data.max()}")
    
    if args.plotvoxel:
        voxel = tuple(args.plotvoxel)
        print(f"Plotting effects for voxel {voxel} (no TWOM processing required)")
        
        twom_file = find_existing_twom_file(args.file_list)
        twom_value = None
        if twom_file:
            print(f"Found existing TWOM file: {twom_file}")
            try:
                twom_img = nib.load(twom_file)
                twom_data = twom_img.get_fdata()
                twom_value = twom_data[voxel[0], voxel[1], voxel[2]]
                print(f"TWOM value at {voxel}: {twom_value:.6f}")
            except Exception as e:
                print(f"Warning: Could not read TWOM value: {e}")
        
        view_voxel_effects(voxel, all_data, files, suffix_to_drop=args.suffix_to_drop)
    else:
        print("Extracting degrees of freedom from input files...")
        cst = process_volume(files, args.maskfile, args.prob_min, args.prob_max, args.weight_type, args.radius, args.radius_units,
                            analyse_all_voxels=args.analyse_all_voxels, 
                            min_subjects=args.min_subjects, 
                            min_cluster_size=args.min_cluster_size)

        example_img_data, example_img = get_data(files[0], return_img=True)
        # Determine processing method for filename
        processing_method = 'allvoxels' if args.analyse_all_voxels else 'relevantvoxels'
        
        save_results(cst, example_img.affine, example_img.header, n_subj, 
                    args.file_list, args.radius, args.radius_units, args.weight_type,
                    savefile=args.output_filename, directory='outputs', method_tag=processing_method)

if __name__ == "__main__":
    main()
