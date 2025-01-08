import scipy.io
import pandas as pd
import numpy as np

# Step 1: Load the .mat file
name = 'geo_recon_corr_vertex_struct'

mat_data = scipy.io.loadmat(f'{name}.mat')

# Remove metadata keys
data_keys = [key for key in mat_data.keys() if not key.startswith('__')]

# Step 2: Process the Data
for key in data_keys:
    data = mat_data[key]
    print(f"Processing key: {key}, Type: {type(data)}, Shape: {getattr(data, 'shape', 'N/A')}")

    if isinstance(data, np.ndarray) and data.shape == (1, 1):
        # Unpack the single-element array
        data = data[0, 0]

    if isinstance(data, np.ndarray):
        # Assuming the first row contains numbers 0 to 46, drop it
        data = data[1:]  # Keep only the rows after the first row

    # Extract the last non-zero elements in each row of the nested arrays
    extracted_data = []  # List to store the last non-zero elements
    for row in data:
        if isinstance(row, np.ndarray):
            flattened_row = row.flatten()  # Flatten the row to a 1D array
            non_zero_elements = flattened_row[flattened_row != 0]  # Filter non-zero elements
            extracted_data.append(non_zero_elements)  # Append 0 if no non-zero elements

    # Convert the list of extracted data into a DataFrame with each row on its own line
    final_data = pd.DataFrame(extracted_data)  # Each extracted value becomes a row

    # Save to CSV
    csv_filename = f'{name}.csv'
    final_data.to_csv(csv_filename, index=False, header=False)
    print(f"Saved processed data to {csv_filename}")