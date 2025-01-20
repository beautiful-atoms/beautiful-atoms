import toml
import os

# Define paths
TOML_FILE = "batoms/blender_manifest.toml"
WHEELS_DIR = "batoms/wheels"

# Load existing TOML data
with open(TOML_FILE, "r") as f:
    data = toml.load(f)

# Get the list of wheel files
wheels = [
    f"./wheels/{file}" for file in os.listdir(WHEELS_DIR) if file.endswith(".whl")
]

# Update the TOML data
data["wheels"] = wheels

# Write the updated TOML file
with open(TOML_FILE, "w") as f:
    toml.dump(data, f)

print(f"Updated {TOML_FILE} with {len(wheels)} wheel files.")
