import gcsfs
import numpy as np
import xarray as xr
from pathlib import Path

# note: the process is pretty manual at the moment. I still haven't figured out how to
# store the base data. in the current implementation, we run the model and store the results
# in a bucket. in the future when we make a change to the model, we can compare our recent
# results against to these results. however, this approach isn't robust in many ways and it's
# very tedious. the next step would be versioning the current model and combining these two
# scripts into a single one. that way, we can generate the same data everytime we run the tests
# and extend the tests as we go. optional: we can check if the test case exists in the bucket.
# if so, we don't rerun the model again to produce test data. that way we can save time as well.

# todo: update dependencies to add h5netcdf

GIVEN_PATH_DOESNT_EXIST = "The given path doesn't exist: {path}"


def compare_netcdf_files(ds1: xr.Dataset, ds2: xr.Dataset) -> None:
    dim_diff = set(ds1.dims) ^ set(ds2.dims)
    if dim_diff:
        print(f"Dimension mismatch: {dim_diff}")
        return

    for dim in ds1.dims:
        if ds1.dims[dim] != ds2.dims[dim]:
            print(
                f"Dimension size mismatch for {dim}: {ds1.dims[dim]} != {ds2.dims[dim]}"
            )

    var_diff = set(ds1.data_vars) ^ set(ds2.data_vars)
    if var_diff:
        print(f"Variable mismatch: {var_diff}")

    common_vars = set(ds1.data_vars) & set(ds2.data_vars)
    for var in common_vars:
        arr1 = ds1[var].values
        arr2 = ds2[var].values

        if arr1.shape != arr2.shape:
            print(f"Shape mismatch for variable {var}: {arr1.shape} != {arr2.shape}")
            continue

        if not np.allclose(arr1, arr2, equal_nan=True):
            print(f"Value mismatch for variable {var}")

    attr_diff = set(ds1.attrs.keys()) ^ set(ds2.attrs.keys())
    if attr_diff:
        print(f"Global attribute mismatch: {attr_diff}")

    for attr in ds1.attrs:
        if ds1.attrs.get(attr) != ds2.attrs.get(attr):
            print(
                f"Global attribute value mismatch for {attr}: {ds1.attrs.get(attr)} != {ds2.attrs.get(attr)}"
            )

    for var in common_vars:
        attr_diff = set(ds1[var].attrs.keys()) ^ set(ds2[var].attrs.keys())
        if attr_diff:
            print(f"Attribute mismatch for variable {var}: {attr_diff}")

        for attr in ds1[var].attrs:
            if ds1[var].attrs.get(attr) != ds2[var].attrs.get(attr):
                print(
                    f"Attribute value mismatch for {var} -> {attr}: {ds1[var].attrs.get(attr)} != {ds2[var].attrs.get(attr)}"
                )

    print("Comparison complete. The files are identical!")


fs = gcsfs.GCSFileSystem(project="rangelands-explo-1571664594580", token=None)

path_to_expected_data = Path("rangelands/RCTM_comparison_test_data/KFS")
path_to_actual_data = Path("rangelands/RCTM_benchmark_sites/Ameriflux/KFS")

assert fs.exists(path_to_expected_data), GIVEN_PATH_DOESNT_EXIST.format(
    path_to_expected_data
)
assert fs.exists(path_to_actual_data), GIVEN_PATH_DOESNT_EXIST.format(
    path_to_actual_data
)

output_folder = "RCTM_output"

expected_output_path = path_to_expected_data / output_folder
actual_output_path = path_to_actual_data / output_folder

assert fs.exists(expected_output_path), GIVEN_PATH_DOESNT_EXIST.format(
    expected_output_path
)
assert fs.exists(actual_output_path), GIVEN_PATH_DOESNT_EXIST.format(actual_output_path)
assert len(fs.ls(expected_output_path)) == len(fs.ls(actual_output_path))

expected_transient = expected_output_path / "transient"
actual_transient = actual_output_path / "transient"

FILES_TO_COMPARE = [
    "C_stock_hist.nc",
    "flux_hist.nc",
]

for file in FILES_TO_COMPARE:
    temp_path = actual_transient / file
    assert fs.exists(temp_path), GIVEN_PATH_DOESNT_EXIST.format(temp_path)

    temp_path = expected_transient / file
    assert fs.exists(temp_path), GIVEN_PATH_DOESNT_EXIST.format(temp_path)


for file in FILES_TO_COMPARE:
    print(f"Reading {file} to compare...")
    actual_ds = xr.open_dataset(
        fs.open(actual_transient / file, block_size=8 * 1024 * 1024), engine="h5netcdf"
    )
    expected_ds = xr.open_dataset(
        fs.open(expected_transient / file, block_size=8 * 1024 * 1024),
        engine="h5netcdf",
    )

    print(f"Comparing {file}...")
    compare_netcdf_files(actual_ds, expected_ds)
