import copy
import os
import yaml

# note: we need to do something about the local directories.
# we don't create $HOME/temp and $HOME/Rangeland-Carbon/RCTM/tests/KFS (?)

# note: we need a better way to store gee_key_json file.

# note: if we can remove the local files, that's fine. otherwise, we need a
# dedicated folder to store them

# todo: turn this code into a command or a module. then invoke it with a single function call
def run_model(config_path):
    # todo: gotta get rid of these
    import sys

    sys.path.insert(0, "/home/dteber/carbon/Rangeland-Carbon")
    import os

    os.environ["RCTMPATH"] = "/home/dteber/carbon/Rangeland-Carbon/RCTM"
    os.environ["PYTHONPATH"] = "/home/dteber/carbon/Rangeland-Carbon"

    from RCTM.pipelines.RCTM_model_pipeline import RCTMPipeline

    config_filename = config_path
    Model_pipe = RCTMPipeline(config_filename=config_filename)
    Model_pipe.run_RCTM()


# question: can we combine these two into a single path
bucket_name = "rangelands"
workflow_base_dir = "RCTM_benchmark_sites/Ameriflux/KFS"
workflow_output_dir = os.path.join(workflow_base_dir, "RCTM_output")
home_dir = os.getenv("HOME")

test_cases = {
    "case_1": {
        "workflows_path": os.path.join(
            home_dir, "carbon", "Rangeland-Carbon/RCTM/tests/KFS"
        ),
        "path_to_temp_dir": os.path.join(
            home_dir, "temp"
        ),  # question: what's the purpose of this? can we get rid of this?
        "bucket_name": bucket_name,
        "gcloud_workflow_base_dir": workflow_base_dir,
        "service_account": "rctm-07dea6fa718fdd0921f124263@rangelands-explo-1571664594580.iam.gserviceaccount.com",
        "gee_key_json": os.path.join(home_dir, "res", "gee_key.json"),
        "gcloud_project": "rangelands-explo-1571664594580",
        "path_to_existing_ee_geom": "projects/rangelands-explo-1571664594580/assets/Ameriflux_RS/KFS/KFS",
        "sfm_imagery_download_date_range": ["2002-01-01", "2023-12-31"],
        "RCTM_spinup_date_range": ["2002-05-01", "2005-12-31"],
        "RCTM_transient_input_gen_date_range": ["2002-01-01", "2023-12-31"],
        "RCTM_transient_run_date_range": ["2002-01-01", "2002-01-25"],
        "spin_years": 0,
        "run_transient": True,
        "init_C_stocks_with_image": True,
        "C_stock_init_image": os.path.join(
            workflow_output_dir, "spinup", "RCTM_C_stocks_spin_output.tif"
        ),
        "transient_C_stock_hist": os.path.join(
            workflow_output_dir, "transient", "C_stock_test.nc"
        ),
        "transient_flux_hist": os.path.join(
            workflow_output_dir, "transient", "flux_hist_test.nc"
        ),
    },
    # define extra test cases here
    # ...
}

with open("benchmark_config.yaml") as stream:
    try:
        template_config = yaml.safe_load(stream)
    except yaml.YAMLError as exception:
        print(exception)

for name, case in test_cases.items():
    print(f"Preparing {name} case")
    test_case_config = copy.deepcopy(template_config)
    for key, val in case.items():
        test_case_config[key] = val

    test_case_config_path = f"{name}.yaml"
    with open(test_case_config_path, "w") as stream:
        yaml.dump(test_case_config, stream, default_flow_style=False)

    print(f"Running the model for {name}")
    run_model(test_case_config_path)
    os.remove(test_case_config_path)
