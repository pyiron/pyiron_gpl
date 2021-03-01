# coding: utf-8
# Copyright (c) ICAMS, Ruhr University Bochum, 2021
# Distributed under the terms of "GPLv3", see the LICENSE file.

import logging
import numpy as np
import os
import pandas as pd
import re
import ruamel.yaml as yaml

from shutil import copyfile

from pyiron_base import GenericJob, InputList, ImportAlarm, Settings
from pyiron_contrib.atomistic.atomistics.job.trainingcontainer import TrainingContainer


s = Settings()

try:
    from pyace import BBasisConfiguration, ACEBBasisSet

    import_alarm = ImportAlarm()

    # set loggers level for WARNING to avoid extra print-outs
    # because pyace do its own logging settings
    loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
    for logger in loggers:
        logger.setLevel(logging.WARNING)

except ImportError as e:
    import_alarm = ImportAlarm("Could not import `pyace` package. The package should be installed for proper "
                               "functionality")


class PaceMakerJob(GenericJob):
    """
    Thin wrapper class of the `pacemaker` - Atomic Cluster Expansion fitting code.

    Current functionality is limited to single-species fit with limited number of basis functions only.
    Please, contact developers if you would like to use fully functional `pacemaker` code.

    Usage example:

    job = fit_pr.create_job(job_type=PaceMakerJob, job_name="fit_job")

    # setup ACE potential form
    job.input["potential"]= {
            # spline mesh settings
            "deltaSplineBins": 0.001,

            # specie
            "element": "Cu",

            # embedding function settings
            "ndensity": 2,
            "fs_parameters": [1, 1, 1, 0.5],
            "npot": "FinnisSinclairShiftedScaled",

            # cutoff function
            "NameOfCutoffFunction": "cos",

            # potential specification
            ## radial basis functions type and parameters
            "radbase": "ChebExpCos",
            "radparameters": [5.25],
            "rcut": cutoff,
            "dcut": 0.01,

            ## max correlation order
            "rankmax": 3,
            ## specification of max n,l for each of the correlation order
            "nradmax": [5,2,1],
            "lmax": [0,2,1], ##NOTE: for order=1 lmax always is 0
    }

    # setup fitting: loss function, optimization settings
    job.input["fit"]= {
        'optimizer': 'BFGS',
        'maxiter': 150,
        'loss': {
            'kappa': 0.5,
            'L1_coeffs': 5e-7, # L1-regularization
            'L2_coeffs': 5e-7, # L2-regularization
            'w1_coeffs': 1,
            'w2_coeffs': 1,
            #radial smoothness regularization
            'w0_rad': 1e-4,
            'w1_rad': 1e-4,
            'w2_rad': 1e-4
        }
    }

    # setup global cutoff for atomic distances
    job.input["cutoff"] = cutoff

    # setup training data, could be:
    # - TrainingContainer
    # - pandas Dataframe
    # - filename of .pckl.gzip pandas Dataframe
    job.structure_data = data_job #

    """

    def __init__(self, project, job_name):
        super().__init__(project, job_name)
        self.__name__ = "PaceMakerJob"
        self.__version__ = "0.1"

        self.input = InputList(table_name="input")
        self.input['cutoff'] = 10.
        self.input['metadata'] = {}
        self.input['data'] = {}  # data_config
        self.input['potential'] = {}  # potential_config
        self.input['fit'] = {}  # fit_config
        self.input['backend'] = {'evaluator': 'pyace'}  # backend_config

        self.structure_data = None
        self._executable = None
        self._executable_activate()

    def _save_structure_dataframe_pckl_gzip(self, df):
        df.rename(columns={"number_of_atoms": "NUMBER_OF_ATOMS",
                           "energy": "energy_corrected",
                           "atoms": "ase_atoms"}, inplace=True)
        df["NUMBER_OF_ATOMS"] = df["NUMBER_OF_ATOMS"].astype(int)
        if "pbc" not in df.columns:
            df["pbc"] = df["ase_atoms"].map(lambda atoms: np.all(atoms.pbc))

        data_file_name = os.path.join(self.working_directory, "df_fit.pckl.gzip")
        logging.info(
            f"Saving training structures dataframe into {data_file_name} with pickle protocol = 4, compression = gzip")
        df.to_pickle(data_file_name, compression="gzip", protocol=4)
        return data_file_name

    def write_input(self):
        # prepare datafile
        if self.structure_data is None:
            raise ValueError(
                "`structure_data` is none, but should be pd.DataFrame, TrainingContainer or valid pickle.gzip filename")
        if isinstance(self.structure_data, pd.DataFrame):
            logging.info("structure_data is pandas.DataFrame")
            data_file_name = self._save_structure_dataframe_pckl_gzip(self.structure_data)
            self.input["data"] = {"filename": data_file_name}
        elif isinstance(self.structure_data, str):  # filename
            if os.path.isfile(self.structure_data):
                logging.info("structure_data is valid file path")
                self.input["data"] = {"filename": self.structure_data}
            else:
                raise ValueError(f"Provided structure_data filename ({self.structure_data}) doesn't exists")
        elif isinstance(self.structure_data, TrainingContainer):
            logging.info("structure_data is TrainingContainer")
            df = self.structure_data.to_pandas()
            data_file_name = self._save_structure_dataframe_pckl_gzip(df)
            self.input["data"] = {"filename": data_file_name}

        metadata_dict = self.input["metadata"]
        metadata_dict["pyiron_job_id"] = str(self.job_id)

        input_yaml_dict = {
            "cutoff": self.input["cutoff"],
            "metadata": metadata_dict,
            'potential': self.input['potential'],
            'data': self.input["data"],
            'fit': self.input["fit"],
            'backend': self.input["backend"],
        }

        if isinstance(self.input["potential"], str):
            pot_file_name = self.input["potential"]
            if os.path.isfile(pot_file_name):
                logging.info("Input potential is filename")
                pot_basename = os.path.basename(pot_file_name)
                copyfile(pot_file_name, os.path.join(self.working_directory, pot_basename))
                input_yaml_dict['potential'] = pot_basename
            else:
                raise ValueError(f"Provided potential filename ({self.input['potential']}) doesn't exists")

        with open(os.path.join(self.working_directory, "input.yaml"), "w") as f:
            yaml.dump(input_yaml_dict, f)

    def _analyse_log(self, logfile="log.txt"):
        log_filename = os.path.join(self.working_directory, logfile)

        with open(log_filename, "r") as f:
            loglines = f.readlines()

        losses = []
        ef_rmses = []

        for l in loglines:
            if "INFO" in l and "Iteration:" in l:
                loss = re.findall("Loss: ([\d.]*)", l)[0]
                losses.append(loss)

                ef_rmse_list = re.findall(
                    "RMSE Energy\(low\): ([0-9.]+) \(([0-9.]+)\) meV/at | Forces\(low\): ([0-9.]+) \(([0-9.]+)\) meV/A",
                    l)

                ef_rmses.append([ef_rmse_list[0][0], ef_rmse_list[0][1], ef_rmse_list[1][-2], ef_rmse_list[1][-1]])

        losses = np.array(losses).astype(float)

        ef_rmses = np.array(ef_rmses).astype(float)
        res_dict = {}
        res_dict["loss"] = losses
        res_dict["rmse_energy"] = ef_rmses[:, 0]
        res_dict["rmse_energy_low"] = ef_rmses[:, 1]
        res_dict["rmse_forces"] = ef_rmses[:, 2]
        res_dict["rmse_forces_low"] = ef_rmses[:, 3]
        return res_dict

    @import_alarm
    def collect_output(self):
        output_potential_filename = self.get_output_potential_filename()
        final_potential_filename = self.get_final_potential_filename()

        copyfile(output_potential_filename, final_potential_filename)
        with open(output_potential_filename, "r") as f:
            yaml_lines = f.readlines()
        final_potential_yaml_string = "".join(yaml_lines)

        # convert resulting potential to CTilde form and save
        bbasis = ACEBBasisSet(output_potential_filename)
        cbasis = bbasis.to_ACECTildeBasisSet()
        cbasis.save(self.get_final_potential_filename_ace())

        with open(self.get_final_potential_filename_ace(), "r") as f:
            ace_lines = f.readlines()
        final_potential_ace_string = "".join(ace_lines)

        elements_name = bbasis.elements_name

        with self.project_hdf5.open("output/potential") as h5out:
            h5out["yaml"] = final_potential_yaml_string
            h5out["ace"] = final_potential_ace_string
            h5out["elements_name"] = elements_name

        log_res_dict = self._analyse_log()

        with self.project_hdf5.open("output/log") as h5out:
            for key, arr in log_res_dict.items():
                h5out[key] = arr

    def get_lammps_potential(self):
        elements_name = self["output/potential/elements_name"]
        elem = " ".join(elements_name)
        pot_file_name = self.get_final_potential_filename_ace()
        pot_dict = {
            'Config': [["pair_style pace\n", f"pair_coeff  * * {pot_file_name} {elem}\n"]],
            'Filename': [""],
            'Model': ["ACE"],
            'Name': [self.job_name],
            'Species': [elements_name]
        }
        ace_potential = pd.DataFrame(pot_dict)

        return ace_potential

    def to_hdf(self, hdf=None, group_name=None):
        super().to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as h5in:
            self.input.to_hdf(h5in)

    def from_hdf(self, hdf=None, group_name=None):
        super().from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as h5in:
            self.input.from_hdf(h5in)

    def get_output_potential_filename(self):
        return os.path.join(self.working_directory, "output_potential.yaml")

    def get_final_potential_filename(self):
        return os.path.join(self.working_directory, self.job_name + ".yaml")

    def get_final_potential_filename_ace(self):
        return os.path.join(self.working_directory, self.job_name + ".ace")

    def get_current_potential_filename(self):
        return os.path.join(self.working_directory, "interim_potential_1.yaml")

    @import_alarm
    def get_current_potential(self):
        current_potential_filename = self.get_current_potential_filename()
        bbasis = BBasisConfiguration(current_potential_filename)
        return bbasis

    @import_alarm
    def get_final_potential(self):
        final_potential_filename = self.get_final_potential_filename()
        bbasis = BBasisConfiguration(final_potential_filename)
        return bbasis
