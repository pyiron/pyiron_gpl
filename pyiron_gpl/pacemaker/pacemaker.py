# coding: utf-8
# Copyright (c) ICAMS, Ruhr University Bochum, 2021
# Distributed under the terms of "GPLv3", see the LICENSE file.

import logging
import numpy as np
import os
import pandas as pd
import re
import ruamel.yaml as yaml
from pyiron_base import GenericJob, GenericParameters
from pyiron_base.settings.generic import Settings
from shutil import copyfile


s = Settings()

try:
    from pyace import BBasisConfiguration, ACEBBasisSet

    HAS_PYACE = True
except ImportError as e:
    print("Could not import `pyace` package. The package should be installed for proper functionality")
    HAS_PYACE = False

# set loggers
loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
for logger in loggers:
    logger.setLevel(logging.WARNING)

class PaceMakerJob(GenericJob):
    def __init__(self, project, job_name):
        super().__init__(project, job_name)
        self.__name__ = "PaceMakerJob"
        self.__version__ = "0.1"

        self.input = GenericParameters(table_name="input")
        self.input['cutoff'] = 10.
        self.input['metadata'] = {}
        self.input['data'] = {}  # data_config
        self.input['potential'] = {}  # potential_config
        self.input['fit'] = {}  # fit_config
        self.input['backend'] = {'evaluator': 'tensorpot'}  # backend_config

        self.structure_data = None

        # self.executable = "pacemaker input.yaml -l log.txt"
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
        logging.info("Saving training structures dataframe into {} with pickle protocol = 4, compression = gzip".format(
            data_file_name))
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
                raise ValueError("Provided structure_data filename ({}) doesn't exists".format(self.structure_data))
        elif hasattr(self.structure_data, "get_pandas"):  # duck-typing check for TrainingContainer
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
                # TODO: check if initial potential is provided (for continuation of fit)
            else:
                raise ValueError("Provided potential filename ({}) doesn't exists".format(self.input["potential"]))

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

    def collect_output(self):
        final_potential_filename = self.get_final_potential_filename()
        with open(final_potential_filename, "r") as f:
            yaml_lines = f.readlines()
        final_potential_yaml_string = "".join(yaml_lines)

        bbasis = ACEBBasisSet(final_potential_filename)
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
            'Config': [["pair_style pace\n", "pair_coeff  * * {} {}\n".format(pot_file_name, elem)]],
            'Filename': [""],
            'Model': ["ACE"],
            'Name': [self.job_name],
            'Species': [elements_name]
        }

        ace_potential = pd.DataFrame(pot_dict)

        return ace_potential

    def to_hdf(self, hdf=None, group_name=None):
        super().to_hdf(
            hdf=hdf,
            group_name=group_name
        )
        with self.project_hdf5.open("input") as h5in:
            self.input.to_hdf(h5in)

    def from_hdf(self, hdf=None, group_name=None):
        super().from_hdf(
            hdf=hdf,
            group_name=group_name
        )
        with self.project_hdf5.open("input") as h5in:
            self.input.from_hdf(h5in)

    def get_final_potential_filename(self):
        return os.path.join(self.working_directory, "output_potential.yaml")

    def get_final_potential_filename_ace(self):
        return os.path.join(self.working_directory, "output_potential.ace")

    def get_current_potential_filename(self):
        return os.path.join(self.working_directory, "interim_potential_1.yaml")

    def get_current_potential(self):
        if HAS_PYACE:
            current_potential_filename = self.get_current_potential_filename()
            bbasis = BBasisConfiguration(current_potential_filename)
            return bbasis
        else:
            raise RuntimeError("`pyace` package is not installed")

    def get_final_potential(self):
        if HAS_PYACE:
            final_potential_filename = self.get_final_potential_filename()
            bbasis = BBasisConfiguration(final_potential_filename)
            return bbasis
        else:
            raise RuntimeError("`pyace` pacakge is not installed")
