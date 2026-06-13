#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import opengate as gate
from opengate.tests import utility
import pathlib
import numpy as np
from box import Box
import matplotlib.pyplot as plt
from opengate.tests import utility
from test100_window_turbo_source_base import (
    build_geometry,
    compare_profiles,
    calculate_profile,
)

paths = utility.get_default_test_paths(__file__, output_folder="test100")


def change_source_parameters(source_back):
    keV = gate.g4_units.keV
    mm = gate.g4_units.mm
    source_back.particle = "gamma"
    source_back.energy.mono = 141 * keV
    source_back.position.type = "cylinder"
    source_back.position.translation = [0, -100, 0]
    source_back.position.radius = 100 * mm
    source_back.position.dz = 100 * mm


def initialize(duration=10):
    sim = gate.Simulation()
    sec = gate.g4_units.second
    sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option3"
    sim.physics_manager.global_production_cuts.all = 1 * gate.g4_units.mm
    # main options
    sim.g4_verbose = False
    sim.g4_verbose_level = 1
    sim.visu = False
    sim.visu_type = "qt"
    sim.number_of_threads = 32
    sim.progress_bar = True
    sim.run_timing_intervals = [[0, duration * sec]]
    sim.add_actor("SimulationStatisticsActor", "Stats")
    return sim


def run_window_turbo_source(skip: bool):
    Bq = gate.g4_units.Bq
    mm = gate.g4_units.mm
    sim = initialize(100)
    sim.g4_verbose = False
    sim.number_of_threads = 20
    sim.random_seed = 1
    sim.progress_bar = False
    # sim.visu = True
    # sim.visu_type = "qt"
    radius_down = 13.6
    build_geometry(
        sim, paths.output / f"window_turbo_skip_{skip}", pin_radius_down=radius_down
    )
    source_back = sim.add_source("WindowTurboSource", "source_back")

    source_back.activity = 1000000 * Bq
    source_back.direction.a1 = -radius_down * mm
    source_back.direction.a2 = radius_down * mm
    source_back.direction.b1 = -radius_down * mm
    source_back.direction.b2 = radius_down * mm
    source_back.direction.plane_distance = 86 * mm
    source_back.direction.plane_phi = np.pi / 2
    source_back.direction.skip_mode = skip
    if skip:
        source_back.direction.max_solid_angle = [0.09743848102142992]
    change_source_parameters(source_back)
    sim.run(start_new_process=True)
    stats = sim.get_actor("Stats")
    print(stats)
    print("-" * 80)


def run_generic_source(duration=100):
    Bq = gate.g4_units.Bq

    sim = initialize(duration=duration)
    build_geometry(sim, "generic")

    # physic list
    # print('Phys lists :', sim.get_available_physicLists())

    source_back = sim.add_source("GenericSource", "source_back")
    source_back.activity = 1000000 * Bq
    change_source_parameters(source_back)
    sim.number_of_threads = 20

    sim.run(start_new_process=True)

    stats = sim.get_actor("Stats")
    print(stats)
    print("-" * 80)


def analyze_next_time(filename_list):
    import uproot
    import os

    plt.figure()
    for i, filename in enumerate(filename_list):
        filename = str(filename)
        with uproot.open(filename) as f:
            # print(filename)
            trees = f.keys()
            tree = f["Singles"]
            data = tree.arrays(library="np")
            gt = data["GlobalTime"]
            # print("Total events:", len(gt))
            # print("max gt:", gt.max())
            gt = np.sort(gt)
            dt = np.diff(gt)
            # if i == 0:
            #     dt *= 10
            dt_hist, edge = np.histogram(dt, bins=400, range=(0, 1e7))
            print("dist_hist.sum():", dt_hist[:30].sum())
            plt.plot(
                edge[:-1],
                dt_hist,
                label=os.path.basename(filename).replace(".root", ""),
            )
    plt.legend()
    plt.xlabel("Time difference (s)")
    plt.ylabel("Counts")
    # plt.yscale("log")
    plt.title("Time difference between consecutive events")
    plt.tight_layout()
    plt.savefig("time_difference_histogram.png")
    plt.close()

    # print("Total energy deposit:", data["TotalEnergyDeposit"].sum())


if __name__ == "__main__":
    pathFile = pathlib.Path(__file__).parent.resolve()
    # run_generic_source()
    # run_window_turbo_source(skip=True)
    # run_window_turbo_source(skip=False)
    analyze_next_time(
        [
            "generic_singles.root",
            paths.output / "window_turbo_skip_True_singles.root",
            paths.output / "window_turbo_skip_False_singles.root",
        ]
    ),
    # run_window_turbo_source()

    # utility.test_ok(compare_result)
