#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import opengate as gate
import uproot
import SimpleITK as sitk

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


def change_source_parameters(source_back, source_1, source_2, source_3):
    keV = gate.g4_units.keV
    mm = gate.g4_units.mm
    source_back.particle = "gamma"
    source_back.energy.mono = 141 * keV
    source_back.position.type = "cylinder"
    source_back.position.translation = [0, -100, 0]
    source_back.position.radius = 100 * mm
    source_back.position.dz = 100 * mm

    source_1.particle = "gamma"
    source_1.energy.mono = 141 * keV
    source_1.position.type = "cylinder"
    source_1.position.translation = [-50, -100, 0]
    source_1.position.radius = 10 * mm
    source_1.position.dz = 100 * mm

    source_2.particle = "gamma"
    source_2.energy.mono = 141 * keV
    source_2.position.type = "cylinder"
    source_2.position.translation = [0, -100, 0]
    source_2.position.radius = 15 * mm
    source_2.position.dz = 100 * mm

    source_3.particle = "gamma"
    source_3.energy.mono = 141 * keV
    source_3.position.type = "cylinder"
    source_3.position.translation = [50, -100, 0]
    source_3.position.radius = 20 * mm
    source_3.position.dz = 100 * mm


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


def add_phantom(sim):
    mm = gate.g4_units.mm
    phantom = sim.add_volume("Tubs", "phantom")
    phantom.rmin = 0
    phantom.rmax = 100 * mm
    phantom.dz = 100 * mm
    phantom.material = "G4_WATER"
    phantom.color = [0, 0, 1, 1]
    phantom.translation = [0, -100 * mm, 0]
    phantom.mother = "world"


def run_window_turbo_source(ff=False, thread_count=4):
    print("running with ff =", ff, "thread_count =", thread_count)
    Bq = gate.g4_units.Bq
    mm = gate.g4_units.mm
    activity = 1e6
    # sim = initialize(781.25)
    sim = initialize(400)
    if ff:
        sim.g4_commands_before_init.append("/process/em/UseGeneralProcess false")
    sim.g4_verbose = False
    sim.number_of_threads = thread_count
    sim.random_seed = 1
    sim.progress_bar = False
    radius_down = 13.6
    build_geometry(
        sim,
        paths.output / f"window_turbo_ff_{ff}",
        pin_radius_down=radius_down,
        count_scatter=not ff,
    )
    add_phantom(sim)
    sim.physics_manager.set_production_cut("phantom", "all", 1 * mm)
    source_back = sim.add_source("WindowTurboSource", "source_back")
    source_1 = sim.add_source("WindowTurboSource", "source_1")
    source_2 = sim.add_source("WindowTurboSource", "source_2")
    source_3 = sim.add_source("WindowTurboSource", "source_3")
    source_1.activity = activity * Bq
    source_2.activity = activity * Bq
    source_3.activity = activity * Bq
    source_back.activity = activity * Bq
    source_back.direction.a1 = -radius_down * mm
    source_back.direction.a2 = radius_down * mm
    source_back.direction.b1 = -radius_down * mm
    source_back.direction.b2 = radius_down * mm
    source_back.direction.plane_distance = 86 * mm
    source_back.direction.plane_phi = np.pi / 2
    source_2.direction = source_back.direction.copy()
    source_3.direction = source_back.direction.copy()
    source_1.direction = source_back.direction.copy()
    change_source_parameters(source_back, source_1, source_2, source_3)
    if ff:
        ff_actor = sim.add_actor("GammaFreeFlightActor", "ff")
        ff_actor.attached_to = "world"
        ff_actor.exclude_volumes = ["head"]
    sim.run(start_new_process=False)
    stats = sim.get_actor("Stats")
    print(stats)
    print("-" * 80)


def run_generic_source(job_count):
    for ijob in range(job_count):
        activity = 1000000
        Bq = gate.g4_units.Bq

        sim = initialize(15.625)
        sim.random_seed = ijob + 1
        build_geometry(sim, f"generic_{ijob}", count_scatter=True)

        mm = gate.g4_units.mm
        add_phantom(sim)
        sim.visu = False
        sim.visu_type = "qt"
        sim.number_of_threads = 32

        # physic list
        # print('Phys lists :', sim.get_available_physicLists())

        source_back = sim.add_source("GenericSource", "source_back")
        source_1 = sim.add_source("GenericSource", "source_1")
        source_2 = sim.add_source("GenericSource", "source_2")
        source_3 = sim.add_source("GenericSource", "source_3")
        source_1.activity = activity * Bq
        source_2.activity = activity * Bq
        source_3.activity = activity * Bq
        source_back.activity = activity * Bq
        change_source_parameters(source_back, source_1, source_2, source_3)
        sim.run(start_new_process=True)

        stats = sim.get_actor("Stats")
        print(stats)
        print("-" * 80)


def analyze_root(file):
    with uproot.open(file) as f:
        if "Singles" not in f:
            return None
        tree = f["Singles"]
        tree_np = tree.arrays(library="np")
        print(
            "Weight stats:",
            tree_np["Weight"].max(),
            tree_np["Weight"].min(),
            tree_np["Weight"].mean(),
        )
        print(tree_np.keys())


def generate_prj_one(path):
    with uproot.open(path) as file:
        if "Singles" not in file:
            return None, None
        tree = file["Singles"]
        tree_np = tree.arrays(library="np")
        is_scatter = tree_np["compton_count"] + tree_np["rayleigh_count"] > 0
        keys = tree.keys()
        print(keys)
        print(
            max(tree_np["compton_count"]),
            min(tree_np["compton_count"]),
            tree_np["compton_count"].mean(),
        )
        print(
            max(tree_np["rayleigh_count"]),
            min(tree_np["rayleigh_count"]),
            tree_np["rayleigh_count"].mean(),
        )
        x_scatter = tree_np["PostPosition_X"][is_scatter]
        z_scatter = tree_np["PostPosition_Z"][is_scatter]
        x_primary = tree_np["PostPosition_X"][~is_scatter]
        z_primary = tree_np["PostPosition_Z"][~is_scatter]
        prj_primary, _, _ = np.histogram2d(
            z_primary, x_primary, bins=100, range=[[-75, 75], [-75, 75]]
        )
        prj_scatter, _, _ = np.histogram2d(
            z_scatter, x_scatter, bins=100, range=[[-75, 75], [-75, 75]]
        )
    return prj_primary, prj_scatter


def generate_prj_generic():
    import glob
    import SimpleITK as sitk

    root_list = glob.glob("generic_*_singles.root")
    output_prj_primary = np.zeros((100, 100))
    output_prj_scatter = np.zeros((100, 100))
    for ifile, filepath in enumerate(root_list):
        prj_primary, prj_scatter = generate_prj_one(filepath)
        if prj_primary is None or prj_scatter is None:
            continue
        output_prj_primary += prj_primary
        output_prj_scatter += prj_scatter
        if ifile > 1:
            break
    image_output_primary = sitk.GetImageFromArray(output_prj_primary.astype(np.float32))
    image_output_scatter = sitk.GetImageFromArray(output_prj_scatter.astype(np.float32))
    image_output_primary.SetSpacing([1.5, 1.5])
    image_output_scatter.SetSpacing([1.5, 1.5])
    sitk.WriteImage(image_output_primary, "generic_primary.mhd")
    sitk.WriteImage(image_output_scatter, "generic_scatter.mhd")


if __name__ == "__main__":
    pathFile = pathlib.Path(__file__).parent.resolve()
    # run_generic_source(50)
    # generate_prj_generic()
    # run_window_turbo_source(True)
    import sys

    run_window_turbo_source(bool(int(sys.argv[1])), int(sys.argv[2]))

    prj_primary, prj_scatter = generate_prj_one(
        paths.output / "window_turbo_ff_False_singles.root"
    )
    image_primary = sitk.GetImageFromArray(prj_primary.astype(np.float32))
    image_scatter = sitk.GetImageFromArray(prj_scatter.astype(np.float32))
    image_primary.SetSpacing([1.5, 1.5])
    image_scatter.SetSpacing([1.5, 1.5])
    sitk.WriteImage(
        image_primary, paths.output / "window_turbo_ff_False_primary_lazy_generated.mhd"
    )
    # sitk.WriteImage(image_scatter, paths.output / "window_turbo_ff_False_scatter.mhd")

    # analyze_root(paths.output / "window_turbo_ff_singles.root")

    # run_window_turbo_source()
    # profile_wt = calculate_profile(paths.output / "window_turbo_counts.mhd")
    # profile_generic = calculate_profile(paths.output_ref / "generic.mhd")
    # compare_result = compare_profiles(
    #     profile_generic,
    #     profile_wt,
    #     tolerance=4.0,
    #     fig_name=paths.output / "profile_comparison.png",
    # )

    # utility.test_ok(compare_result)
