# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25)

project (raptor_test_performance LANGUAGES CXX C)

include (../raptor-test.cmake)

CPMGetPackage (benchmark)

raptor_add_benchmark (bin_influence_benchmark.cpp)
