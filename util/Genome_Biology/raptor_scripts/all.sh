#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -Eeuo pipefail

# Parameters:
# FPR
# (Vanilla) IBF size
# window size
# k-mer size
# cleanup?

/example/runs/raptor_scripts/raptor.sh 0.05 4G 20 20 false
sleep 1
/example/runs/raptor_scripts/raptor.sh 0.015 8G 20 20 false
sleep 1
/example/runs/raptor_scripts/raptor.sh 0.05 1584M 23 19 false
sleep 1
/example/runs/raptor_scripts/raptor.sh 0.015 3069M 23 19 false
sleep 1
/example/runs/raptor_scripts/raptor.sh 0.05 1584M 24 20 false
sleep 1
/example/runs/raptor_scripts/raptor.sh 0.015 3069M 24 20 true
sleep 1
