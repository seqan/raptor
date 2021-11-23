#!/usr/bin/env bash
# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

/example/runs/raptor_scripts/raptor.sh 0.05 4G 20 20
sleep 1
/example/runs/raptor_scripts/raptor.sh 0.015 8G 20 20
sleep 1
/example/runs/raptor_scripts/raptor.sh 0.05 1584M 23 19
sleep 1
/example/runs/raptor_scripts/raptor.sh 0.015 3069M 23 19
sleep 1
/example/runs/raptor_scripts/raptor.sh 0.05 1584M 24 20
sleep 1
/example/runs/raptor_scripts/raptor.sh 0.015 3069M 24 20
sleep 1
