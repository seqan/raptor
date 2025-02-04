# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# Provides the `declare_datasource` function.
include (declare_datasource)

# Makes all files in this directory and subdirectories available to the build, i.e., copying them to <build>/data/.
# `datasources.cmake`, `README.md`, and `REUSE.toml` are ignored.
# You may organise your data in subdirectories, but each file must have a unique name.
include (add_local_data)
