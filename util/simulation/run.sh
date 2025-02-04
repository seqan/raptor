#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -Eeuo pipefail

SCRIPT_ROOT=$(dirname "$(readlink -f "$0")")
source "${SCRIPT_ROOT}/config.sh"

cleanup() {
    status=${status:-$?}
    set +x
    trap '' EXIT ABRT HUP INT PIPE QUIT TERM
    exec 2>&4 1>&3

    if [[ "$status" -ne 0 ]]; then
        echo "[ERROR] The log file can be found at ${__LOG_FILE}"
        echo "[ERROR] Last 10 lines of the log file follow"
        tail -n 10 "${__LOG_FILE}" | xargs -d '\n' -L 1 echo "    "
    else
        echo "[Success] The log file can be found at ${__LOG_FILE}"
        echo "[Success] The output can be found in ${OUTPUT_DIR}"
    fi

    if [[ -v __SIMULATION_TMP_DIR && -d ${__SIMULATION_TMP_DIR} ]]; then
        rm -fdr "${__SIMULATION_TMP_DIR}"
    fi

    exit "$status"
}
sig_cleanup() {
    status=$?
    set +x
    trap '' EXIT
    cleanup
}
trap cleanup EXIT
trap sig_cleanup ABRT HUP INT PIPE QUIT TERM

quiet_loop() {
    if [[ $- =~ x ]]; then
        set +x
    fi
}

__LOG_DIR=${OUTPUT_DIR}/logs
__BUILD_DIR=${OUTPUT_DIR}/build
__BINARY_DIR=${__BUILD_DIR}/bin
__NEEDS_MASON=$(test "${NUMBER_OF_HAPLOTYPES}" -ne 1 && echo ON || echo OFF)

mkdir -p "${__LOG_DIR}"
mkdir -p "${__BUILD_DIR}"
mkdir -p "${__BINARY_DIR}"

__LOG_FILE=${__LOG_DIR}/$(date +"%Y-%m-%d_%H-%M-%S").log
touch "${__LOG_FILE}"
echo "## Log file can be found at ${__LOG_FILE}"
exec 3>&1 4>&2 1>>"${__LOG_FILE}" 2>&1

set -x

echo "## [$(date +"%Y-%m-%d %T")] Building Dependencies" | tee /dev/fd/3
cd "${__BUILD_DIR}"
if [[ ! -d "${REPO_PATH}" && ! -d ${__BUILD_DIR}/raptor ]]; then
    wget -q -O raptor.tar.gz https://github.com/seqan/raptor/archive/refs/heads/main.tar.gz
    tar xf raptor.tar.gz && rm raptor.tar.gz
    REPO_PATH=${__BUILD_DIR}/raptor-main
fi

cmake "${REPO_PATH}/util/iScience" -DCMAKE_BUILD_TYPE=Release \
                                   -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
                                   -DRAPTOR_UTILITY_BUILD_MASON="${__NEEDS_MASON}" \
                                   -DINSTALL_RAPTOR=OFF \
                                   -DSHARG_NO_TDL=ON \
                                   -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
                                   -Wno-dev

__BUILD_TARGETS=(generate_reads simulate_sequence split_sequence)
if [[ "${__NEEDS_MASON}" == "ON" ]]; then
    __BUILD_TARGETS+=(mason)
fi
make -j"${THREADS}" --no-print-directory "${__BUILD_TARGETS[@]}"

set +x

for __NUMBER_OF_BINS in ${NUMBER_OF_BINS}; do
    expand_suffix "__NUMBER_OF_BINS"

    __SIMULATION_DIR=${OUTPUT_DIR}/${__NUMBER_OF_BINS}
    __SIMULATION_BIN_DIR=${__SIMULATION_DIR}/bins
    __SIMULATION_INFO_DIR=${__SIMULATION_DIR}/info
    __SIMULATION_TMP_DIR=${__SIMULATION_DIR}/TMP_$(date +"%Y-%m-%d_%H-%M-%S")
    __SIMULATION_SPLIT_EXTRA_ARGS=()
    if [[ -n "${SAMPLE_FROM}" ]]; then
        __SIMULATION_SPLIT_EXTRA_ARGS+=(--sample "${SAMPLE_FROM}")
    fi
    if [[ "${__NEEDS_MASON}" == "OFF" ]]; then
        __SIMULATION_SPLIT_EXTRA_ARGS+=(--output "${__SIMULATION_BIN_DIR}")
    fi
    __SIMULATION_SPLIT_EXTENSION=$(test "${__NEEDS_MASON}" == "ON"  && echo fasta || echo fa)

    mkdir -p "${__SIMULATION_DIR}"
    mkdir -p "${__SIMULATION_BIN_DIR}"
    mkdir -p "${__SIMULATION_INFO_DIR}"
    mkdir -p "${__SIMULATION_TMP_DIR}"

    set -x
    echo "## [$(date +"%Y-%m-%d %T")] Number of bins: ${__NUMBER_OF_BINS}" | tee /dev/fd/3

    echo "## [$(date +"%Y-%m-%d %T")]   Simulating ${__NUMBER_OF_BINS} bins with a total length of $(to_iec "${REFERENCE_LENGTH}")" | tee /dev/fd/3
    # Simulate reference
    "${__BINARY_DIR}"/simulate_sequence \
        --length "${REFERENCE_LENGTH}" \
        --output "${__SIMULATION_TMP_DIR}"/ref.fasta \
        --seed "${SEED}"

    # Evenly distribute it over bins
    "${__BINARY_DIR}"/split_sequence \
        --input "${__SIMULATION_TMP_DIR}"/ref.fasta \
        --parts "${__NUMBER_OF_BINS}" \
        "${__SIMULATION_SPLIT_EXTRA_ARGS[@]}"

    # We do not need the reference anymore
    rm "${__SIMULATION_TMP_DIR}"/ref.fasta

    if [[ "${__NEEDS_MASON}" == "ON"  ]]; then
        echo "## [$(date +"%Y-%m-%d %T")]   Simulating $(to_iec "${NUMBER_OF_HAPLOTYPES}") haplotypes for each bin" | tee /dev/fd/3
        # Simulate haplotypes for each bin
        for file in "${__SIMULATION_TMP_DIR}"/*.fa
        do
            (
                "${__BINARY_DIR}"/mason_variator \
                    --in-reference "${file}" \
                    --num-haplotypes "${NUMBER_OF_HAPLOTYPES}" \
                    --out-fasta "${__SIMULATION_BIN_DIR}"/"$(basename "${file}" .fa)".fasta \
                    --out-vcf "${__SIMULATION_INFO_DIR}"/"$(basename "${file}" .fa)".vcf \
                    --quiet &>/dev/null
                if [[ "${KEEP_VCF}" == "false" ]]; then
                    rm "${__SIMULATION_INFO_DIR}"/"$(basename "${file}" .fa)".vcf
                else
                    gzip --force "${__SIMULATION_INFO_DIR}"/"$(basename "${file}" .fa)".vcf
                fi
            ) &

            if [[ $(jobs -r -p | wc -l) -ge ${THREADS} ]]; then
                wait -n
            fi
            quiet_loop
        done
        wait
        set -x
    fi

    echo "## [$(date +"%Y-%m-%d %T")]   Simulating $(to_iec "${READ_COUNT}") reads of length ${READ_LENGTH} with ${READ_ERRORS} errors" | tee /dev/fd/3
    __LIST_FILE="${__SIMULATION_DIR}"/filenames.txt
    find "${__SIMULATION_BIN_DIR}" -type f -name "*.${__SIMULATION_SPLIT_EXTENSION}" | sort -V > "${__LIST_FILE}"
    "${__BINARY_DIR}"/generate_reads \
        --input "${__LIST_FILE}" \
        --output "${__SIMULATION_DIR}/reads_e${READ_ERRORS}_${READ_LENGTH}.fastq" \
        --errors "${READ_ERRORS}" \
        --number_of_reads "${READ_COUNT}" \
        --read_length "${READ_LENGTH}" \
        --weights from_file_sizes
    set -x

    rm -fdr "${__SIMULATION_TMP_DIR}"

    cp --no-preserve=mode "${SCRIPT_ROOT}"/config.sh "${__SIMULATION_INFO_DIR}"
    cp --no-preserve=mode "${SCRIPT_ROOT}"/"$(basename "$0")" "${__SIMULATION_INFO_DIR}"

    set +x
done
