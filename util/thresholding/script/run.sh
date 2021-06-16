#!/usr/bin/env bash

# Exit on error
set -e

##################################################################
##################################################################
# Place this file in the build directory.
##################################################################
##################################################################

##################################################################
# Set pattern sizes and taus to run for.
##################################################################
pattern_sizes=( 100 150 250 )
tau_values=( $(seq 0 0.05 1.01) 0.99 0.999 )

##################################################################
# Create directories.
##################################################################
for p in ${pattern_sizes[*]}
do
    mkdir -p test_data/$p
    mkdir -p test_out/negative/$p
    mkdir -p test_out/positive/$p
    mkdir -p test_out/with_precompute/negative/$p
    mkdir -p test_out/with_precompute/positive/$p
done

##################################################################
# Simulate random data.
##################################################################
for p in ${pattern_sizes[*]}
do
	./bin/random_data --query_length $p --min_error 3 --max_error 3 --out test_data/$p
done

##################################################################
# Search for random query.
##################################################################
for p in ${pattern_sizes[*]}
do
    for tau in ${tau_values[*]}
    do ./bin/query --reference test_data/$p/reference.fasta \
               --query test_data/$p/query_random.fasta \
               --out test_out/negative/$p \
               --error 3 \
               --method all \
               --tau $tau
    done
done

##################################################################
# Copy precomputed thresholds.
##################################################################
for p in ${pattern_sizes[*]}
do
	cp test_out/negative/$p/binary_* test_out/with_precompute/negative/$p
done

##################################################################
# Search for random query using the precomputed thresholds.
##################################################################
for p in ${pattern_sizes[*]}
do
    for tau in ${tau_values[*]}
    do ./bin/query --reference test_data/$p/reference.fasta \
               --query test_data/$p/query_random.fasta \
               --out test_out/with_precompute/negative/$p \
               --error 3 \
               --method all \
               --from_file \
               --tau $tau
    done
done

##################################################################
# Search for query with 3 errors.
##################################################################
for p in ${pattern_sizes[*]}
do
    for tau in ${tau_values[*]}
    do ./bin/query --reference test_data/$p/reference.fasta \
               --query test_data/$p/query_e3.fasta \
               --out test_out/positive/$p \
               --error 3 \
               --method all \
               --tau $tau
    done
done

##################################################################
# Copy precomputed thresholds.
##################################################################
for p in ${pattern_sizes[*]}
do
	cp test_out/positive/$p/binary_* test_out/with_precompute/positive/$p
done

##################################################################
# Search for query with 3 errors using the precomputed thresholds.
##################################################################
for p in ${pattern_sizes[*]}
do
    for tau in ${tau_values[*]}
    do ./bin/query --reference test_data/$p/reference.fasta \
               --query test_data/$p/query_random.fasta \
               --out test_out/with_precompute/positive/$p \
               --error 3 \
               --method all \
               --from_file \
               --tau $tau
    done
done
