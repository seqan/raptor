cmake_minimum_required (VERSION 3.8)

add_cli_test (raptor_options_test.cpp)
target_use_datasources (raptor_options_test FILES example_data.tar.gz)

add_cli_test (raptor_test.cpp)
target_use_datasources (raptor_test FILES example_data.tar.gz)
target_use_datasources (raptor_test FILES expected_results.tar.gz)
