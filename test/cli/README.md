# CLI Test

Here are test files for command line interface tests, i.e. the app is executed with defined input files and parameters.
The test then validates whether the output is correct.

Each test fixture should be inherited from the `cli_test` class: It provides the functionality of executing the app,
finding the input files, capturing the output and creating individual test directories.
We provide a new test macro `EXPECT_RANGE_EQ` for the comparison of whole ranges.
It provides more insights about differences than an implementation with `EXPECT_TRUE(std::ranges::equal())`.
The test output files are stored in the directory `test/output`.

To prevent issues when running multiple CLI tests in parallel, each CLI should use its own derived test fixture.
For example, if two CLI tests run a test called `my_test` and both use the `cli_test` test fixture, there may
be issues because both tests will use the same working directory. Instead, define different test fixtures for both
tests, e.g. `struct my_test_fixture : public cli_test {};`.

Attention: The default `make` target does not build tests.
Please invoke the build with `make cli_test` or use `make test` to build and run all kinds of tests.
