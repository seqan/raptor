# Coverage Test

This is the test for the code coverage.
It reaches 100% if each code line is executed at least once through the app tests.

Coverage tests are not yet part of the CMake configuration like it is in SeqAn3.
However, the CI tests code coverage by running gcov and uploading the results to codecov.

The custom build type `Coverage` can be used to configure the project such that coverage reports are possible.
