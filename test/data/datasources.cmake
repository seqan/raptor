cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

declare_datasource (FILE example_data.tar.gz
                    URL https://ftp.imp.fu-berlin.de/pub/seiler/raptor/example_data.tar.gz
                    URL_HASH SHA256=306fbaebb9723c23c215038092e950c95b8932b314730288d6b9baa3d52f9634)


declare_datasource (FILE expected_results.tar.gz
                    URL https://ftp.imp.fu-berlin.de/pub/seiler/raptor/expected_results.tar.gz
                    URL_HASH SHA256=2685ef95ebea074514f4736888b493857f0327514684ef88d798b3f25df5fd5a)
