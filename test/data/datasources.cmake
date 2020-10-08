cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

declare_datasource (FILE example_data.tar.gz
                    URL https://ftp.imp.fu-berlin.de/pub/seiler/raptor/example_data.tar.gz
                    URL_HASH SHA256=7c2e7bdbf573cfe2314c8255080a5d57d966722e8bfb53712e03d87ce463ff15)


declare_datasource (FILE expected_results.tar.gz
                    URL https://ftp.imp.fu-berlin.de/pub/seiler/raptor/expected_results.tar.gz
                    URL_HASH SHA256=2685ef95ebea074514f4736888b493857f0327514684ef88d798b3f25df5fd5a)
