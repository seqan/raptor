cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

declare_datasource (FILE example_data.tar.gz
                    URL https://ftp.imp.fu-berlin.de/pub/seiler/raptor/example_data.tar.gz
                    URL_HASH SHA256=8f1f87adcda881549b8f5c8072c4e2d04c7739a4641644a885dd99b11b2a9179)


declare_datasource (FILE expected_results.tar.gz
                    URL https://ftp.imp.fu-berlin.de/pub/seiler/raptor/expected_results.tar.gz
                    URL_HASH SHA256=2685ef95ebea074514f4736888b493857f0327514684ef88d798b3f25df5fd5a)
