cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# declare_datasource (FILE example_data.tar.gz
#                     URL https://ftp.imp.fu-berlin.de/pub/seiler/raptor/example_data.tar.gz
#                     URL_HASH SHA256=7c2e7bdbf573cfe2314c8255080a5d57d966722e8bfb53712e03d87ce463ff15)


# declare_datasource (FILE expected_results.tar.gz
#                     URL https://ftp.imp.fu-berlin.de/pub/seiler/raptor/expected_results.tar.gz
#                     URL_HASH SHA256=2685ef95ebea074514f4736888b493857f0327514684ef88d798b3f25df5fd5a)

declare_datasource (FILE bin1.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/bin1.fa
                    URL_HASH SHA256=03ac546c2dfa34f2a75283d43dd2293aa702ed2b7dd877dfbb4f966d258201f0)
declare_datasource (FILE bin2.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/bin2.fa
                    URL_HASH SHA256=5bbf2f6fee427d3e967932203987b8c391ac86f0f089101455d529cb7e27b985)
declare_datasource (FILE bin3.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/bin3.fa
                    URL_HASH SHA256=93d80997611587f9635bdf7ed48eff26d995b54899e757026551df90e37464e7)
declare_datasource (FILE bin4.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/bin4.fa
                    URL_HASH SHA256=b4e3be3c001aac8481be2251aab110e4318996d26f2cdc91fb3061c35cbec4a1)

declare_datasource (FILE query.fq
                    URL ${CMAKE_SOURCE_DIR}/test/data/query.fq
                    URL_HASH SHA256=f48eb3f357e23df89e7e15d2f77f9285a428f73c4903eb1c6580271e0dea3d87)
declare_datasource (FILE query_socks.fq
                    URL ${CMAKE_SOURCE_DIR}/test/data/query_socks.fq
                    URL_HASH SHA256=c21cfd169d5b16192486b3175fb33f09eb3e0ae1154cde4aa102b62fc12b6683)

declare_datasource (FILE 1bins19window.ibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/1bins19window.ibf
                    URL_HASH SHA256=799688b64a255f4b4ee3a3b7e8388865f8c6df899b307b5f596a8c65ea60c7ec)
declare_datasource (FILE 1bins23window.ibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/1bins23window.ibf
                    URL_HASH SHA256=f20d6526b14be258dfb6b56cab884783f929cc89a8ca3f03c4522ba415a52561)
declare_datasource (FILE 64bins19window.ibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/64bins19window.ibf
                    URL_HASH SHA256=ae8f2480da5c14658b445e283242b2a3235ad79a4df3017d6e64855b1f6d5683)
declare_datasource (FILE 64bins23window.ibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/64bins23window.ibf
                    URL_HASH SHA256=8604caf9db0066ec615b9582296463fbac1d4616673d8c8184b4062f9d9c9c28)
declare_datasource (FILE 128bins19window.ibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/128bins19window.ibf
                    URL_HASH SHA256=dd6e460cfefb55af6047952869cdc62f06ca67f5cbf046e182a41c184c509887)
declare_datasource (FILE 128bins23window.ibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/128bins23window.ibf
                    URL_HASH SHA256=52b145747bb1404efccd9387197dadf3083d3fbd06aa2e19b21a7364be830ce7)

declare_datasource (FILE 1bins19window0error.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/1bins19window0error.out
                    URL_HASH SHA256=1ff27bf9de0ad5e5487c857327c2494fb42c38e81cf8f86291157cf2c941af86)
declare_datasource (FILE 1bins19window1error.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/1bins19window1error.out
                    URL_HASH SHA256=1ff27bf9de0ad5e5487c857327c2494fb42c38e81cf8f86291157cf2c941af86)
declare_datasource (FILE 1bins23window1error.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/1bins23window1error.out
                    URL_HASH SHA256=1ff27bf9de0ad5e5487c857327c2494fb42c38e81cf8f86291157cf2c941af86)
declare_datasource (FILE 64bins19window0error.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/64bins19window0error.out
                    URL_HASH SHA256=6a8a88f7cba0aae01f8cde75e875753c196704d36cf7430ae9c385ff7843b137)
declare_datasource (FILE 64bins19window1error.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/64bins19window1error.out
                    URL_HASH SHA256=3b65eb1cda3a36b30412931f8eb6e167789d3f97ff45d52d40ccd692dfb6b922)
declare_datasource (FILE 64bins23window1error.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/64bins23window1error.out
                    URL_HASH SHA256=3b65eb1cda3a36b30412931f8eb6e167789d3f97ff45d52d40ccd692dfb6b922)
declare_datasource (FILE 128bins19window0error.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/128bins19window0error.out
                    URL_HASH SHA256=ed7832ba7a7e945e5501aa86520daf2b3dce013a3c67f515faa701f639acc740)
declare_datasource (FILE 128bins19window1error.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/128bins19window1error.out
                    URL_HASH SHA256=a331b0c4ed4947ffa7fa44c1b99d50803d2b1e289365594286b273bd57fe45dd)
declare_datasource (FILE 128bins23window1error.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/128bins23window1error.out
                    URL_HASH SHA256=a331b0c4ed4947ffa7fa44c1b99d50803d2b1e289365594286b273bd57fe45dd)

declare_datasource (FILE 1bins19window0errorsocks.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/1bins19window0errorsocks.out
                    URL_HASH SHA256=b1053f5b749071c54528f7b00cc66f278eb50f5643d46db01774932ea510fbbd)
declare_datasource (FILE 64bins19window0errorsocks.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/64bins19window0errorsocks.out
                    URL_HASH SHA256=d5e883ce66d1fc574cd6c9ae62c2d25a7671ad7be0948c49d969b1de88293def)
declare_datasource (FILE 128bins19window0errorsocks.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/128bins19window0errorsocks.out
                    URL_HASH SHA256=3087b98c1f4e0d0a1a54f0dc11100f25fdfa8eba62d8f24c082ea1c523146618)
