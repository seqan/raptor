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

declare_datasource (FILE 1bins19window.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/1bins19window.index
                    URL_HASH SHA256=ecf2a165fb4c1c1602305dad13aecb084aa6f414fea58aafdf6b52b0e4770181)
declare_datasource (FILE 1bins23window.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/1bins23window.index
                    URL_HASH SHA256=e5d88a3d18d8e9193511a87d147d720e1996397468345caaa74b2ae2debfaa55)
declare_datasource (FILE 64bins19window.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/64bins19window.index
                    URL_HASH SHA256=d333b5a7b93560ec748fdf3d22087c283d43ddbaacca6eb612d13f04d9be5876)
declare_datasource (FILE 64bins23window.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/64bins23window.index
                    URL_HASH SHA256=8bf6c1e91de4401383f299e00a2803cf1256a63048bbbd68c1086f92750b79b7)
declare_datasource (FILE 128bins19window.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/128bins19window.index
                    URL_HASH SHA256=fa0d51f25037315926fbf2d291af2e3ad95c0508630b5a8ed1e291b45324e8c4)
declare_datasource (FILE 128bins23window.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/128bins23window.index
                    URL_HASH SHA256=c8d63563e8e301a21e19a46605fac354fa4d1b1b0b5e7d5a6a066fa024e05e62)
declare_datasource (FILE 1_1.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/1_1.index
                    URL_HASH SHA256=2a00d9fd2cf9865841b559219cb1896d33e8a20e77604906ba551e94fbdbed7f)
declare_datasource (FILE 1_1.index_0
                    URL ${CMAKE_SOURCE_DIR}/test/data/1_1.index_0
                    URL_HASH SHA256=b39aaf45fab2837903add389b3d3efe33557808f8dd823113efbbd4df1e601f3)
declare_datasource (FILE 1_1.index_1
                    URL ${CMAKE_SOURCE_DIR}/test/data/1_1.index_1
                    URL_HASH SHA256=1a882db0e38056ae1b63efdb0cfd49b46d8c26d7a2edcaf62ca27115487bd682)
declare_datasource (FILE 1_1.index_2
                    URL ${CMAKE_SOURCE_DIR}/test/data/1_1.index_2
                    URL_HASH SHA256=1f4cff9952cdc6137500d89b3a3b2e11eef33734b858be9e328dfe4c03d8c4c0)
declare_datasource (FILE 1_1.index_3
                    URL ${CMAKE_SOURCE_DIR}/test/data/1_1.index_3
                    URL_HASH SHA256=b72db8be6455990c4a67c22748b725d501e400f61e5291293c11cd6c3793f93a)

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
