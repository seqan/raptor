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

declare_datasource (FILE bin1.fa.gz
                    URL ${CMAKE_SOURCE_DIR}/test/data/bin1.fa.gz
                    URL_HASH SHA256=8c18f83a6e5c3fcde122052ceb843beadc7eaf31157b68b7a7caa70004be86ab)
declare_datasource (FILE bin2.fa.gz
                    URL ${CMAKE_SOURCE_DIR}/test/data/bin2.fa.gz
                    URL_HASH SHA256=60534c9bb4b593d39cedd5039614eb7eec6ab92961c490bb461232d9626d71a9)
declare_datasource (FILE bin3.fa.gz
                    URL ${CMAKE_SOURCE_DIR}/test/data/bin3.fa.gz
                    URL_HASH SHA256=1139affca816f0b404011a038d8c9d94481e28f84e566a57da4be4abc0b9eb13)
declare_datasource (FILE bin4.fa.gz
                    URL ${CMAKE_SOURCE_DIR}/test/data/bin4.fa.gz
                    URL_HASH SHA256=9986c8665294cc4a5eaf1341990f8680afd2ffa1aa6b4906d269f7ca96c393f2)

declare_datasource (FILE query.fq
                    URL ${CMAKE_SOURCE_DIR}/test/data/query.fq
                    URL_HASH SHA256=f48eb3f357e23df89e7e15d2f77f9285a428f73c4903eb1c6580271e0dea3d87)
declare_datasource (FILE query_empty.fq
                    URL ${CMAKE_SOURCE_DIR}/test/data/query_empty.fq
                    URL_HASH SHA256=80cd7628b6fdcb7dbbe99dd7f05a7a5578b462db9dd67b6d0e6982b03c5d4cd5)
declare_datasource (FILE query_socks.fq
                    URL ${CMAKE_SOURCE_DIR}/test/data/query_socks.fq
                    URL_HASH SHA256=c21cfd169d5b16192486b3175fb33f09eb3e0ae1154cde4aa102b62fc12b6683)

declare_datasource (FILE 1bins19window.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/1bins19window.index
                    URL_HASH SHA256=c25c3d0679c0cda7fa71841765c6aa3371227a63fd12401039627abf8680698c)
declare_datasource (FILE 1bins23window.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/1bins23window.index
                    URL_HASH SHA256=28bc4d11c0badbbd577919836d76dafd748f2642013b840eee4808e76f956081)
declare_datasource (FILE 64bins19window.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/64bins19window.index
                    URL_HASH SHA256=9c3f56c8907d96ab5f2439fcbbc7ca369f1c01aaa08e2f7a354687c6565c1c84)
declare_datasource (FILE 64bins23window.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/64bins23window.index
                    URL_HASH SHA256=87ac0febfd3b16ad67c8eab11beb52180cd545a7ec2a781f0d02f6df832a2e7b)
declare_datasource (FILE 128bins19window.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/128bins19window.index
                    URL_HASH SHA256=4890292fe5379227306774bf59bfb49acde789889296da83edf59001ada7677c)
declare_datasource (FILE 128bins23window.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/128bins23window.index
                    URL_HASH SHA256=558351b046f236877e94f31972933f0e9526ff2a6510215091fe4bfb1ae72157)
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

declare_datasource (FILE 1_1c.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/1_1c.index
                    URL_HASH SHA256=361172918e940944a4fa5ef518d3cc7ec8f024964a99fbfbf1cc6ce1132f637c)
declare_datasource (FILE 64bins23windowc.index
                    URL ${CMAKE_SOURCE_DIR}/test/data/64bins23windowc.index
                    URL_HASH SHA256=fa0cd03253dc8fd9620190a95f8db021c3ffb4db3a5b99ea754b5793879f3cd6)

declare_datasource (FILE 1bins19window.hibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/1bins19window.hibf
                    URL_HASH SHA256=e7d259809af8db71fd78161010e8c822880448ccf3d672f6974532a9cbf4d169)
declare_datasource (FILE 64bins19window.hibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/64bins19window.hibf
                    URL_HASH SHA256=354c85086bbbb497206b3fe50832435d9fe8877e8147eb6f57c1a6ef1f085a52)
declare_datasource (FILE 128bins19window.hibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/128bins19window.hibf
                    URL_HASH SHA256=6d1c4614559d3625a4625dac29627e8b7a6b3e83cf4fc7addde722f95054cebe)
declare_datasource (FILE 1bins23window.hibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/1bins23window.hibf
                    URL_HASH SHA256=56d3a9839230a2c9da943d5377ef1dfa087a0ed7930249aecdd8a6a64b6e1079)
declare_datasource (FILE 64bins23window.hibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/64bins23window.hibf
                    URL_HASH SHA256=d9a95119322915f33faad70ee5b7d39f7815170f609aeebc37d44451482516b6)
declare_datasource (FILE 128bins23window.hibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/128bins23window.hibf
                    URL_HASH SHA256=bdb5f09f0357d577cf1db68145597f97b01d47bb7a1ef859a26cbccd93649939)

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

declare_datasource (FILE 1bins19window0errorhibf.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/1bins19window0errorhibf.out
                    URL_HASH SHA256=db866e80d7327cf7bfb792ef9c6af425354b14e8ebd7058796cec9bb43672525)
declare_datasource (FILE 1bins19window1errorhibf.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/1bins19window1errorhibf.out
                    URL_HASH SHA256=db866e80d7327cf7bfb792ef9c6af425354b14e8ebd7058796cec9bb43672525)
declare_datasource (FILE 64bins19window0errorhibf.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/64bins19window0errorhibf.out
                    URL_HASH SHA256=8fcdd69d9447d1b4020139c046ef838059c9c6f153260658aee67dfae6b89659)
declare_datasource (FILE 64bins19window1errorhibf.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/64bins19window1errorhibf.out
                    URL_HASH SHA256=61ace76230a97278875d8452dd550e7c853c1d042932d2f0d14a1b0ab63e4ff4)
declare_datasource (FILE 128bins19window0errorhibf.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/128bins19window0errorhibf.out
                    URL_HASH SHA256=bf513dc1311a81f6bdb246406c069facf672a880dddff8b48fa2337b7325d19c)
declare_datasource (FILE 128bins19window1errorhibf.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/128bins19window1errorhibf.out
                    URL_HASH SHA256=bf513dc1311a81f6bdb246406c069facf672a880dddff8b48fa2337b7325d19c)

declare_datasource (FILE 1binsempty.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/1binsempty.out
                    URL_HASH SHA256=515b16f2e507b9f4bac1b717c65d3a16b7c0f2bebe0977f7e0a471479dbcff42)
declare_datasource (FILE 64binsempty.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/64binsempty.out
                    URL_HASH SHA256=3ca3257deccc01931390e3b0185d26ab97775a022eb8696930d680ba8701454e)
declare_datasource (FILE 128binsempty.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/128binsempty.out
                    URL_HASH SHA256=e4eee8bce310b11bdb4da718fe3f3129508955b35cb89c7dc0f61b952d48410b)

declare_datasource (FILE 1binsemptyhibf.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/1binsemptyhibf.out
                    URL_HASH SHA256=db866e80d7327cf7bfb792ef9c6af425354b14e8ebd7058796cec9bb43672525)
declare_datasource (FILE 64binsemptyhibf.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/64binsemptyhibf.out
                    URL_HASH SHA256=46ccd96c51f14c30661cd63c0ad370323c429581c877a67d199b72b832398b5b)
declare_datasource (FILE 128binsemptyhibf.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/hibf/128binsemptyhibf.out
                    URL_HASH SHA256=db4af85cad7f68db24a65663015770df250e56becad66044b02657d7f989807e)

configure_datasource (FILE 1bins.pack
                      URL ${CMAKE_SOURCE_DIR}/test/data/hibf/1bins.pack
                      URL_HASH SHA256=1b7c15eff5a2d884b0058d3ecdc2155f0de0cdc2e59e7d9598030eaaecda8464)
configure_datasource (FILE 64bins.pack
                      URL ${CMAKE_SOURCE_DIR}/test/data/hibf/64bins.pack
                      URL_HASH SHA256=f33eff072b77f47e755ddce2472e7e171e63a0a8192ccaf9fc797c97921226ef)
configure_datasource (FILE 128bins.pack
                      URL ${CMAKE_SOURCE_DIR}/test/data/hibf/128bins.pack
                      URL_HASH SHA256=8dcaf32047203556cb9652628a294ffc3faee1b60a5c887687c9cb52e9ca4a91)
