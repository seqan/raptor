cmake_minimum_required (VERSION 3.11)

include (cmake/app_internal_datasources.cmake)

declare_internal_datasource (FILE bin1.fa
                             URL ${CMAKE_SOURCE_DIR}/test/data/bin1.fa
                             URL_HASH SHA256=03ac546c2dfa34f2a75283d43dd2293aa702ed2b7dd877dfbb4f966d258201f0)
declare_internal_datasource (FILE bin2.fa
                             URL ${CMAKE_SOURCE_DIR}/test/data/bin2.fa
                             URL_HASH SHA256=5bbf2f6fee427d3e967932203987b8c391ac86f0f089101455d529cb7e27b985)
declare_internal_datasource (FILE bin3.fa
                             URL ${CMAKE_SOURCE_DIR}/test/data/bin3.fa
                             URL_HASH SHA256=93d80997611587f9635bdf7ed48eff26d995b54899e757026551df90e37464e7)
declare_internal_datasource (FILE bin4.fa
                             URL ${CMAKE_SOURCE_DIR}/test/data/bin4.fa
                             URL_HASH SHA256=b4e3be3c001aac8481be2251aab110e4318996d26f2cdc91fb3061c35cbec4a1)

declare_internal_datasource (FILE bin1.fa.gz
                             URL ${CMAKE_SOURCE_DIR}/test/data/bin1.fa.gz
                             URL_HASH SHA256=8c18f83a6e5c3fcde122052ceb843beadc7eaf31157b68b7a7caa70004be86ab)
declare_internal_datasource (FILE bin2.fa.gz
                             URL ${CMAKE_SOURCE_DIR}/test/data/bin2.fa.gz
                             URL_HASH SHA256=60534c9bb4b593d39cedd5039614eb7eec6ab92961c490bb461232d9626d71a9)
declare_internal_datasource (FILE bin3.fa.gz
                             URL ${CMAKE_SOURCE_DIR}/test/data/bin3.fa.gz
                             URL_HASH SHA256=1139affca816f0b404011a038d8c9d94481e28f84e566a57da4be4abc0b9eb13)
declare_internal_datasource (FILE bin4.fa.gz
                             URL ${CMAKE_SOURCE_DIR}/test/data/bin4.fa.gz
                             URL_HASH SHA256=9986c8665294cc4a5eaf1341990f8680afd2ffa1aa6b4906d269f7ca96c393f2)

declare_internal_datasource (FILE query.fq
                             URL ${CMAKE_SOURCE_DIR}/test/data/query.fq
                             URL_HASH SHA256=f48eb3f357e23df89e7e15d2f77f9285a428f73c4903eb1c6580271e0dea3d87)
declare_internal_datasource (FILE query_empty.fq
                             URL ${CMAKE_SOURCE_DIR}/test/data/query_empty.fq
                             URL_HASH SHA256=80cd7628b6fdcb7dbbe99dd7f05a7a5578b462db9dd67b6d0e6982b03c5d4cd5)
declare_internal_datasource (FILE query_socks.fq
                             URL ${CMAKE_SOURCE_DIR}/test/data/query_socks.fq
                             URL_HASH SHA256=c21cfd169d5b16192486b3175fb33f09eb3e0ae1154cde4aa102b62fc12b6683)

declare_internal_datasource (FILE 1bins19window.index
                             URL ${CMAKE_SOURCE_DIR}/test/data/1bins19window.index
                             URL_HASH SHA256=c25c3d0679c0cda7fa71841765c6aa3371227a63fd12401039627abf8680698c)
declare_internal_datasource (FILE 1bins23window.index
                             URL ${CMAKE_SOURCE_DIR}/test/data/1bins23window.index
                             URL_HASH SHA256=28bc4d11c0badbbd577919836d76dafd748f2642013b840eee4808e76f956081)
declare_internal_datasource (FILE 64bins19window.index
                             URL ${CMAKE_SOURCE_DIR}/test/data/64bins19window.index
                             URL_HASH SHA256=9c3f56c8907d96ab5f2439fcbbc7ca369f1c01aaa08e2f7a354687c6565c1c84)
declare_internal_datasource (FILE 64bins23window.index
                             URL ${CMAKE_SOURCE_DIR}/test/data/64bins23window.index
                             URL_HASH SHA256=87ac0febfd3b16ad67c8eab11beb52180cd545a7ec2a781f0d02f6df832a2e7b)
declare_internal_datasource (FILE 128bins19window.index
                             URL ${CMAKE_SOURCE_DIR}/test/data/128bins19window.index
                             URL_HASH SHA256=4890292fe5379227306774bf59bfb49acde789889296da83edf59001ada7677c)
declare_internal_datasource (FILE 128bins23window.index
                             URL ${CMAKE_SOURCE_DIR}/test/data/128bins23window.index
                             URL_HASH SHA256=558351b046f236877e94f31972933f0e9526ff2a6510215091fe4bfb1ae72157)
declare_internal_datasource (FILE 1_1.index
                             URL ${CMAKE_SOURCE_DIR}/test/data/1_1.index
                             URL_HASH SHA256=2a00d9fd2cf9865841b559219cb1896d33e8a20e77604906ba551e94fbdbed7f)
declare_internal_datasource (FILE 1_1.index_0
                             URL ${CMAKE_SOURCE_DIR}/test/data/1_1.index_0
                             URL_HASH SHA256=b39aaf45fab2837903add389b3d3efe33557808f8dd823113efbbd4df1e601f3)
declare_internal_datasource (FILE 1_1.index_1
                             URL ${CMAKE_SOURCE_DIR}/test/data/1_1.index_1
                             URL_HASH SHA256=1a882db0e38056ae1b63efdb0cfd49b46d8c26d7a2edcaf62ca27115487bd682)
declare_internal_datasource (FILE 1_1.index_2
                             URL ${CMAKE_SOURCE_DIR}/test/data/1_1.index_2
                             URL_HASH SHA256=1f4cff9952cdc6137500d89b3a3b2e11eef33734b858be9e328dfe4c03d8c4c0)
declare_internal_datasource (FILE 1_1.index_3
                             URL ${CMAKE_SOURCE_DIR}/test/data/1_1.index_3
                             URL_HASH SHA256=b72db8be6455990c4a67c22748b725d501e400f61e5291293c11cd6c3793f93a)

declare_internal_datasource (FILE 1_1c.index
                             URL ${CMAKE_SOURCE_DIR}/test/data/1_1c.index
                             URL_HASH SHA256=361172918e940944a4fa5ef518d3cc7ec8f024964a99fbfbf1cc6ce1132f637c)
declare_internal_datasource (FILE 64bins23windowc.index
                             URL ${CMAKE_SOURCE_DIR}/test/data/64bins23windowc.index
                             URL_HASH SHA256=fa0cd03253dc8fd9620190a95f8db021c3ffb4db3a5b99ea754b5793879f3cd6)

declare_internal_datasource (FILE 1bins19window.hibf
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/1bins19window.hibf
                             URL_HASH SHA256=c72e143f7e205fb0eb536ee9db57b3a94b04f2b02fbd2f362919ab714c93994b)
declare_internal_datasource (FILE 64bins19window.hibf
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/64bins19window.hibf
                             URL_HASH SHA256=4604008953ded1b79236bb637a5ba1c29ebe464a874c438c9f2c78c08d5f1083)
declare_internal_datasource (FILE 128bins19window.hibf
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/128bins19window.hibf
                             URL_HASH SHA256=44c15ce5afdd677f7552dfce51d32ebcfce41577836ae128985f6462aea83677)
declare_internal_datasource (FILE 1bins23window.hibf
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/1bins23window.hibf
                             URL_HASH SHA256=01d2e55ab34f8fab023d3a96285c5ec54fffc58fb415a837f7bfdd11ff73913e)
declare_internal_datasource (FILE 64bins23window.hibf
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/64bins23window.hibf
                             URL_HASH SHA256=ce106cf0f93f05584f1ff07663b0a3eec8717f0ca76068233b0edef02a6b459d)
declare_internal_datasource (FILE 128bins23window.hibf
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/128bins23window.hibf
                             URL_HASH SHA256=e12bacafac51a4517c765effb1332b7964c7b22a7d98361847abd92056da296f)
declare_internal_datasource (FILE three_levels.hibf
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/three_levels.hibf
                             URL_HASH SHA256=f63ed7dae469e5b821b0f9285e7a7ac61fe29fc51ec6e3e8808789bf155cc6ee)

declare_internal_datasource (FILE 1bins19window0error.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/1bins19window0error.out
                             URL_HASH SHA256=1ff27bf9de0ad5e5487c857327c2494fb42c38e81cf8f86291157cf2c941af86)
declare_internal_datasource (FILE 1bins19window1error.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/1bins19window1error.out
                             URL_HASH SHA256=1ff27bf9de0ad5e5487c857327c2494fb42c38e81cf8f86291157cf2c941af86)
declare_internal_datasource (FILE 1bins23window1error.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/1bins23window1error.out
                             URL_HASH SHA256=1ff27bf9de0ad5e5487c857327c2494fb42c38e81cf8f86291157cf2c941af86)
declare_internal_datasource (FILE 64bins19window0error.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/64bins19window0error.out
                             URL_HASH SHA256=6a8a88f7cba0aae01f8cde75e875753c196704d36cf7430ae9c385ff7843b137)
declare_internal_datasource (FILE 64bins19window1error.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/64bins19window1error.out
                             URL_HASH SHA256=3b65eb1cda3a36b30412931f8eb6e167789d3f97ff45d52d40ccd692dfb6b922)
declare_internal_datasource (FILE 64bins23window1error.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/64bins23window1error.out
                             URL_HASH SHA256=3b65eb1cda3a36b30412931f8eb6e167789d3f97ff45d52d40ccd692dfb6b922)
declare_internal_datasource (FILE 128bins19window0error.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/128bins19window0error.out
                             URL_HASH SHA256=ed7832ba7a7e945e5501aa86520daf2b3dce013a3c67f515faa701f639acc740)
declare_internal_datasource (FILE 128bins19window1error.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/128bins19window1error.out
                             URL_HASH SHA256=a331b0c4ed4947ffa7fa44c1b99d50803d2b1e289365594286b273bd57fe45dd)
declare_internal_datasource (FILE 128bins23window1error.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/128bins23window1error.out
                             URL_HASH SHA256=a331b0c4ed4947ffa7fa44c1b99d50803d2b1e289365594286b273bd57fe45dd)

declare_internal_datasource (FILE 1bins19window0errorsocks.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/1bins19window0errorsocks.out
                             URL_HASH SHA256=b1053f5b749071c54528f7b00cc66f278eb50f5643d46db01774932ea510fbbd)
declare_internal_datasource (FILE 64bins19window0errorsocks.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/64bins19window0errorsocks.out
                             URL_HASH SHA256=d5e883ce66d1fc574cd6c9ae62c2d25a7671ad7be0948c49d969b1de88293def)
declare_internal_datasource (FILE 128bins19window0errorsocks.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/128bins19window0errorsocks.out
                             URL_HASH SHA256=3087b98c1f4e0d0a1a54f0dc11100f25fdfa8eba62d8f24c082ea1c523146618)

declare_internal_datasource (FILE 1bins19window0errorhibf.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/1bins19window0errorhibf.out
                             URL_HASH SHA256=1f7c7541afd7df07a8e3ef31d83f7e64ba63901fc47989a7275865aa418ad4c4)
declare_internal_datasource (FILE 1bins19window1errorhibf.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/1bins19window1errorhibf.out
                             URL_HASH SHA256=1f7c7541afd7df07a8e3ef31d83f7e64ba63901fc47989a7275865aa418ad4c4)
declare_internal_datasource (FILE 64bins19window0errorhibf.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/64bins19window0errorhibf.out
                             URL_HASH SHA256=d831fe3ce5e769ec89e955d1873ca97f60f910748b5f9dbba034c48d60503040)
declare_internal_datasource (FILE 64bins19window1errorhibf.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/64bins19window1errorhibf.out
                             URL_HASH SHA256=87fecb71ee97ebc5c8d8a3f6797417fe5bb156a9701537ee53000d6935e649f3)
declare_internal_datasource (FILE 128bins19window0errorhibf.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/128bins19window0errorhibf.out
                             URL_HASH SHA256=1143c77b32b24fddda41a9c4cbc0bf2ab4f6e3da9d7db203a35620bfb2419820)
declare_internal_datasource (FILE 128bins19window1errorhibf.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/128bins19window1errorhibf.out
                             URL_HASH SHA256=8d57920d39924df15d69e550cb12485cd0b53ca13252adcf9952916f97386f6a)

declare_internal_datasource (FILE 1binsempty.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/1binsempty.out
                             URL_HASH SHA256=515b16f2e507b9f4bac1b717c65d3a16b7c0f2bebe0977f7e0a471479dbcff42)
declare_internal_datasource (FILE 64binsempty.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/64binsempty.out
                             URL_HASH SHA256=3ca3257deccc01931390e3b0185d26ab97775a022eb8696930d680ba8701454e)
declare_internal_datasource (FILE 128binsempty.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/128binsempty.out
                             URL_HASH SHA256=e4eee8bce310b11bdb4da718fe3f3129508955b35cb89c7dc0f61b952d48410b)

declare_internal_datasource (FILE 1binsemptyhibf.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/1binsemptyhibf.out
                             URL_HASH SHA256=0ad9dd0d12a84e83cdf71343937b6d129fb3b44b9667b4b1345d3f14de9991af)
declare_internal_datasource (FILE 64binsemptyhibf.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/64binsemptyhibf.out
                             URL_HASH SHA256=e0d47758dff99736ae4aa5cd20cfc6435fed3848579c3d99697e40371f59c36f)
declare_internal_datasource (FILE 128binsemptyhibf.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/128binsemptyhibf.out
                             URL_HASH SHA256=fb811d81504d98670451753d81fb47bb09021679c3204c34d64ad752fdb5f9e9)

declare_internal_datasource (FILE three_levels.out
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/three_levels.out
                             URL_HASH SHA256=40da60ce362374330ffe803b3472e2d868f04bf27fdfaa9bb24fd247a67ab225)

declare_internal_datasource (FILE 1bins.pack
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/1bins.pack
                             URL_HASH SHA256=bd3475ee141be102213fbfb7c8c625bb80017e3dd65d7c13f540df2d2c857d43
                             CONFIGURE true)
declare_internal_datasource (FILE 64bins.pack
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/64bins.pack
                             URL_HASH SHA256=f2efc5e2ce4a1aa6235b3be96a67c6dcc9a540dc945ec6bd16574d913ee18ed6
                             CONFIGURE true)
declare_internal_datasource (FILE 128bins.pack
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/128bins.pack
                             URL_HASH SHA256=fcd607d2ed5a1d6a18e1bd332f48b56ee03c666a1d99d817f3f5b0b99490861c
                             CONFIGURE true)
declare_internal_datasource (FILE three_levels.pack
                             URL ${CMAKE_SOURCE_DIR}/test/data/hibf/three_levels.pack
                             URL_HASH SHA256=f4d0061995ba803265396e195154eca1629cdce331354f7a44f1ef1c0a075079
                             CONFIGURE true)

