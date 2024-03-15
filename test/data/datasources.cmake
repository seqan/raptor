cmake_minimum_required (VERSION 3.25)

declare_internal_datasource (FILE bin1.fa
                             URL ${CMAKE_CURRENT_LIST_DIR}/bin1.fa
                             URL_HASH SHA256=03ac546c2dfa34f2a75283d43dd2293aa702ed2b7dd877dfbb4f966d258201f0
)
declare_internal_datasource (FILE bin2.fa
                             URL ${CMAKE_CURRENT_LIST_DIR}/bin2.fa
                             URL_HASH SHA256=5bbf2f6fee427d3e967932203987b8c391ac86f0f089101455d529cb7e27b985
)
declare_internal_datasource (FILE bin3.fa
                             URL ${CMAKE_CURRENT_LIST_DIR}/bin3.fa
                             URL_HASH SHA256=93d80997611587f9635bdf7ed48eff26d995b54899e757026551df90e37464e7
)
declare_internal_datasource (FILE bin4.fa
                             URL ${CMAKE_CURRENT_LIST_DIR}/bin4.fa
                             URL_HASH SHA256=b4e3be3c001aac8481be2251aab110e4318996d26f2cdc91fb3061c35cbec4a1
)

declare_internal_datasource (FILE bin1.fa.gz
                             URL ${CMAKE_CURRENT_LIST_DIR}/bin1.fa.gz
                             URL_HASH SHA256=8c18f83a6e5c3fcde122052ceb843beadc7eaf31157b68b7a7caa70004be86ab
)
declare_internal_datasource (FILE bin2.fa.gz
                             URL ${CMAKE_CURRENT_LIST_DIR}/bin2.fa.gz
                             URL_HASH SHA256=60534c9bb4b593d39cedd5039614eb7eec6ab92961c490bb461232d9626d71a9
)
declare_internal_datasource (FILE bin3.fa.gz
                             URL ${CMAKE_CURRENT_LIST_DIR}/bin3.fa.gz
                             URL_HASH SHA256=1139affca816f0b404011a038d8c9d94481e28f84e566a57da4be4abc0b9eb13
)
declare_internal_datasource (FILE bin4.fa.gz
                             URL ${CMAKE_CURRENT_LIST_DIR}/bin4.fa.gz
                             URL_HASH SHA256=9986c8665294cc4a5eaf1341990f8680afd2ffa1aa6b4906d269f7ca96c393f2
)

declare_internal_datasource (FILE query.fq
                             URL ${CMAKE_CURRENT_LIST_DIR}/query.fq
                             URL_HASH SHA256=f48eb3f357e23df89e7e15d2f77f9285a428f73c4903eb1c6580271e0dea3d87
)
declare_internal_datasource (FILE query_empty.fq
                             URL ${CMAKE_CURRENT_LIST_DIR}/query_empty.fq
                             URL_HASH SHA256=80cd7628b6fdcb7dbbe99dd7f05a7a5578b462db9dd67b6d0e6982b03c5d4cd5
)

declare_internal_datasource (FILE 3.0.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/3.0.index
                             URL_HASH SHA256=59850d63d6fa6c9ac11bdfe8c9246aafa026b8cb4769300cee1f9f40bca46920
)
declare_internal_datasource (FILE 3.0.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/3.0.hibf
                             URL_HASH SHA256=c3a292d9883055c1dea202e43d8a63e3a279e46f17cf9a8e78b3f6f228b9d920
)
declare_internal_datasource (FILE 3.0.partitioned.index_0
                             URL ${CMAKE_CURRENT_LIST_DIR}/3.0.partitioned.index_0
                             URL_HASH SHA256=dd268cf4d111dd32abcea5765bc5a1c5db7a089f1e8354694a60ff15f0f1938d
)
declare_internal_datasource (FILE 3.0.partitioned.index_1
                             URL ${CMAKE_CURRENT_LIST_DIR}/3.0.partitioned.index_1
                             URL_HASH SHA256=c09b3d5d1119f3d220a2994abd903c94c87214191fe6e0f4e41a3c98d7122271
)
declare_internal_datasource (FILE 3.0.partitioned.index_2
                             URL ${CMAKE_CURRENT_LIST_DIR}/3.0.partitioned.index_2
                             URL_HASH SHA256=4c2c2f982bafb7c026f6e0af60702b5c62a2724820b9bc88f3733cce3a700e5b
)
declare_internal_datasource (FILE 3.0.partitioned.index_3
                             URL ${CMAKE_CURRENT_LIST_DIR}/3.0.partitioned.index_3
                             URL_HASH SHA256=7f71057f9257af688da757b60335e25895e31c03d21dff8d8f5844a041d5c2b8
)

declare_internal_datasource (FILE 1bins19window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/1bins19window.index
                             URL_HASH SHA256=44b072342ab71eaa50c87a70802a3289ed861744bf7b93160df6eb82db380058
)
declare_internal_datasource (FILE 1bins23window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/1bins23window.index
                             URL_HASH SHA256=2a24ad4a8c01b2cd16a5dbe7d94df19170face1b0ad59931c37404f1d0fdffda
)
declare_internal_datasource (FILE 64bins19window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/64bins19window.index
                             URL_HASH SHA256=1499735d434ce0aeb27649e80e515814896c270de5b3a3914d82db9bf89a80f8
)
declare_internal_datasource (FILE 64bins23window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/64bins23window.index
                             URL_HASH SHA256=41f1ac29d3fb72f10fdde17aa433a651722f3c703806a2f02dd267d424f80460
)
declare_internal_datasource (FILE 128bins19window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/128bins19window.index
                             URL_HASH SHA256=adbb3a9d681f8e8b734eb9d732ab6600a743a7a22650d163eb5f5faf0f02bf33
)
declare_internal_datasource (FILE 128bins23window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/128bins23window.index
                             URL_HASH SHA256=dac3c7bd384b14a2bfcf8cbee0468e87b64409f5a52588c97c9d38d54063f9fb
)

declare_internal_datasource (FILE 1bins19window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/1bins19window.hibf
                             URL_HASH SHA256=c62b862e21516cce4fed683cac8fb223e50c753315e7555eef13f67733b69b66
)
declare_internal_datasource (FILE 1bins23window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/1bins23window.hibf
                             URL_HASH SHA256=bcab2c0d30870734cd16bc6b1e184967a6c517cac7bca09777c398456c6cb65f
)
declare_internal_datasource (FILE 64bins19window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/64bins19window.hibf
                             URL_HASH SHA256=089561152351e2a6e235d3b41d2be345edfc5172c46599747477244d33fccabd
)
declare_internal_datasource (FILE 64bins23window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/64bins23window.hibf
                             URL_HASH SHA256=f821b875f9853973b8323fad08ea4a07fce3d8aa048f0ebc42e51e95cfd5ebfb
)
declare_internal_datasource (FILE 128bins19window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/128bins19window.hibf
                             URL_HASH SHA256=19ee4c7a5a32e73d329fc9112e43692a8db2749eae23a0601932411041c0d1a4
)
declare_internal_datasource (FILE 128bins23window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/128bins23window.hibf
                             URL_HASH SHA256=a581085d7bf1a5caf73f06518a4591b8b5f0d1c7a0206f8bae68facef0835c63
)
declare_internal_datasource (FILE three_levels.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/three_levels.hibf
                             URL_HASH SHA256=5d0f2454affc05861c42cce7de4021ee4aa59938d1f9e425205b539417b0e5c1
)

declare_internal_datasource (FILE 1bins.pack
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/1bins.pack
                             URL_HASH SHA256=b89a1d3f71760edd7013e4c2336d2725ecc891ff1be50a9e8800ae2a1ba8b295
                             CONFIGURE true
)
declare_internal_datasource (FILE 64bins.pack
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/64bins.pack
                             URL_HASH SHA256=9a4ff8fd71fa9bb733abd8bf66aee07aea05fb12e7ac9d806d80af0f7a3e66fc
                             CONFIGURE true
)
declare_internal_datasource (FILE 128bins.pack
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/128bins.pack
                             URL_HASH SHA256=c0778bade7ff1b891aa842a3a61264d911e12a190247d6436adb27ab570dc1ae
                             CONFIGURE true
)
declare_internal_datasource (FILE three_levels.pack
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/three_levels.pack
                             URL_HASH SHA256=6ed58e2e04a93085168a524226b6d9c4305643c8aa4ec925ccfa16cc158a61e9
                             CONFIGURE true
)
declare_internal_datasource (FILE test.layout
                             URL ${CMAKE_CURRENT_LIST_DIR}/test.layout
                             URL_HASH SHA256=db98c48b482378e4d36131af1c73f8328df180096bf11e1b8bcf025cdb5d643a
                             CONFIGURE true
)
declare_internal_datasource (FILE test_nocfg.layout
                             URL ${CMAKE_CURRENT_LIST_DIR}/test_nocfg.layout
                             URL_HASH SHA256=df071bb98ab4dd06d052fc105e1cc65cb34b581dcc6240ca2fb85b39a0181306
                             CONFIGURE true
)

declare_internal_datasource (FILE bin1.header
                             URL ${CMAKE_CURRENT_LIST_DIR}/minimiser/bin1.header
                             URL_HASH SHA256=c1944ce2230abef4648895864c1f48f2072ae01234e49df74a66dc9a55e0b344
)
declare_internal_datasource (FILE bin1.minimiser
                             URL ${CMAKE_CURRENT_LIST_DIR}/minimiser/bin1.minimiser
                             URL_HASH SHA256=40c06f410d9723a9601444b643ec7124d7152a5c30c33a37d75b3d71d7cff00d
)

declare_internal_datasource (FILE empty.fq
                             URL ${CMAKE_CURRENT_LIST_DIR}/empty.fq
                             URL_HASH SHA256=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
)
declare_internal_datasource (FILE too_short.fq
                             URL ${CMAKE_CURRENT_LIST_DIR}/too_short.fq
                             URL_HASH SHA256=d6e5bf240c450ae800d188cc486763345c209589b178631195e20b0831d6ab3f
)
declare_internal_datasource (FILE query_variance.fq
                             URL ${CMAKE_CURRENT_LIST_DIR}/query_variance.fq
                             URL_HASH SHA256=6db9fd4f1b04fe62e14d03ff4a47aff74cbe4e59cce859b0e83309b14971e819
)

# Helper bash commands for bulk updates

# for bins in 1 64 128; do
#     for window in 19 23; do
#         cp <BUILD_DIR>/output/build_hibf_suite/build_hibf.with_file/${bins}_bins_${window}_window_serial/raptor.index <REPO_DIR>/raptor/test/data/hibf/${bins}bins${window}window.hibf
#     done
# done

# for bins in 1 64 128; do
#     for window in 19 23; do
#         sha256sum <REPO_DIR>/raptor/test/data/hibf/${bins}bins${window}window.hibf
#     done
# done

# cp <BUILD_DIR>/output/build_hibf.three_levels/raptor.index <REPO_DIR>/raptor/test/data/hibf/three_levels.hibf
# sha256sum <REPO_DIR>/raptor/test/data/hibf/three_levels.hibf

# for bins in 1 64 128; do
#     for window in 19 23; do
#         cp <BUILD_DIR>/output/build_ibf_suite/build_ibf.with_file/${bins}_bins_${window}_window_serial/raptor.index <REPO_DIR>/raptor/test/data/${bins}bins${window}window.index
#     done
# done

# for bins in 1 64 128; do
#     for window in 19 23; do
#         sha256sum <REPO_DIR>/raptor/test/data/${bins}bins${window}window.index
#     done
# done
