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
declare_internal_datasource (FILE multi_record_bin.fa
                             URL ${CMAKE_CURRENT_LIST_DIR}/multi_record_bin.fa
                             URL_HASH SHA256=1647571eb3a9bde7212bed54ef5c323b9e02acbea40cbebcd8adee055ba20a6c
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
                             URL_HASH SHA256=6df6b7e912ab308cc18622a6a806aaeef8b2e00dfc8c1d7c24bb6ceac5b91767
)
declare_internal_datasource (FILE 1bins23window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/1bins23window.index
                             URL_HASH SHA256=d37b5918f96d169dc5cf230fed755db3c2e7770eded49902a26ab807801a245d
)
declare_internal_datasource (FILE 64bins19window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/64bins19window.index
                             URL_HASH SHA256=2e31467af179aea5da0fba49382e9a502fad7f973c30b6c6e2f8014bd266cb5d
)
declare_internal_datasource (FILE 64bins23window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/64bins23window.index
                             URL_HASH SHA256=05213d4ad9fd30e1f32f6aec6ee067242369189a39250a21f998491077741449
)
declare_internal_datasource (FILE 128bins19window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/128bins19window.index
                             URL_HASH SHA256=279fc8e3b4e4547b77ed9e97fbeb2e8578e6127673fa349137a0bde29ac95e0f
)
declare_internal_datasource (FILE 128bins23window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/128bins23window.index
                             URL_HASH SHA256=7b821510a57b73eef9686fdc3108180d15d1f5e92a6078b156179877c64ccd5a
)

declare_internal_datasource (FILE 1bins19window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/1bins19window.hibf
                             URL_HASH SHA256=16551c3450b561c0441561cb1c88641f6c8ab4ec293aba1e92196a77decbd028
)
declare_internal_datasource (FILE 1bins23window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/1bins23window.hibf
                             URL_HASH SHA256=6bcf45079bb2ad014eaba31b18a1ad9ef353e51cdb93b0089b401ac99191ad91
)
declare_internal_datasource (FILE 64bins19window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/64bins19window.hibf
                             URL_HASH SHA256=f4db9470d5b3f08e0b26635b014cb96fe61b75ec21195047626a971d3cc52cda
)
declare_internal_datasource (FILE 64bins23window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/64bins23window.hibf
                             URL_HASH SHA256=98ece2a4bbbc14cab5bb867f420701b0cf68701c42436e4d9cb3c10f5e01494f
)
declare_internal_datasource (FILE 128bins19window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/128bins19window.hibf
                             URL_HASH SHA256=11c50229a0740ebf376f12a0e83183bbdecfb61298c5377ac0045d586ef1cba4
)
declare_internal_datasource (FILE 128bins23window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/128bins23window.hibf
                             URL_HASH SHA256=4b37222cd8fb10c787a4c6f7e753739a074df264df763e6d5d368963c6f36514
)
declare_internal_datasource (FILE three_levels.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/three_levels.hibf
                             URL_HASH SHA256=ca2f730d6f60a4e6f4eaf91e6b2f30c4a2aaa09f7b489f44b3e440adf5f5d1dd
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
