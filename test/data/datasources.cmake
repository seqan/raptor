cmake_minimum_required (VERSION 3.18)

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
                             URL_HASH SHA256=0107f74c9fc616e116750186324ba25eb6881c33924149874559e8be740e9213
)
declare_internal_datasource (FILE 1bins23window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/1bins23window.index
                             URL_HASH SHA256=e41c59e98cfddb0c82759ffec8af1d2ccb689b603b58a19280c9d873473addb8
)
declare_internal_datasource (FILE 64bins19window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/64bins19window.index
                             URL_HASH SHA256=9da9bb2aa2f9c4211e8d1a78133ade54192aa221818a910fb66fd3673401e72a
)
declare_internal_datasource (FILE 64bins23window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/64bins23window.index
                             URL_HASH SHA256=f1c7047f80287912b23be6d94888f302a7ea40764849b2cfa661c05e72acb9ac
)
declare_internal_datasource (FILE 128bins19window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/128bins19window.index
                             URL_HASH SHA256=d4b6a15e6d64c35876c3f69cab3c34b013cfd4c81f50f8870f84704bab944d0e
)
declare_internal_datasource (FILE 128bins23window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/128bins23window.index
                             URL_HASH SHA256=227a504976e881c9b0f79e817c105b50f0c2b3e9b8a38acea590593f7d786973
)

declare_internal_datasource (FILE 1bins19window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/1bins19window.hibf
                             URL_HASH SHA256=e4576bf8d32e87c76d1dd8fbd6815c5c600d0dfc9fe181eb8cfd01274c2d88e2
)
declare_internal_datasource (FILE 1bins23window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/1bins23window.hibf
                             URL_HASH SHA256=11b519d88d6926d0fb00cee838b573f04636a0baccfd4c08377c2bccc7b4ec52
)
declare_internal_datasource (FILE 64bins19window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/64bins19window.hibf
                             URL_HASH SHA256=9fdff7c58b4918e9616f5f8eb9de128fd944d45dec81e297966d6e9a5d2ecbcc
)
declare_internal_datasource (FILE 64bins23window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/64bins23window.hibf
                             URL_HASH SHA256=03189c55bf09be0bbe03bb0953852277eaf3a8870be5cfe9bfbe5badf3ecffee
)
declare_internal_datasource (FILE 128bins19window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/128bins19window.hibf
                             URL_HASH SHA256=79f96e60efa557fc027ea98ab0f285e7d75a171b8d1f8c661f0b8e986e2c0b29
)
declare_internal_datasource (FILE 128bins23window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/128bins23window.hibf
                             URL_HASH SHA256=8b7125958b149a78cd97283788a8ccfe1bce7075ff8c522bf477e5da30acc3c8
)
declare_internal_datasource (FILE three_levels.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/three_levels.hibf
                             URL_HASH SHA256=5c80e67e99d1de1ec78701717d16d50c78976253dc3451b27645c4c2567d8d0a
)

declare_internal_datasource (FILE 1bins.pack
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/1bins.pack
                             URL_HASH SHA256=5130c56fd5ba7a392c335fbeea38e16fdef0d895215a2dfa0815d203dcc62158
                             CONFIGURE true
)
declare_internal_datasource (FILE 64bins.pack
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/64bins.pack
                             URL_HASH SHA256=37bd8f436fc0b5d805aa0e65baab8b9ab75c7b902fda2afe604aacf367676a55
                             CONFIGURE true
)
declare_internal_datasource (FILE 128bins.pack
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/128bins.pack
                             URL_HASH SHA256=425467230374a1a50bb464279ef883d3d1bd14df3f81449ba49ec3c9f87ff28a
                             CONFIGURE true
)
declare_internal_datasource (FILE three_levels.pack
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/three_levels.pack
                             URL_HASH SHA256=d0247914a3050cfb6e4c52f9b84ac50b7806d8d75c7586c96d4e184c66499494
                             CONFIGURE true
)
declare_internal_datasource (FILE test.layout
                             URL ${CMAKE_CURRENT_LIST_DIR}/test.layout
                             URL_HASH SHA256=0e9344550ca8de68fe51430b46100f39e8806e590bab7908479f87cd54f9586d
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
