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

declare_internal_datasource (FILE 1bins19window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/1bins19window.index
                             URL_HASH SHA256=ddf786e06302fba9f24b0725b00f039c606b59993894f05ccea1cb4ecd1645a5
)
declare_internal_datasource (FILE 1bins23window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/1bins23window.index
                             URL_HASH SHA256=da2b106b3f92c4006b1664206200a48909d22ea0012122273f5c1ee58e317150
)
declare_internal_datasource (FILE 64bins19window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/64bins19window.index
                             URL_HASH SHA256=faf3f7ac2bb0651313247231a57b7835ff29371fc49546dac18ae38489472a1b
)
declare_internal_datasource (FILE 64bins23window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/64bins23window.index
                             URL_HASH SHA256=28bf2a0e465a197dd6f4c3bbe84e492f589870b6425bcdbad6d060171c7002e0
)
declare_internal_datasource (FILE 128bins19window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/128bins19window.index
                             URL_HASH SHA256=a5fee0afc79f64a4a9c7f0d257b5f076f84753824892bd1ddca0d9abd577b879
)
declare_internal_datasource (FILE 128bins23window.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/128bins23window.index
                             URL_HASH SHA256=500afa4fc9ac0acca90c108b4e9ad664bad1717c92a4bb3d39ce4e1869f684e0
)

declare_internal_datasource (FILE 2.0.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/2.0.index
                             URL_HASH SHA256=c7dcacfa9b170446a279c33e13663c8f45c07afe657d9fb127b6a797afd7767d
)
declare_internal_datasource (FILE 2.0.compressed.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/2.0.compressed.index
                             URL_HASH SHA256=220fa3cf85a3669ffc67eb9664918092367d53f2d76db821db8048e9bce40446
)
declare_internal_datasource (FILE 2.0.partitioned.index_0
                             URL ${CMAKE_CURRENT_LIST_DIR}/2.0.partitioned.index_0
                             URL_HASH SHA256=54b7e93b39c4b67df704b639d0f45a1365e809188a899c6bd11ce34cc9cd77fd
)
declare_internal_datasource (FILE 2.0.partitioned.index_1
                             URL ${CMAKE_CURRENT_LIST_DIR}/2.0.partitioned.index_1
                             URL_HASH SHA256=732c155bf2ed32be2da3d00e5907b58cf7b8adf92f327cf4d27ec160b8456a8b
)
declare_internal_datasource (FILE 2.0.partitioned.index_2
                             URL ${CMAKE_CURRENT_LIST_DIR}/2.0.partitioned.index_2
                             URL_HASH SHA256=27f2d356f54baf061ca6f89022cf623921c4a2414f3eebea5efd65bc943c5070
)
declare_internal_datasource (FILE 2.0.partitioned.index_3
                             URL ${CMAKE_CURRENT_LIST_DIR}/2.0.partitioned.index_3
                             URL_HASH SHA256=6e8358f1757c9667541a14e3ac637fbc9efcbfcf520daa28db2f27e4e4d6ee9d
)

declare_internal_datasource (FILE 64bins23windowc.index
                             URL ${CMAKE_CURRENT_LIST_DIR}/64bins23windowc.index
                             URL_HASH SHA256=a52aa73bb64f98341b9473b651149d130cd2c1f9a889ffd8a7d65a1c134d4eb7
)

declare_internal_datasource (FILE 1bins19window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/1bins19window.hibf
                             URL_HASH SHA256=2a49572c03edded5bd80f6a21d360344f5f39667c223d8a7652295a72673b973
)
declare_internal_datasource (FILE 1bins23window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/1bins23window.hibf
                             URL_HASH SHA256=14eb9a5b77ab98deaf0ec2c126107fce4bf558815acb4ab509d471ee129cf902
)
declare_internal_datasource (FILE 64bins19window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/64bins19window.hibf
                             URL_HASH SHA256=c3a292d9883055c1dea202e43d8a63e3a279e46f17cf9a8e78b3f6f228b9d920
)
declare_internal_datasource (FILE 64bins23window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/64bins23window.hibf
                             URL_HASH SHA256=c11b4381f8bdf086ffd1fb8d37ce361c12adb3475a8eedf103b09e20d4f2401d
)
declare_internal_datasource (FILE 128bins19window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/128bins19window.hibf
                             URL_HASH SHA256=21165f0bef2ac8c3b2a9c243d3afb60febdf18507b5d6154bdb3a074339caa96
)
declare_internal_datasource (FILE 128bins23window.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/128bins23window.hibf
                             URL_HASH SHA256=da8eecf7a1ed4e5aa09c3de09088dacbb6d8549cbc598c4b047df7894080416e
)
declare_internal_datasource (FILE three_levels.hibf
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/three_levels.hibf
                             URL_HASH SHA256=d04b045cf1db05a5ebc53b26e1bd2aaec430742310698dd4605990ab962fa439
)

declare_internal_datasource (FILE 1bins.pack
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/1bins.pack
                             URL_HASH SHA256=bd3475ee141be102213fbfb7c8c625bb80017e3dd65d7c13f540df2d2c857d43
                             CONFIGURE true
)
declare_internal_datasource (FILE 64bins.pack
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/64bins.pack
                             URL_HASH SHA256=f2efc5e2ce4a1aa6235b3be96a67c6dcc9a540dc945ec6bd16574d913ee18ed6
                             CONFIGURE true
)
declare_internal_datasource (FILE 128bins.pack
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/128bins.pack
                             URL_HASH SHA256=fcd607d2ed5a1d6a18e1bd332f48b56ee03c666a1d99d817f3f5b0b99490861c
                             CONFIGURE true
)
declare_internal_datasource (FILE three_levels.pack
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/three_levels.pack
                             URL_HASH SHA256=f4d0061995ba803265396e195154eca1629cdce331354f7a44f1ef1c0a075079
                             CONFIGURE true
)
declare_internal_datasource (FILE 1bins.pack
                             URL ${CMAKE_CURRENT_LIST_DIR}/hibf/1bins.pack
                             URL_HASH SHA256=bd3475ee141be102213fbfb7c8c625bb80017e3dd65d7c13f540df2d2c857d43
                             CONFIGURE true
)
declare_internal_datasource (FILE test.layout
                             URL ${CMAKE_CURRENT_LIST_DIR}/test.layout
                             URL_HASH SHA256=28b1f1e62814aa3fa389704111b6c50eac04b97d4ebbde920fc1db7f84b93aef
                             CONFIGURE true
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

# for bins in 1 64 128; do for window in 19 23; do sha256sum /home/infri/develop/raptor/build/unit/GCC-12.2.0-Debug/output/build_hibf_suite/build_hibf.with_file/${bins}_bins_${window}_window_serial/raptor.index; done; done
# for bins in 1 64 128; do for window in 19 23; do cp /home/infri/develop/raptor/build/unit/GCC-12.2.0-Debug/output/build_hibf_suite/build_hibf.with_file/${bins}_bins_${window}_window_serial/raptor.index /home/infri/develop/raptor/test/data/hibf/${bins}bins${window}window.hibf; done; done
