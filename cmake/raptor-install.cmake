# --------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------------

# This file describes where and which parts of Sharg should be installed to.

cmake_minimum_required (VERSION 3.14)

include (GNUInstallDirs)

# install documentation files in /share/doc
install (FILES "${RAPTOR_CLONE_DIR}/CHANGELOG.md" #
               "${RAPTOR_CLONE_DIR}/LICENSE.md" #
               "${RAPTOR_CLONE_DIR}/README.md" TYPE DOC
)

# install cmake files in /share/cmake
install (FILES "${RAPTOR_CLONE_DIR}/cmake/raptor-config.cmake" "${RAPTOR_CLONE_DIR}/cmake/raptor-config-version.cmake"
         DESTINATION "${CMAKE_INSTALL_DATADIR}/cmake/raptor"
)

# install raptor header files in /include/raptor
install (DIRECTORY "${RAPTOR_INCLUDE_DIR}/raptor" TYPE INCLUDE)

# install submodule header files, e.g. all external dependencies in /home/user/seqan3/submodules/*,
# in /include/seqan3/submodules/<submodule>/include
foreach (submodule_dir ${RAPTOR_DEPENDENCY_INCLUDE_DIRS})
    # e.g. submodule_dir: (1) /home/user/seqan3/submodules/sdsl-lite/include or (2) /usr/include
    # strip /home/user/seqan3/submodules/ and /include part.
    file (RELATIVE_PATH submodule "${RAPTOR_SUBMODULES_DIR}/lib" "${submodule_dir}/..")
    # submodule is either a single module name, like sdsl-lite or a relative path to a folder ../../../usr
    # skip relative folders and only keep submodules that reside in the submodules folder
    if (NOT submodule MATCHES "^\\.\\.") # skip relative folders
        install (DIRECTORY "${submodule_dir}" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/raptor/lib/${submodule}")
    endif ()
endforeach ()
