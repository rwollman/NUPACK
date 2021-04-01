vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO mfornace/openmm-1
    REF ef016e56eb3a2945abfc9ced76b9e2c9a9c0ad96
    SHA512 02d5d5a6e4629a303e382ea14ca16199d0ca462853f7781b68349d5541ab3d0a00ce65eea2b0cc8f47afc476fbef323599563d668ad343fd1b035f19fb571d54
    HEAD_REF master
)

vcpkg_configure_cmake(
    SOURCE_PATH ${SOURCE_PATH}
    PREFER_NINJA
    OPTIONS
        -DBUILD_TESTING=OFF
)

vcpkg_install_cmake()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/docs")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/docs")

# Handle copyright
file(INSTALL ${SOURCE_PATH}/docs-source/licenses/Licenses.txt DESTINATION ${CURRENT_PACKAGES_DIR}/share/${PORT} RENAME copyright)
