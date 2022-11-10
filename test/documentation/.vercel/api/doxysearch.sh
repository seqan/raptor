import "querystring@1.3.0"

build() {
    # Every file generated inside the working directory during build() will be bundled.
    # We need the binary.
    cp ${WORK_PATH}/doxygen-bin/doxysearch.cgi .
    # We need the dynamically linked library.
    cp ${WORK_PATH}/doxygen-bin/libxapian.so libxapian.so.22
    # We need the index.
    cp -r ${WORK_PATH}/doxygen-bin/doxysearch.db .
}

# Extracts the querystring and calls doxysearch.
# Doxysearch returns a content-header that must be removed (tail -n +3).
handler() {
    local path="$(jq -r '.path' < "$1")"
    local string=$(querystring "$path" | querystring_unescape)
    export LD_LIBRARY_PATH="$(pwd):/lib64"
    ./doxysearch.cgi "$string" | tail -n +3
}
