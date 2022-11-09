import "querystring@1.3.0"

build() {
    cp ${WORK_PATH}/doxygen-bin/doxysearch.cgi .
    cp ${WORK_PATH}/doxygen-bin/libxapian.so libxapian.so.22
    cp -r ${WORK_PATH}/doxygen-bin/doxysearch.db .
}

handler() {
    local path
    path="$(jq -r '.path' < "$1")"
    string=$(querystring "$path" | querystring_unescape)
    export LD_LIBRARY_PATH="$(pwd):/lib64"
    ./doxysearch.cgi "$string" | tail -n +3
}
