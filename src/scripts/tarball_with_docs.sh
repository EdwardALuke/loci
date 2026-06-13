#!/bin/bash
###############################################################################
#
# Copyright 2008-2025, Mississippi State University
#
# This file is part of the Loci Framework.
#
# The Loci Framework is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The Loci Framework is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
#
###############################################################################

set -e
set -u

usage() {
    echo "usage: $0 <revision-name> <git-info> <git-branch> <obj-dir>"
}

if [ "$#" -ne 4 ]; then
    usage
    exit 2
fi

revision_name=$1
git_info=$2
git_branch=$3
obj_dir=$4

# LOCI_REVISION_NAME is defined with single quotes in src/conf/version.conf so
# make can safely use branch names in shell commands.  The historical tarball
# recipe lets the shell consume those quotes; mirror that filename convention.
revision_name=${revision_name#\'}
revision_name=${revision_name%\'}

source_root=$(pwd)
if [ "${obj_dir:0:1}" = "/" ]; then
    obj_path=$obj_dir
else
    obj_path=$source_root/$obj_dir
fi

prefix=Loci-$revision_name
archive=$prefix-with-docs.tgz
tutorial_pdf=$obj_path/docs/tutorial/docs/tutorial.pdf
doxygen_html=$obj_path/docs/doxygen/html

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "ERROR: tarball-with-docs must be run from a git checkout."
    exit 1
fi

if [ ! -f "$tutorial_pdf" ]; then
    echo "ERROR: $tutorial_pdf not found."
    echo "Run 'make docs' with pdflatex available before packaging docs."
    exit 1
fi

if [ ! -f "$doxygen_html/index.html" ]; then
    echo "ERROR: $doxygen_html/index.html not found."
    echo "Run 'make docs' with doxygen available before packaging docs."
    exit 1
fi

stage_dir=$(mktemp -d "${TMPDIR:-/tmp}/loci-tarball-with-docs.XXXXXX")
trap 'rm -rf "$stage_dir"' EXIT

version_file=$stage_dir/$prefix/src/version.conf
manifest=$stage_dir/$prefix/docs/generated-docs-manifest.txt

git archive --format=tar --prefix="$prefix/" HEAD | tar -x -C "$stage_dir"

{
    echo "GIT_INFO   = $git_info"
    echo "GIT_BRANCH = $git_branch"
    sed -n '3,$p' src/conf/version.conf
} > "$version_file"

mkdir -p "$stage_dir/$prefix/docs/tutorial/docs"
cp -p "$tutorial_pdf" "$stage_dir/$prefix/docs/tutorial/docs/tutorial.pdf"

rm -rf "$stage_dir/$prefix/docs/doxygen/html"
mkdir -p "$stage_dir/$prefix/docs/doxygen"
cp -aL "$doxygen_html" "$stage_dir/$prefix/docs/doxygen/html"

if command -v doxygen >/dev/null 2>&1; then
    doxygen_version=$(doxygen --version)
else
    doxygen_version="not found"
fi

if command -v pdflatex >/dev/null 2>&1; then
    pdflatex_version=$(pdflatex --version | sed -n '1p')
else
    pdflatex_version="not found"
fi

{
    echo "Generated documentation manifest"
    echo
    echo "Archive: $archive"
    echo "Source commit: $(git rev-parse HEAD)"
    echo "Git info: $git_info"
    echo "Git branch: $git_branch"
    echo "Generated at UTC: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "Doxygen: $doxygen_version"
    echo "pdflatex: $pdflatex_version"
    echo
    echo "Included generated artifacts:"
    echo "- docs/tutorial/docs/tutorial.pdf"
    echo "- docs/doxygen/html/"
} > "$manifest"

tar -czf "$archive" -C "$stage_dir" "$prefix"
echo "Created $archive"
