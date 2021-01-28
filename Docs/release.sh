#!/bin/bash

# Script for preparing a release of SIESTA.
# Script-author:
#  Nick R. Papior, 2016
#
# We encourage the use such as this script to ease and unify
# the release mechanism for SIESTA.
#
# To create the 4.0 release, this command should be used:
#   ./release.sh --prev-tag v3.2 --tag v4.0 \
#      --out siesta-4.0
# which creates the file:
#   ../siesta-releases/siesta-4.0.tar.gz
#
# We encourage that the prev-tag is ALWAYS the previous
# stable release tag (also for beta releases!).
# This ensures a consistent CHANGES file shipped with
# the tar-file.
# To clarify the different cases here we show a few
# release steps.
#
# Release a new version of siesta.
#
# In this case the previous stable release was: 4.0
# The new release is: 4.1
# The command would be:
#   ./release.sh --prev-tag v4.0 --tag v4.1 \
#      --out siesta-4.1
#
# Release a second beta release of siesta.
#
# In this case the previous stable release was: 4.0
# The new release is: 4.1-b2
# The command would be:
#   ./release.sh --prev-tag v4.0 --tag v4.1-b2 \
#      --out siesta-4.1-b2
# A subsequent beta release would be:
#   ./release.sh --prev-tag v4.0 --tag v4.1-b3 \
#      --out siesta-4.1-b3

echo ""
echo "Welcome to the SIESTA release script"
echo ""

# Create wrapper for pushd/popd to remove output
function pushd {
    command pushd "$@" > /dev/null
}
function popd {
    command popd "$@" > /dev/null
}


# This may actually not be changed, but is useful
# throughout.
main_dir=$(git rev-parse --show-toplevel)
# Ensure we are at the root of this file
pushd $main_dir


# First we setup the default options
# _reldir is the main directory of the SIESTA shared
# repository (or at least it should be)
_reldir=$(dirname $main_dir)/siesta-releases
_sign=1
_head=0
_tag=
# Get default output file (siesta-<>.tar.gz)
_out=

function _git_tag_cmd {
    git tag --list 'v*' --sort=-v:refname
}

# Default to latest tag
_def_tag=$(_git_tag_cmd | head -1)
_tag=$_def_tag

function has_tag {
    local tag=$1
    shift
    local ret=1
    for T in $(_git_tag_cmd)
    do
	if [ "$tag" == "$T" ]; then
	    ret=0
	fi
    done
    printf "%d" "$ret"
}

function help {
    local ret=$1
    shift

    echo "$0 creates a release of the SIESTA code at the tip-tag"
    echo ""
    echo "The following options may be used to control the release archive."
    echo ""
    echo "  --tag      instead of selecting the latest tag, choose this tag as the"
    echo "             reference tag for creating a release archive [$_def_tag]"
    echo "  --head     creates a release without signing, and sets tag to HEAD and out to"
    echo "             a unique identifier [git describe --abbrev], useful for making release for"
    echo "             sample groups or individuals without a proper release"
    echo "  --out      the default output file is siesta-<tag>.tar.gz. Do *not* specify tar.gz, "
    echo "  --no-sign  do not sign the output files (useful for test-runs)"
    echo "  --help|-h  show this help."
    exit $ret
}


# Function for signing a specific file
function sign {
    local file=$1
    shift
    if [ -e $file ]; then
	if [ $_sign -eq 1 ]; then
	    echo "Sign file: gpg --armor --sign --detach-sig $file"
	    gpg --armor --sign --detach-sig $file
	    store $file.asc
	    rm -f $file.asc
	fi
    fi
}

# Calculate, md5, sha256, sha512 and store it in the sha256.txt file
_checksum=checksums.txt
function checksums {
    local file=$1
    shift
    if [ -e $file ]; then
	# We move to the directory to ensure the
	# file names are as short as possible.
	openssl md5 $file >> $_reldir/$_checksum
	openssl sha256 $file >> $_reldir/$_checksum
	openssl sha512 $file >> $_reldir/$_checksum
    fi
}
    


# Save file in the $_tag-files folder
function store {
    local file=$1
    shift
    if [ -e $file ]; then
	cp $file $_reldir/$_tag-files/
    fi
}


# First we check whether the release-manager has
# specified the tag that is going to be shipped.
while [ $# -gt 0 ]; do
    opt=$1
    shift
    case $opt in
	--tag)
	    _tag=$1
	    # Check that the tag exists
	    if [ $(has_tag $_tag) -eq 1 ]; then
		echo ""
		echo "$0 could not find tag: '$_tag' in the tags list."
		echo "A release tag *MUST* exist before a release can be created"
		echo ""
		echo "The available tags are:"
		_git_tag_cmd | sed 's/^\(.*\)/  \1/g'
		exit 1
	    fi
	    shift
	    ;;

	--head)
	    _head=1
	    ;;

	--no-sign)
	    _sign=0
	    ;;
	
	--out)
	    _out=${1//.tar.gz/}
	    shift
	    ;;
	
	--help|-h)
	    help 0
	    ;;
	
    esac
done

if [ $_head -eq 1 ]; then
    _tag=$(git describe --abbrev)
    _sign=0
    _out=siesta-$_tag
fi

_tag_no_v=${_tag//v/}


# Get default output file (siesta-<>.tar.gz)
if [ -z "$_out" ]; then
    _out=siesta-${_tag_no_v}
    _out=${_out//-release/}
    _out=${_out//-rel/}
    _out=${_out//release-/}
    _out=${_out//rel-/}
fi

# Extract the release date for the tag
# Have to do this while in a git repository
_date=$(date -d "$(git log -n1 --format='%ci' $_tag)" +"%B %d, %Y")

echo "Chosen release tag is: $_tag"
echo "Creating out file: $_out.tar.gz"
echo "Release date: $_date"
echo ""
echo "Waiting 1 second before creating release... (Ctrl^C kills the sequence)"
sleep 1

# The current procedure of releasing a SIESTA
# version is the following:
#  0. Create the necessary directories
#  1. Create CHANGES and sign-store them.
#  2. Make the proper documentation
#  3. Go out of the top-directory to create the repository.

_check_dir_fail=0
function check_dir {
    local d=$1
    shift
    if [ -d $d ]; then
	echo ""
	echo "Please delete this folder before re-running the script:"
	echo "  rm -rf $d"
	_check_dir_fail=1
    fi
}

# Top release directories, in case it is not already
# created.
mkdir -p $_reldir

# Check directories
check_dir $_reldir/$_tag-files
check_dir $_reldir/$_out
if [ $_check_dir_fail -ne 0 ]; then
    exit 1
fi
# ensure the checksum file is "clean" on entry.
rm -f $_reldir/$_checksum

# The final directory with ALL release files, possibly signed.
mkdir -p $_reldir/$_tag-files


# Create a temporary work-directory
git archive -o $_reldir/$_out.tar.gz --prefix=$_out/ $_tag $main_dir
pushd $_reldir
tar xfz $_out.tar.gz ; rm $_out.tar.gz
popd

# Go into the release directory where all work will be done
pushd $_reldir

# Go into the source directory
pushd $_out

# Update the version.info file
printf "%s" "${_tag_no_v}" > version.info

# Create documentation
pushd Docs

# Update manual information that is version/date dependent
sed -s -i -e "s/\\date{.*}/\\date{$_date}/" siesta.tex tbtrans.tex
# Version tags in the pdf-title
sed -s -i -e "s/\\providecommand\\softwareversion{.*}/\\providecommand\\softwareversion{$_tag_no_v}/" siesta.tex tbtrans.tex

# First create the screen variants...
make final-screen
# Then the regular documentation (for print)
make final

# Clean-up the non-pdf files (currently
# this is not available through the make clean)
make clean
# Also do not ship the release script
rm release.sh
# Remove all .git related files
rm -rf .git*

# Create signatures and move files
for f in *.pdf ; do
    sign $f
    store $f
    checksums $f
done

# Go out of the documentation directory...
popd



# Return from the branch release
popd

# Tar the file
# In some future cases will the tag and out name
# be equivalent and we shouldn't delete what we should process!

# Create the archive... (without any .git files)
tar -czf $_out.tar.gz $_out
sign $_out.tar.gz
store $_out.tar.gz
checksums $_out.tar.gz
rm $_out.tar.gz

# Be sure to store the checksums.txt file
sign $_checksum
store $_checksum
rm $_checksum

popd
