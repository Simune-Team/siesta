#!/bin/bash

# Script for preparing a release
# of SIESTA.
# Script-author:
#  Nick R. Papior, 2016
#
# We encourage the use such as this:
# To create the 4.0 release, we would advocate performing
# this command:
#   ./release.sh --prev-tag siesta-3.2 --tag 4.0-release \
#      --out siesta-4.0
# which creates the file:
#   siesta-releases/siesta-4.0.tar.gz
# We encourage that the prev-tag is ALWAYS the previous
# stable release tag.
# This ensures a consistent changes file shipped with
# the tar-file.
# I.e. for the 4.1 release, the command would be:
#   ./release.sh --prev-tag 4.0-release --tag 4.1-release \
#      --out siesta-4.1


# This may actually not be changed, but is useful
# throughout.
main_dir=`bzr root`
pushd $main_dir

# First we setup the default options
_reldir=siesta-releases

# Get default tags...
_prev_tag=`bzr tags --sort=time | tail -2 | head -1 | awk '{print $1}'`
_tag=`bzr tags --sort=time | tail -1 | awk '{print $1}'`
# Get default output file (siesta-<>.tar.gz)
_out=siesta-${_tag//-release/}
_out=${_out//-rel/}
_out=${_out//release-/}
_out=${_out//rel-/}


function _has_tag {
    local tag=$1
    shift
    local ret=1
    for T in `bzr tags | awk '{print $1}'`
    do
	if [ "$tag" == "$T" ]; then
	    ret=0
	fi
    done
    printf "%d" "$ret"
}

function _help {
    local ret=$1
    shift

    echo "$0 creates a release of the SIESTA code at the tip-tag"
    echo ""
    echo "These options may be used to control the archive."
    echo ""
    echo "  --prev-tag instead of selecting the second-latest tag, choose this tag as the"
    echo "             reference tag for creating a diff with regards to the --tag tag [bzr tags --sort=time]"
    echo "  --tag      instead of selecting the latest tag, choose this tag as the"
    echo "             reference tag for creating a release archive [bzr tags --sort=time]"
    echo "  --out      the default output file is siesta-<tag>.tar.gz."
    echo "  --help|-h  show this help."
    exit $ret
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
	    if [ $(_has_tag $_tag) -eq 1 ]; then
		echo "$0 could not find tag: '$_tag' in the tags list."
		echo "The available tags are:"
		bzr tags | awk '{print " ",$1}'
		exit 1
	    fi
	    shift
	    ;;
	--prev-tag)
	    _prev_tag=$1
	    # Check that the tag exists
	    if [ $(_has_tag $_prev_tag) -eq 1 ]; then
		echo "$0 could not find tag: '$_prev_tag' in the tags list."
		echo "The available tags are:"
		bzr tags | awk '{print " ",$1}'
		exit 1
	    fi
	    shift
	    ;;
	
	--out)
	    _out=$1
	    shift
	    ;;
	
	--help|-h)
	    _help 0
	    ;;
	
    esac
done


# 
	    
# The current procedure of releasing a SIESTA
# version is the following:
# 0. Create a temporary directory
# 1. Make the proper documentation
# 2. Create diff files (both detailed and non-detailed)
# 3. Go out of the top-directory to create the repository.


#   Create a temporary work-directory
cd ../
# If the release already exists... Tell the user and quit
if [ -d $_reldir/$_tag ]; then
    echo "The release has already been processed."
    echo "Delete this folder:"
    echo "  rm -rf $(dirname $main_dir)/$_reldir/$_tag"
    exit 1
fi
bzr branch $main_dir -r $_tag release-manager-$_tag
mkdir -p $_reldir
# Move the branch into the directory.
mv release-manager-$_tag $_reldir/$_tag

# Go into the release directory where all work will be done
pushd $_reldir
pushd $_tag
rel_dir=`bzr root`


#   Create documentation
pushd Docs
#   First create the screen variants...
make final-screen
#   Then the regular documentation (for print)
make final
#   Only retain the pdf files
make clean

#   Go out of the documentation directory...
popd

# Create the two different CHANGES file
{
    echo "##############################################"
    echo " Changes between $_prev_tag and $_tag"
    echo "##############################################"
    echo ""
    bzr log -r$_prev_tag..$_tag
} > CHANGES
{
    echo "##############################################"
    echo " Detailed Changes between $_prev_tag and $_tag"
    echo "##############################################"
    echo ""
    bzr log -n 0 -r$_prev_tag..$_tag
} > CHANGES_DETAILED


# Return from the branch release
popd

# Tar the file
# In some future cases will the tag and out name
# be equivalent and we shouldn't delete what we should process!
if [ $_out != $_tag ]; then
    rm -rf $_out
    cp -rf $_tag $_out
fi

# Create the archive... (without any .bzr files)
tar --exclude '.bzr*' -czf $_out.tar.gz $_out

if [ $_out != $_tag ]; then
    rm -rf $_out
fi

popd
