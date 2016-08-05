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
shared_dir=$(dirname $main_dir)
# Ensure we are at the root of this file
pushd $main_dir


# First we setup the default options
# _reldir is the main directory of the SIESTA shared
# repository (or at least it should be)
_reldir=$(dirname $main_dir)/siesta-releases
_sign=1

# Get the previous major release tag
_prev_tag=`bzr tags --sort=time | grep -e '^v' | grep -v '-' | tail -2 | head -1 | awk '{print $1}'`
# Get the latest tag
_tag=`bzr tags --sort=time | grep -e '^v' | tail -1 | awk '{print $1}'`

# Get default output file (siesta-<>.tar.gz)
_out=siesta-${_tag//v/}
_out=${_out//-release/}
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


# Function for signing a specific file
function _sign {
    local file=$1
    if [ -e $file ]; then
	if [ $_sign -eq 1 ]; then
	    gpg --armor --sign --detach-sig $file
	fi
	_store $file
	if [ $_sign -eq 1 ]; then
	    _store $file.asc
	    rm -f $file.asc
	fi
    fi
}

# Save file in the $_tag-files folder
function _store {
    local file=$1
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

	--no-sign)
	    _sign=0
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


echo "Chosen tags are:"
echo " previous tag: $_prev_tag"
echo " current tag : $_tag"
sleep 1
	    
# The current procedure of releasing a SIESTA
# version is the following:
# 0. Create the necessary directories
# 1. Create CHANGES and sign-store them.
# 2. Make the proper documentation
# 3. Go out of the top-directory to create the repository.

# Top release directories, in case it is not already
# created.
mkdir -p $_reldir
# The final directory with ALL release files, possibly
# signed.
mkdir -p $_reldir/$_tag-files

# If the release already exists... Tell the user and quit
if [ -d $_reldir/$_out ]; then
    echo "The release has already been processed."
    echo "Delete this folder:"
    echo "  rm -rf $(dirname $main_dir)/$_reldir/$_tag"
    exit 1
fi


#   Create a temporary work-directory
bzr export -r $_tag $_reldir/$_out $main_dir

# Create the changes files
# This is necessary to do here as the previous
# siesta versions may not have the tags related.
{
    echo "##############################################"
    echo " Changes between $_prev_tag and $_tag"
    echo "##############################################"
    echo ""
    bzr log -r$_prev_tag..$_tag
} > $_reldir/$_out/CHANGES
{
    echo "##############################################"
    echo " Detailed Changes between $_prev_tag and $_tag"
    echo "##############################################"
    echo ""
    bzr log -n 0 -r$_prev_tag..$_tag
} > $_reldir/$_out/CHANGES_DETAILED


# Go into the release directory where all work will be done
pushd $_reldir

pushd $_out


#   Create documentation
pushd Docs
#   First create the screen variants...
make final-screen
#   Then the regular documentation (for print)
make final

#   Clean-up the non-pdf files (currently
#   this is not available through the make clean)
rm -f siesta*.[^pt]* siesta*.toc
rm -f tbtrans*.[^pt]* tbtrans*.toc

# Create signatures and move files
for f in *.pdf ; do
    _sign $f
done

#   Go out of the documentation directory...
popd



# Return from the branch release
popd

# Tar the file
# In some future cases will the tag and out name
# be equivalent and we shouldn't delete what we should process!

# Create the archive... (without any .bzr files)
tar -czf $_out.tar.gz $_out
_sign $_out.tar.gz
rm $_out.tar.gz

popd
