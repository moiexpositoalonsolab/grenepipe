#!/bin/bash

####################################################################################################
#    Helper functions
####################################################################################################

function print_separator() {
    echo
    echo -en "\e[32m==========================================================================="
    echo     "========================="
    echo -e  "          ${1}"
    echo -en "================================================================================="
    echo -e  "===================\e[0m"
    echo
}

function get_user_confirmation() {
    if [[ ${1} != "" ]]; then
        text=$1
    else
        text="Do you want to continue?"
    fi

    response="x"
    while [[ ${response} != "y" ]]; do
        read -p "${text} [y/n] " -n 1 -r response
        echo

        if [[ $response == "y" ]]; then
            return 1
        fi
        if [[ $response == "n" ]]; then
            return 0
        fi
    done
}

####################################################################################################
#    Check preconditions
####################################################################################################

print_separator "Check preconditions"

# Change to top level of git repo.
# This ensures that the script can be called from any directory.
cd `git rev-parse --show-toplevel`

# Check repo mame.
# If applied to a different repo, this script might do weird stuff, so better check.
base_name=`git rev-parse --show-toplevel | xargs basename`
if [[ ${base_name} != grenepipe* ]]; then
    echo -e "\e[31mWrong repository. Expect grenepipe.\e[0m"
    exit
else
    echo "Repository: ${base_name}"
fi

# Get and check current branch.
branch=`git rev-parse --abbrev-ref HEAD`
if [[ ${branch} != "master" ]]; then
    echo -e "\e[31mNot on master branch. Aborting.\e[0m"
    exit
else
    echo "Branch: ${branch}"
fi

# Check changed files.
changed_files=`git diff --name-only HEAD`
if [[ ${changed_files} != "" ]]; then
    echo -e "\e[31mUncommitted files. Aborting.\e[0m"
    exit
else
    echo "No uncommitted files."
fi

# Check stash.
stashed_files=`git stash list`
if [[ ${stashed_files} != "" ]]; then
    # echo -e "\e[31mStashed files. Aborting.\e[0m"
    # exit

    echo "There are staged files:"
    echo ${stashed_files}
    echo

    get_user_confirmation
    cont=$?
    if [[ $cont == 0 ]]; then
        echo -e "\e[31mAborted.\e[0m"
        exit
    fi
else
    echo "No stashed files."
fi

# Check for unmerged branches.
unmerged=`git branch --no-merged`
if [[ ${unmerged} != "" ]]; then
    echo "There are unmerged branches:"
    echo ${unmerged}
    echo

    get_user_confirmation
    cont=$?
    if [[ $cont == 0 ]]; then
        echo -e "\e[31mAborted.\e[0m"
        exit
    fi
else
    echo "No unmerged branches."
fi

####################################################################################################
#    Version
####################################################################################################

print_separator "Version"

# Output current version.
last_tag=`git describe --abbrev=0 --tags`
echo -en "Latest version tag:    \e[34m${last_tag}\e[0m\n\n"

# Read version tag.
read -p "Enter new version tag: v" version
vversion="v${version}"
# echo ${version}
echo

# Replace version line in common.smk file.
echo "Replace version in workflow/rules/initialize.smk"
sed -i.bak -e "s/grenepipe_version = \".*\" #GRENEPIPE_VERSION#/grenepipe_version = \"${version}\" #GRENEPIPE_VERSION#/g" workflow/rules/initialize.smk
rm workflow/rules/initialize.smk.bak

####################################################################################################
#    Commit and Tag
####################################################################################################

print_separator "Commit and Tag"

# Confirm changes.
echo -e "\e[34mCurrent git status:\e[0m\n"
git status
echo -e "\n\e[34mAbout to commit changes and to create version tag ${vversion}\e[0m\n"

get_user_confirmation
cont=$?
if [[ $cont == 0 ]]; then
    echo -e "\e[31mAborted.\e[0m"
    exit
fi
echo

# Commit and Tag.
echo "Make commit, create tag..."
git commit --allow-empty -am "Release ${vversion}"
git tag -a "${vversion}" -m "Release ${vversion}"
echo -e "...done."

####################################################################################################
#    Publish to GitHub
####################################################################################################

print_separator "Publish to GitHub"

# Show changes.
echo -e "\e[34mCurrent git status:\e[0m\n"
git status
echo

# Confirm.
get_user_confirmation "Do you want to push the commit and tag to GitHub?"
cont=$?

# Push.
if [[ $cont == 1 ]]; then
    echo
    git push --all --follow-tags
fi

####################################################################################################
#    Prepare GitHub Release
####################################################################################################

print_separator "Prepare GitHub Release"

# Get current (the new) tag.
new_tag=`git describe --abbrev=0 --tags`
if [[ ${new_tag} != ${vversion} ]]; then
    echo "Weird, version ${vversion} and tag ${new_tag} differ..."
    echo
fi

# Get all important commits after the last tag and format them for Markdown.
echo -e "\e[34m### Notable Changes ###\e[0m\n"
git log ${last_tag}..${new_tag} --oneline | cut -d " " -f 1 --complement | egrep -iv "^(Minor|Merge|Release)" | sed "s/^/  \* /g"
echo -e "\nUse this list for the release message.\n"

# Ask user for github page.
get_user_confirmation "Do you want to open the GitHub release page?"
cont=$?

# Open github release page.
if [[ $cont == 1 ]]; then
    github_release_url="https://github.com/moiexpositoalonsolab/grenepipe/releases/new"
    if which xdg-open > /dev/null
    then
        xdg-open ${github_release_url}
    elif which gnome-open > /dev/null
    then
        gnome-open ${github_release_url}
    else
        echo "Cannot open page."
    fi
    sleep 1
fi

####################################################################################################
#    Finished
####################################################################################################

print_separator "Finished"
