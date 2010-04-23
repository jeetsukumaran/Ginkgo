#! /bin/sh

function echo_git_info {
    if [[ $(which git 2> /dev/null) ]]
    then
        local STATUS
        STATUS=$(git status 2>/dev/null)
        if [[ -z $STATUS ]]
        then
            return
        fi
        echo "`git symbolic-ref HEAD 2> /dev/null | cut -b 12-`-`git log --pretty=format:\"%h\" -1`"
    fi
}

echo "// Auto-generated build information."
echo
echo "// Date and time of build."
echo "#define BUILDTIMESTAMP \"`date`\""
echo
echo "// Branch and SHA-1 of current HEAD."
echo "#define BUILDDESC \"`echo_git_info`\""

