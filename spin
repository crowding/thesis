#!/bin/bash

list="git ls-files --cached"
build="make"
view="true"
error="true"

read -r -d '' usage <<EOF
$0: [-b builder] [-v viewer] [-l lister] [-v viewer] [...]
 -b builder
    A command to run to build a target whenever a change is detected.
    Default is "$build"
 -v viewer
    A command to view or test the results of compilation.
    Default is "$view"
 -e error
    A command to be run if the builder or viewer returns an error.
    Default is "$error"
 -l A command that lists the files to watch.
    Default is "$list"

Example usage: Remake the dissertation whenever I save an edit

$0 -b "make thesis.pdf" -v "open thesis.pdf -a Skim.app" -e "\$EDITOR thesis.log"
EOF

while getopts ":hb:v:e:l:" flag; do
    case $flag in
        h ) echo help; echo "$usage"; exit 0 ;;
        b ) build=$OPTARG ;;
        v ) view=$OPTARG ;;
        l ) list=$OPTARG ;;
        e ) error=$OPTARG ;;
        \? ) echo "Invalid option -$OPTARG" >& 2; exit 1 ;;
        : ) echo "$OPTARG requires an option" >& 2; exit 1 ;;
    esac
done

#kqwait comes from https://github.com/sschober/kqwait
while ( $build ) && ( $view ) || ( $error ); do
    kqwait $( $list )
done
