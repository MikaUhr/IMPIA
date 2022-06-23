#!/bin/sh

# multijoin - join multiple files

out=$1; shift
tmp=$2; shift

join_rec() {
    mv ${out} ${tmp}
    if [ $# -eq 1 ]; then
        python join_tpm.py ${tmp} "$1" "${out}"
    else
        echo '$#' $#
        f=$1
        echo '$1:' $f
        shift
        python join_tpm.py ${tmp} "$f" ${out} | join_rec "$@"
    fi
}

if [ $# -le 2 ]; then
    #echo '$#' $#
    python join_tpm.py "$@" "${out}"
else
    f1=$1; f2=$2;
    echo "f1:" $f1
    echo "f2:" $f2
    shift 2
    python join_tpm.py "$f1" "$f2" "${out}" | join_rec "$@"
fi
