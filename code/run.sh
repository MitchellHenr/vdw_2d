#! /usr/local/bin/bash
set -e

FILE="$( dirname "${BASH_SOURCE[0]}" )/vdw.py"

# ATOMS=( H He Na )
# MATERIALS=( MoS2 graphene )

# for atom in "${ATOMS[@]}"
# do
# 	for material in "${MATERIALS[@]}"
# 	do
# 		echo $material $atom
# 		python3 $FILE -g $atom $material
# 	done
# done
python3 $FILE -g H MoS2
