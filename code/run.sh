#! /usr/bin/env bash
set -e

FILE="$( dirname "${BASH_SOURCE[0]}" )/vdw.py"

ATOMS=( H He Na )
MATERIALS=( MoS2 graphene )

for atom in "${ATOMS[@]}"
do
	for material in "${MATERIALS[@]}"
	do
		echo $material $atom
		python3 $FILE -g $atom $material
	done
done

say "Your van der Waals forces are ready"
